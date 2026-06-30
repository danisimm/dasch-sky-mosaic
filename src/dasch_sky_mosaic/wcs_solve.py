"""Solve WCS for plate JPGs using astroalign.

Given a plate photo (JPG) and a FITS mosaic, aligns the photo to the FITS
WCS via star pattern matching and returns a WCS header for the photo.

Import solve_plate_wcs() for library use. Downloads and WTML generation
are handled by pipeline.py.
"""
from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.wcs import FITSFixedWarning, WCS

LOG = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# FITS loading
# ---------------------------------------------------------------------------

def _load_plate_fits(fits_path: Path) -> tuple[dict[str, Any], WCS, tuple[int, int], Any]:
    """Return (wcs_header, wcs, shape, image_data) from a FITS mosaic."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FITSFixedWarning)
        with fits.open(fits_path, memmap=True) as hdul:
            for ext in hdul:
                raw = getattr(ext, "data", None)
                if raw is None or not hasattr(raw, "ndim") or raw.ndim != 2:
                    continue
                header = ext.header.copy()
                shape = (int(raw.shape[0]), int(raw.shape[1]))
                data_copy = raw.copy()
                wcs = WCS(header)
                wcs_header = dict(wcs.to_header(relax=True))
                return wcs_header, wcs, shape, data_copy
    raise RuntimeError(f"no 2D image HDU found in {fits_path}")


# ---------------------------------------------------------------------------
# Astroalign backend
# ---------------------------------------------------------------------------

def _normalize_image(img: Any) -> Any:
    p1, p99 = np.percentile(img, (1.0, 99.0))
    if p99 <= p1:
        return img.astype(np.float32) * 0
    img = img.clip(p1, p99)
    return ((img - p1) / (p99 - p1)).astype(np.float32)


def _astroalign_find_transform(
    alignment_source: Any,
    fits_data: Any,
) -> tuple[Any | None, int]:
    """Find a similarity transform mapping plate-photo pixels to FITS pixels.

    Returns (transform, n_matches). Photos have dark stars so the source is
    inverted before detection.
    """
    try:
        import astroalign as aa
    except ImportError:
        LOG.warning("astroalign not installed; cannot solve WCS via this backend.")
        return None, 0

    try:
        # Photos have dark stars; invert so peaks are detectable by astroalign
        src_gray = 1.0 - _normalize_image(np.array(alignment_source, dtype=np.float32))
        fits_norm = _normalize_image(np.array(fits_data, dtype=np.float32))

        transform, (src_pts, _) = aa.find_transform(
            src_gray,
            fits_norm,
            max_control_points=150,
            detection_sigma=4,
            min_area=3,
        )
        n = int(len(src_pts))
        LOG.info(
            "Astroalign: matches=%d scale=%.5f rot_deg=%.3f",
            n,
            float(transform.scale),
            float(np.degrees(transform.rotation)),
        )
        return transform, n
    except Exception:
        LOG.exception("Astroalign failed.")
        return None, 0


def _compose_wcs_from_affine(
    reference_wcs: WCS,
    affine_a: Any,
    affine_t: Any,
) -> dict[str, Any]:
    """Compose a photo→FITS affine transform into a FITS WCS header dict."""
    if reference_wcs.wcs.has_cd():
        ref_cd = np.array(reference_wcs.wcs.cd, dtype=np.float64)
    else:
        ref_cd = (
            np.array(reference_wcs.wcs.get_pc(), dtype=np.float64)
            @ np.diag(np.array(reference_wcs.wcs.cdelt, dtype=np.float64))
        )

    ones = np.ones(2, dtype=np.float64)
    affine_a = np.array(affine_a, dtype=np.float64)
    affine_t = np.array(affine_t, dtype=np.float64)
    cd_src = ref_cd @ affine_a
    b = affine_t + ones - affine_a @ ones
    crpix_src = np.linalg.solve(affine_a, np.array(reference_wcs.wcs.crpix, dtype=np.float64) - b)

    header = reference_wcs.to_header(relax=True)
    for key in list(header.keys()):
        if key.startswith("PC") or key.startswith("CD"):
            del header[key]
    for key in ("CDELT1", "CDELT2"):
        if key in header:
            del header[key]

    header["CRPIX1"] = float(crpix_src[0])
    header["CRPIX2"] = float(crpix_src[1])
    header["CD1_1"] = float(cd_src[0, 0])
    header["CD1_2"] = float(cd_src[0, 1])
    header["CD2_1"] = float(cd_src[1, 0])
    header["CD2_2"] = float(cd_src[1, 1])
    return dict(header)


def _solve_wcs_astroalign(
    photo_path: Path,
    fits_path: Path,
) -> tuple[dict[str, Any], dict[str, Any], int, int]:
    """Align photo JPG to FITS WCS using astroalign.

    Returns (wcs_header, alignment_meta, image_width, image_height).
    """
    from PIL import Image

    wcs_header_fits, reference_wcs, shape, fits_data = _load_plate_fits(fits_path)
    target_h, target_w = shape

    with Image.open(photo_path) as opened:
        src = opened.convert("L")
    src_w, src_h = src.size

    def _cover_resize(candidate: Any) -> tuple[Any, float]:
        c_w, c_h = candidate.size
        scale = max(target_w / max(1, c_w), target_h / max(1, c_h))
        resized = candidate.resize(
            (max(1, int(round(c_w * scale))), max(1, int(round(c_h * scale)))),
            resample=Image.Resampling.BICUBIC,
        )
        return resized, float(scale)

    parity_candidates = [
        ("source", src),
        ("source_mirror", src.transpose(Image.Transpose.FLIP_LEFT_RIGHT)),
    ]
    best_header: dict[str, Any] | None = None
    best_meta: dict[str, Any] | None = None
    best_matches = -1

    for label, candidate in parity_candidates:
        resized, resize_scale = _cover_resize(candidate)
        x0 = max(0, min((resized.width - target_w) // 2, resized.width - target_w))
        y0 = max(0, min((resized.height - target_h) // 2, resized.height - target_h))
        match_img = resized.crop((x0, y0, x0 + target_w, y0 + target_h))

        transform, n_matches = _astroalign_find_transform(match_img, fits_data)
        if transform is None:
            continue

        if label == "source_mirror":
            base_a = np.array([[-1.0, 0.0], [0.0, 1.0]], dtype=np.float64) * resize_scale
            base_t = np.array([src_w - 1.0, 0.0], dtype=np.float64) * resize_scale
        else:
            base_a = np.eye(2, dtype=np.float64) * resize_scale
            base_t = np.zeros(2, dtype=np.float64)
        base_t -= np.array([x0, y0], dtype=np.float64)

        total_a = np.array(transform.params[:2, :2], dtype=np.float64) @ base_a
        total_t = (
            np.array(transform.params[:2, 2], dtype=np.float64)
            + np.array(transform.params[:2, :2], dtype=np.float64) @ base_t
        )
        header = _compose_wcs_from_affine(reference_wcs, total_a, total_t)

        if n_matches > best_matches:
            best_header = header
            best_matches = n_matches
            best_meta = {
                "backend": "astroalign",
                "parity": label,
                "matches": n_matches,
                "transform_scale": float(transform.scale),
                "transform_rotation_deg": float(np.degrees(transform.rotation)),
            }

    if best_header is None:
        raise RuntimeError(f"Astroalign found no match for {photo_path.name}")

    LOG.info("Astroalign parity pick: %s (matches=%d)", best_meta["parity"], best_matches)
    return best_header, best_meta, src_w, src_h


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def solve_plate_wcs(
    plate_id: str,
    photo_path: Path,
    fits_path: Path,
) -> dict[str, Any]:
    """Solve WCS for a plate JPG and return WCS data.

    Returns a dict with:
      - plate_id
      - wcs_header: dict of FITS WCS keywords for the JPG
      - image_width_px, image_height_px
      - alignment_meta: solver diagnostics
    """
    wcs_header, alignment_meta, image_w, image_h = _solve_wcs_astroalign(photo_path, fits_path)
    return {
        "plate_id": plate_id,
        "wcs_header": wcs_header,
        "image_width_px": image_w,
        "image_height_px": image_h,
        "alignment_meta": alignment_meta,
    }
