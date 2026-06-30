"""Solve WCS for plate JPGs using astroalign.

Given a plate photo (JPG) and a FITS mosaic, aligns the photo to the FITS
WCS via star pattern matching and returns a WCS header for the photo.

Import solve_plate_wcs() for library use. Downloads and WTML generation
are handled by pipeline.py.
"""
from __future__ import annotations

import logging
import threading
import warnings
from pathlib import Path
from typing import Any

import astroalign as aa
import numpy as np
from astropy.io import fits
from astropy.wcs import FITSFixedWarning, WCS
from PIL import Image

LOG = logging.getLogger(__name__)


def _load_plate_fits(fits_path: Path) -> tuple[WCS, tuple[int, int], Any]:
    """Return (wcs, shape, image_data) from a FITS mosaic.

    DASCH mosaics always have a single PrimaryHDU at index 0.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FITSFixedWarning)
        with fits.open(fits_path, memmap=True) as hdul:
            hdu = hdul[0]
            raw = getattr(hdu, "data", None)
            if raw is None:
                raise RuntimeError(f"PRIMARY HDU has no data in {fits_path}")
            header = getattr(hdu, "header").copy()
            shape = (int(raw.shape[0]), int(raw.shape[1]))
            data = raw.copy()
    wcs = WCS(header)
    return wcs, shape, data


def _normalize_image(img: Any) -> Any:
    p1, p99 = np.percentile(img, (1.0, 99.0))
    if p99 <= p1:
        return img.astype(np.float32) * 0
    img = img.clip(p1, p99)
    return ((img - p1) / (p99 - p1)).astype(np.float32)


_ASTROALIGN_TIMEOUT_S = 5.0


def _astroalign_find_transform(
    plate_img: Any,
    fits_norm: Any,
) -> tuple[Any | None, int]:
    """Find a similarity transform mapping plate-photo pixels to FITS pixels.

    Returns (transform, n_matches). Photos have dark stars so the source is
    inverted before detection. fits_norm should be pre-normalized.

    Gives up after _ASTROALIGN_TIMEOUT_S seconds so a wrong-parity attempt
    doesn't block for minutes exhausting all triangle combinations.
    """
    src_gray = 1.0 - _normalize_image(np.array(plate_img, dtype=np.float32))

    _result: list[Any] = [None]
    _exc: list[BaseException | None] = [None]

    def _run() -> None:
        try:
            transform, (src_pts, _) = aa.find_transform(
                src_gray,
                fits_norm,
                max_control_points=100,
                detection_sigma=5,
                min_area=5,
            )
            _result[0] = (transform, src_pts)
        except Exception as exc:
            _exc[0] = exc

    thread = threading.Thread(target=_run, daemon=True)
    thread.start()
    thread.join(timeout=_ASTROALIGN_TIMEOUT_S)

    if thread.is_alive():
        LOG.info("Astroalign timed out after %.1fs", _ASTROALIGN_TIMEOUT_S)
        return None, 0

    if _exc[0] is not None:
        LOG.info("Astroalign failed: %s", _exc[0])
        return None, 0

    transform, src_pts = _result[0]
    n = int(len(src_pts))
    LOG.info(
        "Astroalign: matches=%d scale=%.5f rot_deg=%.3f",
        n,
        float(transform.scale),
        float(np.degrees(transform.rotation)),
    )
    return transform, n


def _cover_resize(
    candidate: Any,
    target_w: int,
    target_h: int,
) -> tuple[Any, float]:
    """Scale candidate image to cover (target_w, target_h), return (resized, scale)."""
    c_w, c_h = candidate.size
    scale = max(target_w / max(1, c_w), target_h / max(1, c_h))
    resized = candidate.resize(
        (max(1, int(round(c_w * scale))), max(1, int(round(c_h * scale)))),
        resample=Image.Resampling.BICUBIC,
    )
    return resized, float(scale)


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


def solve_plate_wcs(
    plate_id: str,
    photo_path: Path,
    fits_path: Path,
) -> dict[str, Any]:
    """Align a plate JPG to a FITS mosaic WCS and return the result.

    Returns a dict with:
      - plate_id
      - wcs_header: dict of FITS WCS keywords for the JPG
      - image_width_px, image_height_px
      - alignment_meta: solver diagnostics
    """
    reference_wcs, shape, fits_data = _load_plate_fits(fits_path)
    target_h, target_w = shape

    fits_norm = _normalize_image(np.array(fits_data, dtype=np.float32))

    with Image.open(photo_path) as opened:
        src = opened.convert("L")
    src_w, src_h = src.size

    resized, resize_scale = _cover_resize(src, target_w, target_h)
    x0 = max(0, min((resized.width - target_w) // 2, resized.width - target_w))
    y0 = max(0, min((resized.height - target_h) // 2, resized.height - target_h))
    match_img = resized.crop((x0, y0, x0 + target_w, y0 + target_h))

    parity_candidates = [
        ("source_mirror", match_img.transpose(Image.Transpose.FLIP_LEFT_RIGHT)),
        ("source", match_img),
    ]
    best_header: dict[str, Any] | None = None
    best_meta: dict[str, Any] | None = None

    for label, candidate in parity_candidates:
        transform, n_matches = _astroalign_find_transform(candidate, fits_norm)
        if transform is None:
            continue

        if label == "source_mirror":
            base_a = np.array([[-1.0, 0.0], [0.0, 1.0]], dtype=np.float64) * resize_scale
            base_t = np.array([target_w - 1.0 + x0, -y0], dtype=np.float64)
        else:
            base_a = np.eye(2, dtype=np.float64) * resize_scale
            base_t = np.array([-x0, -y0], dtype=np.float64)

        total_a = np.array(transform.params[:2, :2], dtype=np.float64) @ base_a
        total_t = (
            np.array(transform.params[:2, 2], dtype=np.float64)
            + np.array(transform.params[:2, :2], dtype=np.float64) @ base_t
        )
        best_header = _compose_wcs_from_affine(reference_wcs, total_a, total_t)
        best_meta = {
            "backend": "astroalign",
            "parity": label,
            "matches": n_matches,
            "transform_scale": float(transform.scale),
            "transform_rotation_deg": float(np.degrees(transform.rotation)),
        }
        break  # stop after first successful parity (mirror tried first)

    if best_header is None or best_meta is None:
        raise RuntimeError(f"Astroalign found no match for {photo_path.name}")

    LOG.info("Astroalign parity pick: %s (matches=%d)", best_meta["parity"], best_meta["matches"])

    return {
        "plate_id": plate_id,
        "wcs_header": best_header,
        "image_width_px": src_w,
        "image_height_px": src_h,
        "alignment_meta": best_meta,
    }
