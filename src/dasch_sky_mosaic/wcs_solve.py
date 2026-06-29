"""Solve WCS for plate JPGs and produce placement data for WTML generation.

Two solver backends are available (selected via --solver):
  astroalign  - Aligns the JPG to the FITS mosaic WCS using star pattern matching.
                Requires both a JPG and a FITS download per plate (~2 files).
  scamp       - [STUB] Source extraction + SCAMP astrometric solve.
                Requires only a JPG (~1 file). Not yet implemented.

For each plate the output is a JSON record containing the WCS header and the
WWT SkyImage placement fields needed to build a WTML link.

Run directly:
    python -m dasch_sky_mosaic.wcs_solve <plate_id> [plate_id ...] [options]

Or import solve_plate_wcs() for library use.
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.wcs import FITSFixedWarning, WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from dasch_sky_mosaic.call_s3 import download_photo_from_s3
from dasch_sky_mosaic.call_sg import download_fits_via_sg, download_photo_via_sg
from dasch_sky_mosaic.wtml_gen import derive_skyimage_placement

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
    import numpy as np
    p1, p99 = np.percentile(img, (1.0, 99.0))
    if p99 <= p1:
        return img.astype(np.float32) * 0
    img = img.clip(p1, p99)
    return ((img - p1) / (p99 - p1)).astype(np.float32)


def _astroalign_find_transform(
    alignment_source: Any,
    fits_data: Any,
) -> tuple[Any | None, int, str]:
    """Find a similarity transform mapping plate-photo pixels to FITS pixels.

    Returns (transform, n_matches, mode_label).
    """
    try:
        import astroalign as aa
    except ImportError:
        LOG.warning("astroalign not installed; cannot solve WCS via this backend.")
        return None, 0, "unavailable"

    try:
        src_rgb = np.array(alignment_source, dtype=np.float32)
        src_gray = _normalize_image(src_rgb[:, :, 0])
        fits_norm = _normalize_image(np.array(fits_data, dtype=np.float32))

        best_transform: Any = None
        best_label = ""
        best_matches = 0

        for label, source in (("native", src_gray), ("inverted", 1.0 - src_gray)):
            try:
                transform, (src_pts, _) = aa.find_transform(
                    source,
                    fits_norm,
                    max_control_points=150,
                    detection_sigma=4,
                    min_area=3,
                )
            except Exception:
                continue
            n = int(len(src_pts))
            if n > best_matches:
                best_matches = n
                best_transform = transform
                best_label = label

        if best_transform is None or best_matches < 6:
            return None, 0, "none"

        LOG.info(
            "Astroalign: mode=%s matches=%d scale=%.5f rot_deg=%.3f",
            best_label, best_matches,
            float(best_transform.scale),
            float(np.degrees(best_transform.rotation)),
        )
        return best_transform, best_matches, best_label
    except Exception:
        LOG.exception("Astroalign failed.")
        return None, 0, "error"


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
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Align photo JPG to FITS WCS using astroalign.

    Returns (wcs_header, alignment_meta).
    wcs_header is a dict of FITS-compatible WCS keywords for the JPG.
    """
    from PIL import Image

    wcs_header_fits, reference_wcs, shape, fits_data = _load_plate_fits(fits_path)
    target_h, target_w = shape

    with Image.open(photo_path) as opened:
        src = opened.convert("RGB")
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

        transform, n_matches, mode_label = _astroalign_find_transform(match_img, fits_data)
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
                "mode": mode_label,
                "matches": n_matches,
                "transform_scale": float(transform.scale),
                "transform_rotation_deg": float(np.degrees(transform.rotation)),
            }

    if best_header is not None and best_meta is not None:
        LOG.info("Astroalign parity pick: %s (matches=%d)", best_meta["parity"], best_matches)
        return best_header, best_meta

    # Fallback: simple centered scaling, no rotation.
    LOG.warning("Astroalign found no match for %s; using fallback centered scaling.", photo_path.name)
    resized, resize_scale = _cover_resize(src)
    x0 = max(0, min((resized.width - target_w) // 2, resized.width - target_w))
    y0 = max(0, min((resized.height - target_h) // 2, resized.height - target_h))
    fallback_a = np.eye(2, dtype=np.float64) * resize_scale
    fallback_t = np.array([-x0, -y0], dtype=np.float64)
    fallback_header = _compose_wcs_from_affine(reference_wcs, fallback_a, fallback_t)
    fallback_meta = {
        "backend": "astroalign",
        "parity": "source",
        "mode": "fallback_centered",
        "matches": 0,
    }
    return fallback_header, fallback_meta


def _solve_wcs_scamp(
    photo_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """[STUB] Solve WCS using SExtractor + SCAMP.

    This backend requires only the JPG (no FITS download) and may be faster
    at scale. Implement once the SCAMP output format is confirmed.
    """
    raise NotImplementedError(
        "SCAMP backend not yet implemented. Use --solver astroalign for now."
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def solve_plate_wcs(
    plate_id: str,
    photo_path: Path,
    fits_path: Path | None,
    solver: str = "astroalign",
) -> dict[str, Any]:
    """Solve WCS for a plate JPG and return placement data.

    Returns a dict with:
      - plate_id
      - wcs_header: dict of FITS WCS keywords for the JPG
      - placement: WWT SkyImage fields (from derive_skyimage_placement)
      - image_width_px, image_height_px
      - alignment_meta: solver-specific diagnostics
    """
    from PIL import Image

    with Image.open(photo_path) as img:
        image_w, image_h = img.size

    if solver == "astroalign":
        if fits_path is None:
            raise ValueError("astroalign solver requires a FITS file path")
        wcs_header, alignment_meta = _solve_wcs_astroalign(photo_path, fits_path)
    elif solver == "scamp":
        wcs_header, alignment_meta = _solve_wcs_scamp(photo_path)
    else:
        raise ValueError(f"Unknown solver: {solver!r}. Choose 'astroalign' or 'scamp'.")

    placement = derive_skyimage_placement(wcs_header, image_w, image_h)

    return {
        "plate_id": plate_id,
        "wcs_header": wcs_header,
        "placement": placement,
        "image_width_px": image_w,
        "image_height_px": image_h,
        "alignment_meta": alignment_meta,
    }


# ---------------------------------------------------------------------------
# Batch runner
# ---------------------------------------------------------------------------

def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Solve WCS for one or more DASCH plate JPGs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("plate_ids", nargs="+", help="Plate IDs to process (e.g. a11740 br00130)")
    p.add_argument(
        "--solver",
        choices=["astroalign", "scamp"],
        default="astroalign",
        help="WCS solver backend to use",
    )
    p.add_argument(
        "--download-method",
        choices=["s3", "sg"],
        default="s3",
        help="How to download plate JPGs: s3 (direct, cheaper) or sg (Starglass API)",
    )
    p.add_argument("--photo-dir", default="data/cache/photos", help="Directory to cache downloaded JPGs")
    p.add_argument("--fits-dir", default="data/cache/mosaic_package", help="Directory to cache downloaded FITS")
    p.add_argument("--output-dir", default="data/output/wcs", help="Directory for output JSON files")
    p.add_argument("--binning", type=int, choices=[1, 16], default=16, help="FITS mosaic binning level")
    p.add_argument("--api-base", default="auto", help="Starglass API base URL (for --download-method sg)")
    p.add_argument("--api-key", default=None, help="Starglass API key (or set STARGLASS_API_KEY env var)")
    p.add_argument("--overwrite", action="store_true", help="Re-download and re-solve even if cached")
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    from dotenv import load_dotenv
    load_dotenv()

    args = _parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s %(message)s",
        stream=sys.stderr,
    )

    api_key = args.api_key or os.environ.get("STARGLASS_API_KEY")
    photo_dir = Path(args.photo_dir)
    fits_dir = Path(args.fits_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results: list[dict[str, Any]] = []

    for plate_id in args.plate_ids:
        plate_id = plate_id.strip().lower()
        LOG.info("Processing plate %s", plate_id)
        out_path = output_dir / f"{plate_id}_wcs.json"

        if out_path.exists() and not args.overwrite:
            LOG.info("Skipping %s (cached at %s)", plate_id, out_path)
            results.append(json.loads(out_path.read_text()))
            continue

        try:
            # Download JPG
            if args.download_method == "s3":
                photo_path = download_photo_from_s3(plate_id, photo_dir, overwrite=args.overwrite)
            else:
                photo_path = download_photo_via_sg(
                    plate_id, photo_dir,
                    api_base=args.api_base, api_key=api_key,
                    overwrite=args.overwrite,
                )

            # Download FITS if needed
            fits_path: Path | None = None
            if args.solver == "astroalign":
                fits_path = download_fits_via_sg(
                    plate_id, fits_dir,
                    binning=args.binning,
                    api_base=args.api_base,
                    api_key=api_key,
                )

            result = solve_plate_wcs(
                plate_id=plate_id,
                photo_path=photo_path,
                fits_path=fits_path,
                solver=args.solver,
            )
            out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
            results.append(result)
            LOG.info("Saved WCS for %s → %s", plate_id, out_path)

        except Exception as exc:
            LOG.error("Failed for %s: %s", plate_id, exc)
            results.append({"plate_id": plate_id, "error": str(exc)})

    # Print a short summary
    n_ok = sum(1 for r in results if "error" not in r)
    n_fail = len(results) - n_ok
    print(f"Done: {n_ok}/{len(results)} plates solved, {n_fail} failed.")


if __name__ == "__main__":
    main()
