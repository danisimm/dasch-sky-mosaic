"""Orchestrates downloading and WCS solving for DASCH plates.

Run directly:
    python -m dasch_sky_mosaic.pipeline <plate_id> [plate_id ...]

Edit the configuration block at the bottom of this file before running.
WCS results are saved as JSON to WCS_CACHE_DIR. WTML generation is handled
separately (e.g. by a Lambda function reading those JSON files from S3).
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from dasch_sky_mosaic.call_s3 import download_fits_from_s3, download_photo_from_s3
from dasch_sky_mosaic.call_sg import download_fits_via_sg, download_photo_via_sg
from dasch_sky_mosaic.wcs_solve import solve_plate_wcs

LOG = logging.getLogger(__name__)


def run_pipeline(
    plate_ids: list[str],
    photo_dir: Path = Path("data/cache/photos"),
    fits_dir: Path = Path("data/cache/mosaic_package"),
    download_method: str = "s3",
    wcs_cache_dir: Path | None = None,
) -> list[dict[str, Any]]:
    """Download plate files and solve WCS for each plate.

    Returns a list of records, each containing:
      - plate_id
      - wcs_header: dict of FITS WCS keywords for the JPG
      - image_width_px, image_height_px
      - alignment_meta: solver diagnostics

    Args:
        plate_ids: Plate IDs to process.
        photo_dir: Directory to cache downloaded JPGs.
        fits_dir: Directory to cache downloaded FITS.
        download_method: 's3' (direct, no rate limits) or 'sg' (Starglass API).
        wcs_cache_dir: If set, saves a WCS JSON file per plate here.
    """
    records: list[dict[str, Any]] = []

    for plate_id in plate_ids:
        plate_id = plate_id.strip().lower()
        LOG.info("Processing %s", plate_id)

        try:
            if download_method == "s3":
                photo_path = download_photo_from_s3(plate_id, photo_dir)
                fits_path = download_fits_from_s3(plate_id, fits_dir)
            elif download_method == "sg":
                photo_path = download_photo_via_sg(plate_id, photo_dir)
                fits_path = download_fits_via_sg(plate_id, fits_dir)
            else:
                raise ValueError(f"Unknown download_method: {download_method!r}. Choose 's3' or 'sg'.")

            result = solve_plate_wcs(plate_id=plate_id, photo_path=photo_path, fits_path=fits_path)

            if wcs_cache_dir is not None:
                wcs_cache_dir.mkdir(parents=True, exist_ok=True)
                (wcs_cache_dir / f"{plate_id}_wcs.json").write_text(
                    json.dumps(result, indent=2), encoding="utf-8"
                )

            records.append(result)
            LOG.info("Solved %s", plate_id)

        except Exception as exc:
            LOG.error("Failed for %s: %s", plate_id, exc)

    return records


if __name__ == "__main__":
    import argparse

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    # ── Configuration ────────────────────────────────────────────────────────
    DOWNLOAD_METHOD = "s3"                        # "s3" or "sg"
    PHOTO_DIR     = Path("data/cache/photos")
    FITS_DIR      = Path("data/cache/mosaic_package")
    WCS_CACHE_DIR = Path("data/output/wcs")       # set to None to skip
    # ─────────────────────────────────────────────────────────────────────────

    parser = argparse.ArgumentParser(description="Download and WCS-solve DASCH plates.")
    parser.add_argument("plate_ids", nargs="+", help="Plate IDs to process (e.g. b02312 a11740)")
    args = parser.parse_args()

    records = run_pipeline(
        plate_ids=args.plate_ids,
        photo_dir=PHOTO_DIR,
        fits_dir=FITS_DIR,
        download_method=DOWNLOAD_METHOD,
        wcs_cache_dir=WCS_CACHE_DIR,
    )
    print(f"Done. {len(records)} plates solved.")
