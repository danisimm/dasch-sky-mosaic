"""Direct S3 download helpers for plate photos and FITS mosaics.

DASCH files live in the `dasch-prod-user` bucket. This module provides
functions to download them directly from S3, bypassing the Starglass API.
This is cheaper (fewer API calls) and useful for batch processing.

Known S3 key structure:
  Photo (JPG):  plates/{plate_id}/{plate_id}_pphoto_all.jpg
  FITS mosaic:  plates/{plate_id}/{plate_id}_mosaic_{binning}.fits  (TBD - verify key)
"""
from __future__ import annotations

import logging
from pathlib import Path

import requests

LOG = logging.getLogger(__name__)

S3_BUCKET = "dasch-prod-user"
S3_REGION = "us-east-1"
S3_BASE_URL = f"https://{S3_BUCKET}.s3.{S3_REGION}.amazonaws.com"


def _s3_photo_key(plate_id: str) -> str:
    return f"plates/{plate_id}/{plate_id}_pphoto_all.jpg"


def _s3_photo_url(plate_id: str) -> str:
    return f"{S3_BASE_URL}/{_s3_photo_key(plate_id)}"


def _download_stream(url: str, dest: Path, overwrite: bool = True) -> int:
    """Stream a URL to disk. Returns file size in bytes."""
    if dest.exists() and not overwrite:
        return dest.stat().st_size
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=180.0) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_content(chunk_size=65536):
                if chunk:
                    f.write(chunk)
    return dest.stat().st_size


def download_photo_from_s3(
    plate_id: str,
    dest_dir: Path,
    overwrite: bool = False,
) -> Path:
    """Download a plate JPG directly from S3.

    The file is cached at dest_dir/{plate_id}_pphoto_all.jpg.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{plate_id}_pphoto_all.jpg"

    if dest.exists() and not overwrite:
        LOG.info("Reusing cached photo for %s: %s", plate_id, dest.name)
        return dest

    url = _s3_photo_url(plate_id)
    LOG.info("Downloading photo for %s from S3", plate_id)
    _download_stream(url, dest, overwrite=overwrite)
    return dest


def download_fits_from_s3(
    plate_id: str,
    dest_dir: Path,
    binning: int = 16,
    overwrite: bool = False,
) -> Path:
    """Download a FITS mosaic directly from S3.

    NOTE: The S3 key structure for FITS mosaics has not been confirmed.
    This function is a stub. Update the key template once the correct path
    is known.
    """
    # TODO: confirm the actual S3 key structure for FITS mosaics
    raise NotImplementedError(
        "FITS S3 key structure not yet confirmed. "
        "Use call_sg.download_fits_via_sg() for now."
    )
