"""Direct S3 download helpers for plate photos and FITS mosaics.

DASCH files live in the `dasch-prod-user` bucket. This module provides
functions to download them directly from S3, bypassing the Starglass API.
This is cheaper (fewer API calls) and useful for batch processing.

Known S3 key structure (from mosaic_package API response):
  Photo (JPG):  plates/{plate_id}/{plate_id}_pphoto_all.jpg
  FITS mosaic:  plates/{plate_id}/{plate_id}_mosaic_{mosnum:02d}_16ww.fit.fz

  Always uses binning=16 (the binning=1 files are too large).
  Use download_fits_via_sg() in call_sg if you need the exact presigned URL from the API.
"""
from __future__ import annotations

import logging
from pathlib import Path

from dasch_sky_mosaic.utils import _download_stream

LOG = logging.getLogger(__name__)

S3_BUCKET = "dasch-prod-user"
S3_REGION = "us-east-1"
S3_BASE_URL = f"https://{S3_BUCKET}.s3.{S3_REGION}.amazonaws.com"


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

    url = f"{S3_BASE_URL}/plates/{plate_id}/{plate_id}_pphoto_all.jpg"
    LOG.info("Downloading photo for %s from S3", plate_id)
    _download_stream(url, dest)
    return dest


def download_fits_from_s3(
    plate_id: str,
    dest_dir: Path,
    overwrite: bool = False,
) -> Path:
    """Download a binning=16 FITS mosaic directly from S3 and decompress to .fits.

    Tries mosnum=0 then mosnum=1 since the correct value varies per plate.
    Use download_fits_via_sg() if you need the guaranteed-correct key from the API.
    """
    from dasch_sky_mosaic.call_sg import _decompress_fz, _unpacked_path

    dest_dir.mkdir(parents=True, exist_ok=True)

    for mosnum in (0, 1):
        filename = f"{plate_id}_mosaic_{mosnum:02d}_16ww.fit.fz"
        dest_fz = dest_dir / filename
        dest_fits = _unpacked_path(dest_fz)
        if dest_fits.exists() and not overwrite:
            LOG.info("Reusing cached FITS for %s: %s", plate_id, dest_fits.name)
            return dest_fits
        url = f"{S3_BASE_URL}/plates/{plate_id}/{filename}"
        try:
            LOG.info("Downloading FITS for %s from S3 (mosnum=%d)", plate_id, mosnum)
            _download_stream(url, dest_fz)
            return _decompress_fz(dest_fz)
        except Exception:
            if mosnum == 1:
                raise
            LOG.debug("mosnum=0 not found for %s, trying mosnum=1", plate_id)
