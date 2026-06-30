"""Direct S3 download helpers for plate photos and FITS mosaics.

DASCH files live in the `dasch-prod-user` bucket. This module provides
functions to download them directly from S3, bypassing the Starglass API.
Advantages over call_sg: no Starglass rate limiting, and no reliance on
presigned URLs (which expire ~15 minutes after the API call).

Always uses binning=16 (the binning=1 files are too large).
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Any

import boto3

from dasch_sky_mosaic.utils import _decompress_fz, _unpacked_path

logging.basicConfig(level=logging.INFO, format="%(message)s")
LOG = logging.getLogger(__name__)

S3_BUCKET = "dasch-prod-user"
S3_REGION = "us-east-1"

_CREDS_CSV = Path(__file__).parents[2] / "credentials" / "danisimm_plates_accessKeys.csv"

_s3_client: Any = None


def _get_s3_client() -> Any:
    global _s3_client
    if _s3_client is None:
        with open(_CREDS_CSV, newline="", encoding="utf-8-sig") as f:
            row = next(csv.DictReader(f))
        _s3_client = boto3.client(
            "s3",
            aws_access_key_id=row["Access key ID"],
            aws_secret_access_key=row["Secret access key"],
            region_name=S3_REGION,
        )
    return _s3_client


def _s3_download(key: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    response = _get_s3_client().get_object(Bucket=S3_BUCKET, Key=key)
    with dest.open("wb") as f:
        for chunk in response["Body"].iter_chunks(chunk_size=65536):
            f.write(chunk)


def download_photo_from_s3(
    plate_id: str,
    dest_dir: Path,
    overwrite: bool = False,
) -> Path:
    """Download a plate JPG directly from S3."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{plate_id}_pphoto_all.jpg"

    if dest.exists() and not overwrite:
        LOG.info("Reusing cached photo for %s", plate_id)
        return dest

    LOG.info("Downloading photo for %s from S3", plate_id)
    _s3_download(f"plates/{plate_id}/{plate_id}_pphoto_all.jpg", dest)
    return dest


def download_fits_from_s3(
    plate_id: str,
    dest_dir: Path,
    overwrite: bool = False,
) -> Path:
    """Download a binning=16 FITS mosaic directly from S3 and decompress to .fits.

    Lists objects under the plate prefix and matches with a regex to find the
    correct key regardless of rotation/flag suffixes (e.g. _16r270ww, _16ww).
    """
    dest_dir.mkdir(parents=True, exist_ok=True)

    prefix = f"plates/{plate_id}/{plate_id}_mosaic_"
    resp = _get_s3_client().list_objects_v2(Bucket=S3_BUCKET, Prefix=prefix)
    keys = [obj["Key"] for obj in resp.get("Contents", [])]
    matching = [k for k in keys if k.endswith("ww.fit.fz")]

    if not matching:
        raise FileNotFoundError(f"no binning=16 FITS mosaic found for {plate_id} in S3")

    key = matching[0]
    filename = key.split("/")[-1]
    dest_fz = dest_dir / filename
    dest_fits = _unpacked_path(dest_fz)

    if dest_fits.exists() and not overwrite:
        LOG.info("Reusing cached FITS for %s: %s", plate_id, dest_fits.name)
        return dest_fits

    LOG.info("Downloading FITS for %s from S3: %s", plate_id, filename)
    _s3_download(key, dest_fz)
    return _decompress_fz(dest_fz)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Download DASCH plate files from S3")
    parser.add_argument("plate_id", help="Plate ID (e.g. ab12345)")
    parser.add_argument("download", choices=["photo", "fits", "both"], help="What to download")
    parser.add_argument("--dest", default="output", help="Destination directory (default: output)")
    args = parser.parse_args()

    dest_dir = Path(args.dest)
    if args.download in ("photo", "both"):
        path = download_photo_from_s3(args.plate_id, dest_dir)
        print(f"Photo saved to: {path}")
    if args.download in ("fits", "both"):
        path = download_fits_from_s3(args.plate_id, dest_dir)
        print(f"FITS saved to: {path}")
