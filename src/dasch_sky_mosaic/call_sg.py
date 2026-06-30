"""Starglass API client and download helpers.

All communication with the Starglass/DASCH HTTP API lives here.
Both plate JPGs and FITS mosaics can be fetched via this module.
"""
from __future__ import annotations

import csv
import io
import logging
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast
from urllib.parse import quote, unquote, urlsplit

import requests
from astropy.time import Time

from dasch_sky_mosaic.utils import _decompress_fz, _download_stream, _unpacked_path

try:
    from erfa import ErfaWarning
except ImportError:
    from astropy.utils.exceptions import ErfaWarning

logging.basicConfig(level=logging.INFO, format="%(message)s")
LOG = logging.getLogger(__name__)

PUBLIC_API_BASE = "https://api.starglass.cfa.harvard.edu/public"
FULL_API_BASE = "https://api.starglass.cfa.harvard.edu/full"

_CREDS_FILE = Path(__file__).parents[2] / "credentials" / "starglass_api_key.txt"


def _load_api_key() -> str | None:
    try:
        return _CREDS_FILE.read_text(encoding="utf-8").strip() or None
    except FileNotFoundError:
        return None
MAX_QUERY_DEC_DEG = 85.0

# ---------------------------------------------------------------------------
# Shared progress helper
# ---------------------------------------------------------------------------

def _should_log_progress(index: int, total: int) -> bool:
    """Log each step for small totals; ~10 checkpoints for larger runs."""
    if total <= 20:
        return True
    stride = max(1, total // 10)
    return index == 1 or index == total or (index % stride == 0)

# ---------------------------------------------------------------------------
# ExposureHit (raw queryexps result row)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ExposureHit:
    plate_id: str
    exp_date_jd: float | None
    solnum: int | None
    exposure_num: int | None
    wcssource: str

    @property
    def has_imaging(self) -> bool:
        return self.solnum is not None and self.wcssource == "imwcs"

def _parse_int_field(row: dict[str, str], key: str) -> int | None:
    raw = row.get(key, "").strip()
    if not raw or raw.lower() in ("nan", "none", "null", "--"):
        return None
    try:
        return int(float(raw))
    except (ValueError, TypeError):
        return None

# ---------------------------------------------------------------------------
# Starglass HTTP client
# ---------------------------------------------------------------------------

class StarglassClient:
    def __init__(
        self,
        api_base: str = "auto",
        api_key: str | None = None,
        timeout_seconds: float = 120.0,
    ) -> None:
        self._timeout = timeout_seconds
        self._session = requests.Session()
        if api_base == "auto":
            self._base_url = FULL_API_BASE if api_key else PUBLIC_API_BASE
        elif api_base == "public":
            self._base_url = PUBLIC_API_BASE
        elif api_base == "full":
            self._base_url = FULL_API_BASE
        else:
            self._base_url = api_base.rstrip("/")
        self._session.headers.update({"accept": "application/json"})
        if api_key:
            self._session.headers.update({"x-api-key": api_key})

    @property
    def base_url(self) -> str:
        return self._base_url

    def query_exposures(self, ra_deg: float, dec_deg: float) -> list[ExposureHit]:
        response = self._session.post(
            f"{self._base_url}/dasch/dr7/queryexps",
            json={"ra_deg": ra_deg, "dec_deg": dec_deg},
            timeout=self._timeout,
        )
        response.raise_for_status()
        hits: list[ExposureHit] = []
        for row in csv.DictReader(io.StringIO("\n".join(response.json()))):
            series = row.get("series", "").strip().lower()
            platenum = row.get("platenum", "").strip()
            if not series or not platenum:
                continue
            try:
                plate_id = f"{series}{int(platenum):05d}"
            except ValueError:
                continue
            raw_date = row.get("expdate", "").strip()
            exp_date_jd: float | None = None
            if raw_date:
                try:
                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore", category=ErfaWarning)
                        exp_date_jd = Time(raw_date).jd  # type: ignore
                except Exception:
                    pass
            hits.append(ExposureHit(
                plate_id=plate_id,
                exp_date_jd=exp_date_jd,
                solnum=_parse_int_field(row, "solnum"),
                exposure_num=_parse_int_field(row, "expnum"),
                wcssource=row.get("wcssource", "").strip().lower(),
            ))
        return hits

    def get_plate(self, plate_id: str) -> dict[str, Any]:
        response = self._session.get(
            f"{self._base_url}/plates/p/{quote(plate_id)}",
            timeout=self._timeout,
        )
        response.raise_for_status()
        return cast(dict[str, Any], response.json())

    def get_mosaic_package(self, plate_id: str, binning: int) -> dict[str, Any]:
        response = self._session.post(
            f"{self._base_url}/dasch/dr7/mosaic_package",
            json={"plate_id": plate_id, "binning": binning},
            timeout=self._timeout,
        )
        response.raise_for_status()
        return cast(dict[str, Any], response.json())

# ---------------------------------------------------------------------------
# FITS download via Starglass mosaic_package
# ---------------------------------------------------------------------------

def download_fits_via_sg(
    plate_id: str,
    dest_dir: Path,
) -> Path:
    """Download a binning=16 FITS mosaic via the Starglass mosaic_package endpoint."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    api_key = _load_api_key()
    client = StarglassClient(api_base="auto", api_key=api_key)
    package = client.get_mosaic_package(plate_id=plate_id, binning=16)
    base_url = str(package.get("baseFitsUrl") or "").strip()
    if not base_url:
        raise RuntimeError(f"no FITS URL in mosaic_package response for {plate_id}")

    remote_name = Path(unquote(urlsplit(base_url).path)).name
    if remote_name:
        local_name = remote_name if remote_name.lower().startswith(plate_id.lower()) else f"{plate_id}_{remote_name}"
    else:
        local_name = f"{plate_id}_bin16.fits"

    local_path = dest_dir / local_name
    if local_path.suffix.lower() == ".fz":
        unpacked = _unpacked_path(local_path)
        if unpacked.exists():
            LOG.info("Reusing cached FITS for %s: %s", plate_id, unpacked.name)
            return unpacked

    if not local_path.exists():
        LOG.info("Downloading FITS mosaic for %s via Starglass", plate_id)
        _download_stream(base_url, local_path)
    else:
        LOG.info("Reusing cached .fz for %s: %s", plate_id, local_path.name)

    return _decompress_fz(local_path)

# ---------------------------------------------------------------------------
# Photo (JPG) download via Starglass plates API
# ---------------------------------------------------------------------------

def download_photo_via_sg(
    plate_id: str,
    dest_dir: Path,
) -> Path:
    """Download a plate JPG via the Starglass plates API."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{plate_id}_pphoto_all.jpg"

    if dest.exists():
        LOG.info("Reusing cached photo for %s", plate_id)
        return dest

    api_key = _load_api_key()
    client = StarglassClient(api_base="auto", api_key=api_key)
    payload = client.get_plate(plate_id)
    images = payload.get("plate_images", [])
    if not images or not images[0].get("url"):
        raise RuntimeError(f"no full-plate photo available for {plate_id}")

    LOG.info("Downloading photo for %s via Starglass", plate_id)
    _download_stream(str(images[0]["url"]), dest)
    return dest


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Download DASCH plate files via Starglass API")
    parser.add_argument("plate_id", help="Plate ID (e.g. ab12345)")
    parser.add_argument("download", choices=["photo", "fits", "both"], help="What to download")
    parser.add_argument("--dest", default="data/cache", help="Destination directory (default: data/cache)")
    args = parser.parse_args()

    dest_dir = Path(args.dest)
    if args.download in ("photo", "both"):
        path = download_photo_via_sg(args.plate_id, dest_dir)
        print(f"Photo saved to: {path}")
    if args.download in ("fits", "both"):
        path = download_fits_via_sg(args.plate_id, dest_dir)
        print(f"FITS saved to: {path}")
