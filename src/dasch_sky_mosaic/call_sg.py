"""Starglass API client and download helpers.

All communication with the Starglass/DASCH HTTP API lives here.
Both plate JPGs and FITS mosaics can be fetched via this module.
"""
from __future__ import annotations

import csv
import io
import logging
import os
import shutil
import subprocess
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast
from urllib.parse import quote, unquote, urlparse, urlsplit

import numpy as np
import requests
from astropy.io import fits
from astropy.time import Time

try:
    from erfa import ErfaWarning
except ImportError:
    from astropy.utils.exceptions import ErfaWarning

LOG = logging.getLogger(__name__)

PUBLIC_API_BASE = "https://api.starglass.cfa.harvard.edu/public"
FULL_API_BASE = "https://api.starglass.cfa.harvard.edu/full"
MAX_QUERY_DEC_DEG = 85.0

# The queryexps API has returned obs date under several key names across versions.
_OBS_DATE_KEYS = ("obsDate", "obs_date", "expdate", "expDate")


# ---------------------------------------------------------------------------
# Date utilities
# ---------------------------------------------------------------------------

def jd_to_iso(jd: float) -> str:
    return Time(jd, format="jd", scale="utc").isot  # type: ignore


def parse_obs_date_jd(raw_value: str) -> float | None:
    """Parse an obs_date value from the queryexps API into a Julian Date.

    Accepts a JD float string or an ISO datetime string. Returns None if blank
    or unparseable.
    """
    v = raw_value.strip() if raw_value else ""
    if not v or v.lower() in ("nan", "none", "null", "--"):
        return None
    try:
        jd = float(v)
        if jd > 1_000_000.0:
            return jd
    except ValueError:
        pass
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ErfaWarning)
            return Time(v).jd  # type: ignore
    except Exception:
        return None


def parse_cli_date_jd(value: str | None) -> float | None:
    """Parse a user-supplied date string ('1950-12-31' or JD) into a JD float."""
    if not value:
        return None
    jd = parse_obs_date_jd(value)
    if jd is None:
        raise ValueError(f"Cannot parse date {value!r}. Use ISO format, e.g. 1950-12-31.")
    return jd


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
# Internal CSV / row parsing helpers
# ---------------------------------------------------------------------------

def _parse_csv_records(lines: list[str]) -> list[dict[str, str]]:
    if not lines:
        return []
    return list(csv.DictReader(io.StringIO("\n".join(lines))))


def _plate_id_from_row(row: dict[str, str]) -> str:
    if row.get("plateId"):
        return row["plateId"].strip().lower()
    series = row.get("series", "").strip().lower()
    platenum = row.get("platenum", "").strip()
    if not series or not platenum:
        raise KeyError("queryexps response missing plate identity fields")
    return f"{series}{int(platenum):05d}"


def _extract_obs_date_jd(row: dict[str, str]) -> float | None:
    for key in _OBS_DATE_KEYS:
        raw = row.get(key, "")
        if raw:
            result = parse_obs_date_jd(raw)
            if result is not None:
                return result
    return None


def _extract_solnum(row: dict[str, str]) -> int | None:
    raw = row.get("solnum", "").strip()
    if not raw or raw.lower() in ("nan", "none", "null", "--"):
        return None
    try:
        return int(float(raw))
    except (ValueError, TypeError):
        return None


def _extract_exposure_num(row: dict[str, str]) -> int | None:
    for key in ("exposure_num", "exposureNum", "expnum", "exp_num"):
        raw = row.get(key, "").strip()
        if not raw or raw.lower() in ("nan", "none", "null", "--"):
            continue
        try:
            return int(float(raw))
        except (ValueError, TypeError):
            continue
    return None


def _normalize_api_base(api_base: str, api_key: str | None) -> str:
    if api_base == "auto":
        return FULL_API_BASE if api_key else PUBLIC_API_BASE
    if api_base == "public":
        return PUBLIC_API_BASE
    if api_base == "full":
        return FULL_API_BASE
    return api_base.rstrip("/")


# ---------------------------------------------------------------------------
# ExposureHit (raw queryexps result row)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ExposureHit:
    plate_id: str
    obs_date_jd: float | None
    solnum: int | None
    exposure_num: int | None
    wcssource: str

    @property
    def has_imaging(self) -> bool:
        return self.solnum is not None and self.wcssource == "imwcs"


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
        self._base_url = _normalize_api_base(api_base, api_key)
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
        for record in _parse_csv_records(response.json()):
            try:
                plate_id = _plate_id_from_row(record)
            except (KeyError, ValueError):
                continue
            hits.append(ExposureHit(
                plate_id=plate_id,
                obs_date_jd=_extract_obs_date_jd(record),
                solnum=_extract_solnum(record),
                exposure_num=_extract_exposure_num(record),
                wcssource=record.get("wcssource", "").strip().lower(),
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
# Generic streaming download
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# FITS download via Starglass mosaic_package
# ---------------------------------------------------------------------------

def _funpack_if_needed(path: Path) -> Path:
    """Return an uncompressed FITS path, funpacking .fz files when necessary."""
    if path.suffix.lower() != ".fz":
        return path

    out_path = path.with_suffix("")
    if out_path.suffix.lower() == ".fit":
        out_path = out_path.with_suffix(".fits")
    elif out_path.suffix.lower() != ".fits":
        out_path = out_path.with_name(out_path.name + ".fits")

    if out_path.exists():
        return out_path

    funpack_exe = shutil.which("funpack")
    if funpack_exe:
        LOG.info("Decompressing %s with funpack", path.name)
        result = subprocess.run(
            [funpack_exe, "-O", os.fspath(out_path), os.fspath(path)],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode == 0 and out_path.exists():
            return out_path
        detail = (result.stderr or result.stdout or f"exit {result.returncode}").strip()
        LOG.warning("funpack failed for %s (%s); falling back to astropy", path.name, detail)

    LOG.info("Decompressing %s with astropy", path.name)
    with fits.open(path, memmap=False) as hdul:
        image_hdu = None
        for hdu in hdul:
            raw = getattr(hdu, "data", None)
            if raw is None:
                continue
            arr = np.asarray(raw)
            if arr.ndim < 2:
                continue
            if arr.ndim > 2:
                arr = np.squeeze(arr)
                if arr.ndim != 2:
                    continue
            image_hdu = fits.PrimaryHDU(
                data=np.array(arr, copy=True),
                header=cast(Any, hdu).header.copy(),
            )
            break
        if image_hdu is None:
            raise RuntimeError(f"cannot decompress {path}: no 2D image HDU found")
        fits.HDUList([image_hdu]).writeto(out_path, overwrite=True)
    return out_path


def download_fits_via_sg(
    plate_id: str,
    dest_dir: Path,
    binning: int,
    api_base: str,
    api_key: str | None,
) -> Path:
    """Download a single FITS mosaic via the Starglass mosaic_package endpoint."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    client = StarglassClient(api_base=api_base, api_key=api_key)
    package = client.get_mosaic_package(plate_id=plate_id, binning=binning)
    base_url = str(
        package.get("base_fits_url") or package.get("baseFitsUrl") or ""
    ).strip()
    if not base_url:
        raise RuntimeError(f"no FITS URL in mosaic_package response for {plate_id}")

    remote_name = Path(unquote(urlsplit(base_url).path)).name
    if remote_name:
        local_name = remote_name if remote_name.lower().startswith(plate_id.lower()) else f"{plate_id}_{remote_name}"
    else:
        local_name = f"{plate_id}_bin{binning}.fits"

    local_path = dest_dir / local_name
    if not local_path.exists():
        LOG.info("Downloading FITS mosaic for %s via Starglass", plate_id)
        _download_stream(base_url, local_path)
    else:
        LOG.info("Reusing cached FITS for %s: %s", plate_id, local_path.name)

    return _funpack_if_needed(local_path)


def download_fits_batch_via_sg(
    plate_ids: list[str],
    dest_dir: Path,
    binning: int,
    api_base: str,
    api_key: str | None,
) -> dict[str, Path]:
    """Download FITS mosaics for a list of plate IDs. Returns {plate_id: local_path}."""
    paths: dict[str, Path] = {}
    for idx, plate_id in enumerate(plate_ids, start=1):
        if _should_log_progress(idx, len(plate_ids)):
            LOG.info(
                "FITS download progress: %d/%d (%.0f%%)",
                idx, len(plate_ids), 100.0 * idx / max(1, len(plate_ids)),
            )
        try:
            paths[plate_id] = download_fits_via_sg(plate_id, dest_dir, binning, api_base, api_key)
        except Exception as exc:
            LOG.warning("FITS download failed for %s: %s", plate_id, exc)
    return paths


# ---------------------------------------------------------------------------
# Photo (JPG) download via Starglass plates API
# ---------------------------------------------------------------------------

def _choose_photo_entry(plate_payload: dict[str, Any]) -> dict[str, Any] | None:
    """Select the best full-plate non-thumbnail image from plate metadata."""
    candidates = [
        c
        for c in (plate_payload.get("plate_images", []) or [])
        if not bool(c.get("thumbnail", False)) and str(c.get("portion", "")).lower() == "all"
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda e: int(e.get("thumbnail_ratio") or 0))


def download_photo_via_sg(
    plate_id: str,
    dest_dir: Path,
    api_base: str,
    api_key: str | None,
    overwrite: bool = False,
) -> Path:
    """Download a plate JPG via the Starglass plates API."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    client = StarglassClient(api_base=api_base, api_key=api_key)
    payload = client.get_plate(plate_id)
    chosen = _choose_photo_entry(payload)
    if chosen is None or not chosen.get("url"):
        raise RuntimeError(f"no full-plate photo available for {plate_id}")

    photo_url = str(chosen["url"])
    leaf = Path(urlparse(photo_url).path).name
    filename = leaf if leaf else f"{plate_id}_pphoto_all.jpg"
    dest = dest_dir / filename

    if dest.exists() and not overwrite:
        LOG.info("Reusing cached photo for %s: %s", plate_id, dest.name)
        return dest

    LOG.info("Downloading photo for %s via Starglass", plate_id)
    _download_stream(photo_url, dest, overwrite=overwrite)
    return dest
