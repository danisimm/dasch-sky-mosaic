# This file was generated with the assistance of GitHub Copilot (AI).
from __future__ import annotations

import csv
import io
import logging
import math
import os
import warnings
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast
from urllib.parse import quote

import numpy as np
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from daschlab import open_session

try:
    from erfa import ErfaWarning
except ImportError:  # pragma: no cover
    from astropy.utils.exceptions import ErfaWarning

LOG = logging.getLogger(__name__)

PUBLIC_API_BASE = "https://api.starglass.cfa.harvard.edu/public"
FULL_API_BASE = "https://api.starglass.cfa.harvard.edu/full"
MAX_QUERY_DEC_DEG = 85.0

# The public API returns observation time as either obsDate/obs_date (JD float) or
# expdate/expDate (ISO datetime string). Both case variants are included defensively
# since the exact naming has varied across API versions.
_OBS_DATE_KEYS = ("obsDate", "obs_date", "expdate", "expDate")


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Region:
    ra_deg: float
    dec_deg: float
    width_deg: float
    height_deg: float


@dataclass(frozen=True)
class CandidatePlate:
    plate_id: str
    obs_date_jd: float
    n_wcs_solutions: int
    n_exposures: int
    preferred_solution_num: int | None
    selected_at_points: int


@dataclass(frozen=True)
class BuildConfig:
    region: Region
    as_of_jd: float | None
    earliest_jd: float | None
    session_root: Path
    output_fits: Path
    epoch_fits: Path | None
    manifest_json: Path
    pixel_scale_arcsec: float | None
    projection: str
    binning: int
    query_step_deg: float
    api_base: str
    api_key: str | None
    subtract_background: bool
    allow_multi_solution_plates: bool
    delete_base_mosaics: bool
    overwrite: bool
    max_plates: int | None
    from_manifest: Path | None = None


# ---------------------------------------------------------------------------
# Date utilities
# ---------------------------------------------------------------------------

def jd_to_iso(jd: float) -> str:
    return Time(jd, format="jd", scale="utc").isot  # type: ignore


def parse_obs_date_jd(raw_value: str) -> float | None:
    """Parse an obs_date value from the queryexps CSV response into a Julian Date.

    The value may be an ISO-format date string (Astropy serialisation) or already
    a Julian Date float.  Returns None if the value is blank or cannot be parsed.
    """
    v = raw_value.strip() if raw_value else ""
    if not v or v.lower() in ("nan", "none", "null", "--"):
        return None
    # Try plain float first – real JDs for DASCH plates (~1880-1990) are above 2e6
    try:
        jd = float(v)
        if jd > 1_000_000.0:
            return jd
    except ValueError:
        pass
    # Fall back to Astropy string parsing
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ErfaWarning, message=r".*dubious year.*")
            warnings.filterwarnings("ignore", category=ErfaWarning, message=r".*time is after end of day.*")
            warnings.filterwarnings("ignore", category=ErfaWarning, message=r".*both of next two.*")
            return Time(v).jd  # type: ignore
    except Exception:
        return None


def parse_cli_date_jd(value: str | None) -> float | None:
    """Parse a user-supplied date string into a Julian Date float.

    Accepts ISO format (e.g. '1950-12-31') or a raw Julian Date number.
    Raises ValueError for unrecognised input.
    """
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
# Starglass HTTP API client
# ---------------------------------------------------------------------------

def _normalize_api_base(api_base: str, api_key: str | None) -> str:
    if api_base == "auto":
        return FULL_API_BASE if api_key else PUBLIC_API_BASE
    if api_base == "public":
        return PUBLIC_API_BASE
    if api_base == "full":
        return FULL_API_BASE
    return api_base.rstrip("/")


def _parse_csv_records(lines: list[str]) -> list[dict[str, str]]:
    if not lines:
        return []
    buffer = io.StringIO("\n".join(lines))
    return list(csv.DictReader(buffer))


def _plate_id_from_row(row: dict[str, str]) -> str:
    if row.get("plateId"):
        return row["plateId"].strip().lower()
    series = row.get("series", "").strip().lower()
    platenum = row.get("platenum", "").strip()
    if not series or not platenum:
        raise KeyError("queryexps response did not include enough information to construct a plate_id")
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


@dataclass(frozen=True)
class ExposureHit:
    plate_id: str
    obs_date_jd: float | None  # Julian Date of the exposure midpoint; None if unknown
    solnum: int | None         # WCS solution number; None means no scan-based WCS
    exposure_num: int | None   # Exposure number; None when unavailable
    wcssource: str             # "imwcs" | "logbook" | "none"
    raw: dict[str, str]

    @property
    def has_imaging(self) -> bool:
        """True when this exposure has a WCS solution derived from the plate image."""
        return self.solnum is not None and self.wcssource == "imwcs"


class StarglassClient:
    def __init__(
        self,
        api_base: str = "auto",
        api_key: str | None = None,
        timeout_seconds: float = 120.0,
    ) -> None:
        self._timeout_seconds = timeout_seconds
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
            timeout=self._timeout_seconds,
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
                raw=record,
            ))
        return hits

    def get_plate(self, plate_id: str) -> dict[str, Any]:
        response = self._session.get(
            f"{self._base_url}/plates/p/{quote(plate_id)}",
            timeout=self._timeout_seconds,
        )
        response.raise_for_status()
        return cast(dict[str, Any], response.json())


# ---------------------------------------------------------------------------
# Plate discovery
# ---------------------------------------------------------------------------

def _normalize_ra_deg(value: float) -> float:
    return value % 360.0


def _clip_dec_deg(value: float) -> float:
    return max(-89.999999, min(89.999999, value))


def _grid_axis(half_width_deg: float, spacing_deg: float) -> list[float]:
    if half_width_deg == 0:
        return [0.0]
    n_steps = max(1, math.ceil((2.0 * half_width_deg) / spacing_deg))
    axis = np.linspace(-half_width_deg, half_width_deg, num=n_steps + 1)
    return [float(v) for v in axis]


def iter_query_points(region: Region, query_step_deg: float) -> list[tuple[float, float]]:
    center = SkyCoord(ra=region.ra_deg * u.deg, dec=region.dec_deg * u.deg, frame="icrs")
    offset_frame = center.skyoffset_frame()
    lon_offsets = _grid_axis(region.width_deg / 2.0, query_step_deg)
    lat_offsets = _grid_axis(region.height_deg / 2.0, query_step_deg)
    points: list[tuple[float, float]] = []

    for lat_offset in lat_offsets:
        for lon_offset in lon_offsets:
            point = SkyCoord(lon=lon_offset * u.deg, lat=lat_offset * u.deg, frame=offset_frame).transform_to("icrs")
            point_ra = cast(Any, point.ra)
            point_dec = cast(Any, point.dec)
            points.append((_normalize_ra_deg(float(point_ra.deg)), _clip_dec_deg(float(point_dec.deg))))

    unique: dict[tuple[int, int], tuple[float, float]] = {}
    for ra_deg, dec_deg in points:
        if dec_deg > MAX_QUERY_DEC_DEG:
            continue
        unique[(round(ra_deg * 1000), round(dec_deg * 1000))] = (ra_deg, dec_deg)
    return list(unique.values())


def discover_candidate_plates(config: BuildConfig) -> list[CandidatePlate]:
    """Select plates using only queryexps data — no per-plate get_plate() calls.

    For each grid point the single most recently observed plate with
    wcssource=='imwcs' that satisfies the date bounds is chosen as the
    winner for that point.  Only the union of those per-point winners is
    ever downloaded.
    """
    client = StarglassClient(api_base=config.api_base, api_key=config.api_key)
    query_points = iter_query_points(config.region, config.query_step_deg)
    LOG.info("Querying %d sky positions via %s", len(query_points), client.base_url)

    plate_exposures: dict[str, list[tuple[float | None, int | None, int | None]]] = {}
    point_winner: dict[tuple[int, int], str] = {}
    point_winner_solnum: dict[tuple[int, int], int | None] = {}

    for idx, (ra_deg, dec_deg) in enumerate(query_points, start=1):
        if _should_log_progress(idx, len(query_points)):
            LOG.info(
                "Discovery progress: %d/%d query points (%.0f%%)",
                idx,
                len(query_points),
                (100.0 * idx) / max(1, len(query_points)),
            )
        hits = client.query_exposures(ra_deg=ra_deg, dec_deg=dec_deg)
        eligible = [h for h in hits if h.has_imaging]
        if config.as_of_jd is not None:
            eligible = [h for h in eligible if h.obs_date_jd is not None and h.obs_date_jd <= config.as_of_jd]
        if config.earliest_jd is not None:
            eligible = [h for h in eligible if h.obs_date_jd is not None and h.obs_date_jd >= config.earliest_jd]
        for h in eligible:
            plate_exposures.setdefault(h.plate_id, []).append((h.obs_date_jd, h.solnum, h.exposure_num))
        best = max(
            eligible,
            key=lambda h: h.obs_date_jd if h.obs_date_jd is not None else -1.0,
            default=None,
        )
        if best is not None:
            key = (round(ra_deg * 1000), round(dec_deg * 1000))
            point_winner[key] = best.plate_id
            point_winner_solnum[key] = best.solnum

    selection_counts = Counter(point_winner.values())
    LOG.info(
        "Discovered %d candidate plates; %d selected as most-recent at >= 1 query point",
        len(plate_exposures),
        len(selection_counts),
    )

    candidates: list[CandidatePlate] = []
    for plate_id, point_count in selection_counts.items():
        exposures = plate_exposures.get(plate_id, [])
        solnums = {s for _, s, _ in exposures if s is not None}
        n_wcs = len(solnums)
        if not config.allow_multi_solution_plates and n_wcs > 1:
            continue

        explicit_exposure_nums = {exp for _, _, exp in exposures if exp is not None}
        if explicit_exposure_nums:
            n_exposures = len(explicit_exposure_nums)
        else:
            # queryexps exposure IDs are not always present; approximate via
            # distinct (obs_date_jd, solnum) combinations.
            proxy_ids = {
                (round(jd, 6), sol)
                for jd, sol, _ in exposures
                if jd is not None and sol is not None
            }
            n_exposures = max(1, len(proxy_ids))

        if n_exposures > 1:
            LOG.info(
                "Skipping %s: plate has %d exposures",
                plate_id,
                n_exposures,
            )
            continue

        valid_jds = [jd for jd, _, _ in exposures if jd is not None]
        if not valid_jds:
            continue

        preferred_counter = Counter(
            point_winner_solnum[key]
            for key, winner_plate in point_winner.items()
            if winner_plate == plate_id and point_winner_solnum.get(key) is not None
        )
        preferred_solution_num = preferred_counter.most_common(1)[0][0] if preferred_counter else None

        candidates.append(CandidatePlate(
            plate_id=plate_id,
            obs_date_jd=max(valid_jds),
            n_wcs_solutions=n_wcs,
            n_exposures=n_exposures,
            preferred_solution_num=preferred_solution_num,
            selected_at_points=point_count,
        ))

    candidates.sort(key=lambda c: c.obs_date_jd)

    if config.max_plates is not None and len(candidates) > config.max_plates:
        candidates = sorted(candidates, key=lambda c: -c.obs_date_jd)[: config.max_plates]
        candidates.sort(key=lambda c: c.obs_date_jd)

    LOG.info(
        "Final plate count: %d (allow_multi_solution_plates=%s)",
        len(candidates),
        config.allow_multi_solution_plates,
    )
    return candidates


# ---------------------------------------------------------------------------
# Mosaic downloads via daschlab
# ---------------------------------------------------------------------------

def download_mosaic_paths(
    session_root: Path,
    plate_ids: list[str],
    binning: int,
    api_base: str,
    api_key: str | None,
) -> dict[str, Path]:
    """Download DASCH FITS mosaics and return a plate_id -> local Path map."""
    session_root.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    session = open_session(str(session_root))
    for idx, plate_id in enumerate(plate_ids, start=1):
        if _should_log_progress(idx, len(plate_ids)):
            LOG.info(
                "Download progress: %d/%d plate mosaics (%.0f%%)",
                idx,
                len(plate_ids),
                (100.0 * idx) / max(1, len(plate_ids)),
            )
        LOG.info("Fetching DASCH mosaic for %s", plate_id)
        relpath = session.mosaic(plate_id, binning)
        if relpath is None:
            LOG.warning("No mosaic available for %s; skipping", plate_id)
            continue
        paths[plate_id] = Path(session.path(relpath))

    return paths
