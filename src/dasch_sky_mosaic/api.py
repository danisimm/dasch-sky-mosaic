from __future__ import annotations

import csv
import io
import warnings
from dataclasses import dataclass

import requests
from astropy.time import Time

try:
    from erfa import ErfaWarning
except ImportError:  # pragma: no cover
    from astropy.utils.exceptions import ErfaWarning

PUBLIC_API_BASE = "https://api.starglass.cfa.harvard.edu/public"
FULL_API_BASE = "https://api.starglass.cfa.harvard.edu/full"

# The queryexps API can expose observation time as either obsDate/obs_date
# The public API returns observation time as either obsDate/obs_date (JD float) or
# expdate/expDate (ISO datetime string). Both case variants are included defensively
# since the exact naming has varied across API versions.
_OBS_DATE_KEYS = ("obsDate", "obs_date", "expdate", "expDate")


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
            return Time(v).jd # type: ignore
    except Exception:
        return None


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


@dataclass(frozen=True)
class ExposureHit:
    plate_id: str
    obs_date_jd: float | None  # Julian Date of the exposure midpoint; None if unknown
    solnum: int | None         # WCS solution number; None means no scan-based WCS
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
                wcssource=record.get("wcssource", "").strip().lower(),
                raw=record,
            ))
        return hits