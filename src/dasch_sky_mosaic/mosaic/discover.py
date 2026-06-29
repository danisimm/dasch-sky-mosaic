"""Plate discovery via the DASCH queryexps API.

Given a sky region and date bounds, queries a grid of sky positions and
selects the best set of plates to cover the region.
"""
from __future__ import annotations

import logging
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from dasch_sky_mosaic.call_sg import (
    MAX_QUERY_DEC_DEG,
    ExposureHit,
    StarglassClient,
    _should_log_progress,
)

LOG = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
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
# Grid helpers
# ---------------------------------------------------------------------------

def _normalize_ra_deg(value: float) -> float:
    return value % 360.0


def _clip_dec_deg(value: float) -> float:
    return max(-89.999999, min(89.999999, value))


def _query_point_key(ra_deg: float, dec_deg: float) -> tuple[int, int]:
    return (round(ra_deg * 1000), round(dec_deg * 1000))


def _grid_axis(half_width_deg: float, spacing_deg: float) -> list[float]:
    if half_width_deg == 0 or spacing_deg >= 2.0 * half_width_deg:
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

    for lat in lat_offsets:
        for lon in lon_offsets:
            pt = SkyCoord(lon=lon * u.deg, lat=lat * u.deg, frame=offset_frame).transform_to("icrs")
            points.append((
                _normalize_ra_deg(float(cast(Any, pt.ra).deg)),
                _clip_dec_deg(float(cast(Any, pt.dec).deg)),
            ))

    unique: dict[tuple[int, int], tuple[float, float]] = {}
    for ra_deg, dec_deg in points:
        if dec_deg > MAX_QUERY_DEC_DEG:
            continue
        unique[_query_point_key(ra_deg, dec_deg)] = (ra_deg, dec_deg)
    return list(unique.values())


# ---------------------------------------------------------------------------
# Plate discovery
# ---------------------------------------------------------------------------

def discover_candidate_plates(config: BuildConfig) -> list[CandidatePlate]:
    """Select plates using only queryexps data — no per-plate API calls.

    For each grid point the single most recently observed plate with
    wcssource=='imwcs' that satisfies the date bounds is chosen.
    The union of those per-point winners becomes the candidate set.
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
                "Discovery: %d/%d points (%.0f%%)",
                idx, len(query_points), 100.0 * idx / max(1, len(query_points)),
            )
        hits = client.query_exposures(ra_deg=ra_deg, dec_deg=dec_deg)
        eligible: list[ExposureHit] = []
        for hit in hits:
            if not hit.has_imaging:
                continue
            obs_jd = hit.obs_date_jd
            if obs_jd is None:
                continue
            if config.as_of_jd is not None and obs_jd > config.as_of_jd:
                continue
            if config.earliest_jd is not None and obs_jd < config.earliest_jd:
                continue
            eligible.append(hit)

        for h in eligible:
            plate_exposures.setdefault(h.plate_id, []).append((h.obs_date_jd, h.solnum, h.exposure_num))

        best = max(
            eligible,
            key=lambda h: h.obs_date_jd if h.obs_date_jd is not None else -1.0,
            default=None,
        )
        if best is not None:
            key = _query_point_key(ra_deg, dec_deg)
            point_winner[key] = best.plate_id
            point_winner_solnum[key] = best.solnum

    selection_counts = Counter(point_winner.values())
    LOG.info(
        "Discovered %d candidate plates; %d selected as most-recent at >= 1 point",
        len(plate_exposures), len(selection_counts),
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
            proxy_ids = {
                (round(jd, 6), sol)
                for jd, sol, _ in exposures
                if jd is not None and sol is not None
            }
            n_exposures = max(1, len(proxy_ids))

        if n_exposures > 1:
            LOG.info("Skipping %s: plate has %d exposures", plate_id, n_exposures)
            continue

        valid_jds = [jd for jd, _, _ in exposures if jd is not None]
        if not valid_jds:
            continue

        preferred_counter = Counter(
            point_winner_solnum[key]
            for key, winner in point_winner.items()
            if winner == plate_id and point_winner_solnum.get(key) is not None
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

    LOG.info("Final plate count: %d", len(candidates))
    return candidates
