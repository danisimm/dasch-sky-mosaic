"""Discover and download plate photos (JPGs) for a sky region.

Reuses the same plate discovery logic as the mosaic workflow, then
downloads the full-plate JPG scan for each discovered plate.
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dasch_sky_mosaic.call_sg import (
    _should_log_progress,
    download_photo_via_sg,
)
from dasch_sky_mosaic.mosaic.discover import (
    BuildConfig,
    CandidatePlate,
    Region,
    discover_candidate_plates,
)

LOG = logging.getLogger(__name__)


@dataclass(frozen=True)
class PlatePhotoConfig:
    region: Region
    as_of_jd: float | None
    earliest_jd: float | None
    query_step_deg: float
    api_base: str
    api_key: str | None
    allow_multi_solution_plates: bool
    max_plates: int | None
    output_dir: Path
    manifest_json: Path
    overwrite: bool


def _build_discovery_config(config: PlatePhotoConfig) -> BuildConfig:
    """Wrap PlatePhotoConfig in a BuildConfig for reuse of discover_candidate_plates."""
    placeholder = config.output_dir / "_unused_placeholder.fits"
    return BuildConfig(
        region=config.region,
        as_of_jd=config.as_of_jd,
        earliest_jd=config.earliest_jd,
        session_root=config.output_dir.parent,
        output_fits=placeholder,
        epoch_fits=None,
        manifest_json=config.manifest_json,
        pixel_scale_arcsec=None,
        projection="TAN",
        binning=16,
        query_step_deg=config.query_step_deg,
        api_base=config.api_base,
        api_key=config.api_key,
        subtract_background=False,
        allow_multi_solution_plates=config.allow_multi_solution_plates,
        delete_base_mosaics=False,
        overwrite=True,
        max_plates=config.max_plates,
        from_manifest=None,
    )


def discover_and_download_plate_photos(config: PlatePhotoConfig) -> dict[str, Any]:
    """Discover plates covering the configured region and download their JPGs.

    Returns a JSON-serialisable manifest dict with download status per plate.
    """
    discovery_cfg = _build_discovery_config(config)
    candidates: list[CandidatePlate] = discover_candidate_plates(discovery_cfg)
    if not candidates:
        raise RuntimeError("no candidate plates found for the requested sky region/date constraints")

    config.output_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []

    for idx, plate in enumerate(candidates, start=1):
        if _should_log_progress(idx, len(candidates)):
            LOG.info(
                "Photo download: %d/%d (%.0f%%)",
                idx, len(candidates), 100.0 * idx / max(1, len(candidates)),
            )
        try:
            # Check existence before download to correctly report cached vs new.
            candidate_path = config.output_dir / f"{plate.plate_id}_pphoto_all.jpg"
            already_cached = candidate_path.exists() and not config.overwrite
            local_path = download_photo_via_sg(
                plate_id=plate.plate_id,
                dest_dir=config.output_dir,
                api_base=config.api_base,
                api_key=config.api_key,
                overwrite=config.overwrite,
            )
            status = "cached" if already_cached else "downloaded"
            records.append({
                "plate_id": plate.plate_id,
                "expdate_jd": plate.exp_date_jd,
                "status": status,
                "local_photo_path": str(local_path),
                "bytes": local_path.stat().st_size,
            })
        except Exception as exc:
            LOG.warning("Photo download failed for %s: %s", plate.plate_id, exc)
            records.append({
                "plate_id": plate.plate_id,
                "expdate_jd": plate.exp_date_jd,
                "status": "download_failed",
                "error": str(exc),
            })

    manifest: dict[str, Any] = {
        "workflow": "plate-photo-download",
        "region": {
            "ra_deg": config.region.ra_deg,
            "dec_deg": config.region.dec_deg,
            "width_deg": config.region.width_deg,
            "height_deg": config.region.height_deg,
        },
        "as_of_jd": config.as_of_jd,
        "earliest_jd": config.earliest_jd,
        "query_step_deg": config.query_step_deg,
        "n_candidates": len(candidates),
        "n_downloaded": sum(1 for r in records if r.get("status") in {"downloaded", "cached"}),
        "photos": records,
    }

    config.manifest_json.parent.mkdir(parents=True, exist_ok=True)
    if config.manifest_json.exists() and not config.overwrite:
        raise FileExistsError(f"refusing to overwrite existing file: {config.manifest_json}")
    config.manifest_json.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
