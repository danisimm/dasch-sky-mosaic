from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import requests

from dasch_sky_mosaic.fetch import (
    BuildConfig,
    CandidatePlate,
    Region,
    StarglassClient,
    discover_candidate_plates,
    jd_to_iso,
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
    """Build a fetch.BuildConfig for plate discovery only.

    We reuse the existing discovery machinery and provide inert placeholders for
    FITS-related outputs that are irrelevant for this workflow.
    """
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


def _choose_photo_entry(
    plate_payload: dict[str, Any],
) -> dict[str, Any] | None:
    # Restrict to full plate images only; no jacket images and no thumbnails.
    candidates = [
        c
        for c in (plate_payload.get("plate_images", []) or [])
        if not bool(c.get("thumbnail", False)) and str(c.get("portion", "")).lower() == "all"
    ]

    if not candidates:
        return None
    return max(candidates, key=lambda entry: int(entry.get("thumbnail_ratio") or 0))


def _filename_from_url(url: str, fallback_stem: str) -> str:
    path = urlparse(url).path
    leaf = Path(path).name
    if not leaf:
        return f"{fallback_stem}.jpg"
    return leaf


def _download_file(url: str, dest: Path, overwrite: bool) -> int:
    if dest.exists() and not overwrite:
        return dest.stat().st_size

    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=180.0) as response:
        response.raise_for_status()
        with dest.open("wb") as f_out:
            for chunk in response.iter_content(chunk_size=65536):
                if chunk:
                    f_out.write(chunk)
    return dest.stat().st_size


def discover_and_download_plate_photos(config: PlatePhotoConfig) -> dict[str, Any]:
    discovery_cfg = _build_discovery_config(config)
    candidates: list[CandidatePlate] = discover_candidate_plates(discovery_cfg)

    if not candidates:
        raise RuntimeError("no candidate plates found for requested sky region/date constraints")

    client = StarglassClient(api_base=config.api_base, api_key=config.api_key)
    config.output_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, Any]] = []

    for idx, plate in enumerate(candidates, start=1):
        LOG.info(
            "Photo progress: %d/%d plate metadata fetches (%.0f%%)",
            idx,
            len(candidates),
            (100.0 * idx) / max(1, len(candidates)),
        )
        payload = client.get_plate(plate.plate_id)
        chosen = _choose_photo_entry(payload)

        if chosen is None or not chosen.get("url"):
            records.append(
                {
                    "plate_id": plate.plate_id,
                    "obs_date": jd_to_iso(plate.obs_date_jd),
                    "obs_date_jd": plate.obs_date_jd,
                    "status": "no_full_plate_photo",
                }
            )
            continue

        photo_url = str(chosen["url"])
        filename = _filename_from_url(photo_url, fallback_stem=f"{plate.plate_id}_plate_all")
        local_path = config.output_dir / filename
        existed_before = local_path.exists()
        try:
            n_bytes = _download_file(photo_url, local_path, overwrite=config.overwrite)
            status = "cached" if existed_before and not config.overwrite else "downloaded"
        except Exception as exc:  # pragma: no cover
            records.append(
                {
                    "plate_id": plate.plate_id,
                    "obs_date": jd_to_iso(plate.obs_date_jd),
                    "obs_date_jd": plate.obs_date_jd,
                    "status": "download_failed",
                    "error": str(exc),
                }
            )
            continue

        records.append(
            {
                "plate_id": plate.plate_id,
                "obs_date": jd_to_iso(plate.obs_date_jd),
                "obs_date_jd": plate.obs_date_jd,
                "status": status,
                "portion": chosen.get("portion"),
                "image_type": chosen.get("image_type"),
                "thumbnail": bool(chosen.get("thumbnail", False)),
                "local_photo_path": str(local_path),
                "bytes": n_bytes,
            }
        )

    manifest: dict[str, Any] = {
        "workflow": "plate-photo-download",
        "region": {
            "ra_deg": config.region.ra_deg,
            "dec_deg": config.region.dec_deg,
            "width_deg": config.region.width_deg,
            "height_deg": config.region.height_deg,
        },
        "as_of_date": jd_to_iso(config.as_of_jd) if config.as_of_jd else None,
        "earliest_date": jd_to_iso(config.earliest_jd) if config.earliest_jd else None,
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
