from __future__ import annotations

import json
import logging
import math
import os
import shutil
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from reproject import reproject_interp

from dasch_sky_mosaic.background import (
    _estimate_overlap_template_refinement,
    _evaluate_residual_surface,
    _fit_plate_background,
    _fit_overlap_residual_surface,
    _plate_interior_mask,
    _solve_global_bg_offsets,
)
from dasch_sky_mosaic.fetch import (
    BuildConfig,
    CandidatePlate,
    Region,
    _should_log_progress,
    discover_candidate_plates,
    download_mosaic_paths,
    jd_to_iso,
)

LOG = logging.getLogger(__name__)
_RESIDUAL_SURFACE_DEGREE = 1


def _build_output_wcs(region: Region, projection: str, pixel_scale_arcsec: float) -> tuple[WCS, tuple[int, int]]:
    pixel_scale_deg = pixel_scale_arcsec / 3600.0
    naxis1 = max(1, math.ceil(region.width_deg / pixel_scale_deg))
    naxis2 = max(1, math.ceil(region.height_deg / pixel_scale_deg))

    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [(naxis1 + 1) / 2.0, (naxis2 + 1) / 2.0]
    wcs.wcs.crval = [region.ra_deg, region.dec_deg]
    wcs.wcs.cdelt = np.array([-pixel_scale_deg, pixel_scale_deg])
    wcs.wcs.ctype = [f"RA---{projection}", f"DEC--{projection}"]
    return wcs, (naxis2, naxis1)


def _native_pixel_scale_arcsec(mosaic_path: Path) -> float:
    with fits.open(mosaic_path, memmap=True) as hdul:
        primary_hdu = hdul[0]  # type: ignore[index]
        header = primary_hdu.header.copy()  # type: ignore[attr-defined]
        wcs = WCS(header)
    scales_deg = proj_plane_pixel_scales(wcs)
    scale_arcsec = float(np.mean(np.abs(scales_deg)) * 3600.0)
    if not np.isfinite(scale_arcsec) or scale_arcsec <= 0:
        raise RuntimeError(f"cannot determine native pixel scale from {mosaic_path}")
    return scale_arcsec


def _prepare_output_path(path: Path, overwrite: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite existing file: {path}")


def _write_fits(path: Path, data: np.ndarray, header: fits.Header, overwrite: bool) -> None:
    _prepare_output_path(path, overwrite)
    fits.PrimaryHDU(data=data.astype(np.float32), header=header).writeto(path, overwrite=overwrite)


def _candidates_from_manifest(manifest_path: Path) -> tuple[list[CandidatePlate], dict[str, Path]]:
    """Load plate candidates and local paths from a prior build manifest, bypassing discovery and download."""
    import json as _json
    data = _json.loads(manifest_path.read_text(encoding="utf-8"))
    candidates: list[CandidatePlate] = []
    mosaic_paths: dict[str, Path] = {}
    for p in data.get("plates", []):
        plate_id = p["plate_id"]
        local_path = Path(p["local_mosaic_path"])
        if not local_path.exists():
            LOG.warning("Cached mosaic not found, skipping: %s", local_path)
            continue
        candidates.append(CandidatePlate(
            plate_id=plate_id,
            obs_date_jd=p["obs_date_jd"],
            n_wcs_solutions=p.get("n_wcs_solutions", 1),
            selected_at_points=p.get("selected_at_points", 1),
        ))
        mosaic_paths[plate_id] = local_path
    if not candidates:
        raise RuntimeError(f"no valid plate entries found in manifest {manifest_path}")
    LOG.info("Loaded %d plates from manifest %s", len(candidates), manifest_path.name)
    return candidates, mosaic_paths


def build_mosaic(config: BuildConfig) -> dict[str, Any]:
    if config.from_manifest is not None:
        candidates, mosaic_paths = _candidates_from_manifest(config.from_manifest)
    else:
        candidates = discover_candidate_plates(config)
        if not candidates:
            raise RuntimeError("no DASCH plates matched the requested sky region and date constraints")

        mosaic_paths = download_mosaic_paths(
            session_root=config.session_root,
            plate_ids=[c.plate_id for c in candidates],
            binning=config.binning,
        )
        candidates = [c for c in candidates if c.plate_id in mosaic_paths]
        if not candidates:
            raise RuntimeError("no mosaics could be downloaded for the selected plates")

    pixel_scale_arcsec = config.pixel_scale_arcsec
    if pixel_scale_arcsec is None:
        reference = candidates[-1]
        reference_path = mosaic_paths[reference.plate_id]
        pixel_scale_arcsec = _native_pixel_scale_arcsec(reference_path)
        LOG.info(
            "Using native plate scale %.3f arcsec/pixel from %s",
            pixel_scale_arcsec,
            reference_path.name,
        )

    output_wcs, shape_out = _build_output_wcs(
        region=config.region,
        projection=config.projection,
        pixel_scale_arcsec=pixel_scale_arcsec,
    )

    # Pass 1: reproject all plates and collect data for global background matching.
    all_reprojected: list[np.ndarray] = []
    all_good: list[np.ndarray] = []
    all_templates: list[np.ndarray] = []
    all_obs_jds: list[float] = []
    all_names: list[str] = []

    for idx, candidate in enumerate(candidates, start=1):
        if _should_log_progress(idx, len(candidates)):
            LOG.info(
                "Reprojection progress: %d/%d plate mosaics (%.0f%%)",
                idx,
                len(candidates),
                (100.0 * idx) / max(1, len(candidates)),
            )
        mosaic_path = mosaic_paths[candidate.plate_id]
        LOG.info(
            "Reprojecting %s  obs %s  (JD %.1f)",
            mosaic_path.name,
            jd_to_iso(candidate.obs_date_jd),
            candidate.obs_date_jd,
        )
        with fits.open(mosaic_path, memmap=True) as hdul:
            plate_header = hdul[0].header.copy() # type: ignore
            image = np.array(hdul[0].data, dtype=np.float32, copy=True) # type: ignore

        # Build interior support from raw plate values before any intensity transforms.
        interior_native = _plate_interior_mask(image, binning=config.binning).astype(np.float32)

        image[~np.isfinite(image)] = np.nan
        if config.subtract_background:
            bg_model = _fit_plate_background(image, degree=2)
            image -= bg_model
            bg_template_native = bg_model - np.nanmedian(bg_model)
            LOG.info(
                "Background model subtracted for %s (centre value %.1f)",
                mosaic_path.name,
                float(bg_model[bg_model.shape[0] // 2, bg_model.shape[1] // 2]),
            )
        else:
            bg_template_native = np.zeros_like(image, dtype=np.float32)

        reprojected, footprint = reproject_interp(
            (image, plate_header),
            output_projection=output_wcs,
            shape_out=shape_out,
            return_footprint=True,
        )

        interior_proj, _ = reproject_interp(
            (interior_native, plate_header),
            output_projection=output_wcs,
            shape_out=shape_out,
            return_footprint=True,
        )
        template_proj, _ = reproject_interp(
            (bg_template_native, plate_header),
            output_projection=output_wcs,
            shape_out=shape_out,
            return_footprint=True,
        )

        # Reprojected binary masks are interpolated, so require strong-but-not-perfect support.
        good = np.isfinite(reprojected) & (footprint > 0.2) & np.isfinite(interior_proj) & (interior_proj > 0.1)

        all_reprojected.append(reprojected)
        all_good.append(good)
        all_templates.append(template_proj)
        all_obs_jds.append(candidate.obs_date_jd)
        all_names.append(mosaic_path.name)

    LOG.info("Solving global additive background offsets across %d plates", len(candidates))
    bg_offsets = _solve_global_bg_offsets(
        all_reprojected,
        all_good,
        all_names,
    )

    # Pass 2: apply corrections and composite (oldest to newest; newest wins).
    mosaic = np.full(shape_out, np.nan, dtype=np.float32)
    epoch_map = np.full(shape_out, np.nan, dtype=np.float64)

    for i, (reprojected, good, obs_jd, offset) in enumerate(zip(all_reprojected, all_good, all_obs_jds, bg_offsets)):
        corrected = (reprojected.astype(np.float64) + offset).astype(np.float32)

        # Refine correction against already-placed neighbors using overlap only.
        # This directly targets seam agreement and damps over/under vignette subtraction.
        overlap = good & np.isfinite(mosaic)
        if config.subtract_background and i > 0 and int(np.count_nonzero(overlap)) >= 500:
            gain, add_offset, matched = _estimate_overlap_template_refinement(
                reference=mosaic,
                candidate=corrected,
                template=all_templates[i],
                overlap=overlap,
            )
            if matched:
                LOG.info(
                    "Overlap vignette refinement for %s: n=%d gain=%.4f offset=%.2f",
                    all_names[i],
                    int(np.count_nonzero(overlap)),
                    gain,
                    add_offset,
                )
                corrected = (
                    corrected.astype(np.float64)
                    + gain * all_templates[i].astype(np.float64)
                    + add_offset
                ).astype(np.float32)

        if i > 0 and int(np.count_nonzero(overlap)) >= 500:
            coeffs = _fit_overlap_residual_surface(
                reference=mosaic,
                candidate=corrected,
                overlap=overlap,
                degree=_RESIDUAL_SURFACE_DEGREE,
            )
            if coeffs is not None:
                surface = _evaluate_residual_surface(mosaic.shape, coeffs, degree=_RESIDUAL_SURFACE_DEGREE)
                corrected = (corrected.astype(np.float64) + surface.astype(np.float64)).astype(np.float32)
                LOG.info(
                    "Residual surface correction for %s: degree=%d",
                    all_names[i],
                    _RESIDUAL_SURFACE_DEGREE,
                )

        mosaic[good] = corrected[good]
        epoch_map[good] = obs_jd

    header = output_wcs.to_header(relax=True)
    header["BUNIT"] = "relative"
    header["DASHASOF"] = (jd_to_iso(config.as_of_jd) if config.as_of_jd else "open", "Most-recent as-of date")
    header["DASCHEA"] = (jd_to_iso(config.earliest_jd) if config.earliest_jd else "open", "Earliest allowed obs date")
    header["DASCHNPL"] = (len(candidates), "Number of plates used")
    header["DASCHBIN"] = (config.binning, "DASCH mosaic binning level")
    header["DASCHPSC"] = (pixel_scale_arcsec, "Output pixel scale arcsec/pixel")
    header["DASCHQST"] = (config.query_step_deg, "Grid spacing deg used for plate discovery")
    header["DASCHBG"] = (int(config.subtract_background), "Per-plate background subtracted before stitch")
    header["DASCHRS"] = (1, "Residual overlap surface correction enabled")
    header["DASCHRDG"] = (_RESIDUAL_SURFACE_DEGREE, "Residual overlap surface polynomial degree")
    header.add_history("Built from DASCH value-added mosaics fetched through daschlab.")
    header.add_history("Each pixel shows the most recently acquired plate covering that point.")

    _write_fits(config.output_fits, mosaic, header, overwrite=config.overwrite)

    if config.epoch_fits is not None:
        epoch_header = header.copy()
        epoch_header["BUNIT"] = "JD"
        epoch_header.add_comment("Pixel values are the Julian Date of the plate used at each position.")
        _write_fits(config.epoch_fits, epoch_map.astype(np.float32), epoch_header, overwrite=config.overwrite)

    if config.delete_base_mosaics:
        base_dir = config.session_root / "base_mosaics"
        if base_dir.is_dir():
            shutil.rmtree(base_dir)
            LOG.info("Deleted base mosaics from %s to free disk space", base_dir)

    manifest: dict[str, Any] = {
        "region": asdict(config.region),
        "as_of_date": jd_to_iso(config.as_of_jd) if config.as_of_jd else None,
        "earliest_date": jd_to_iso(config.earliest_jd) if config.earliest_jd else None,
        "projection": config.projection,
        "pixel_scale_arcsec": pixel_scale_arcsec,
        "binning": config.binning,
        "query_step_deg": config.query_step_deg,
        "subtract_background": config.subtract_background,
        "residual_surface_enabled": True,
        "residual_surface_degree": _RESIDUAL_SURFACE_DEGREE,
        "allow_multi_solution_plates": config.allow_multi_solution_plates,
        "plates": [
            {
                "plate_id": c.plate_id,
                "obs_date": jd_to_iso(c.obs_date_jd),
                "obs_date_jd": c.obs_date_jd,
                "n_wcs_solutions": c.n_wcs_solutions,
                "selected_at_points": c.selected_at_points,
                "local_mosaic_path": os.fspath(mosaic_paths[c.plate_id]),
            }
            for c in candidates
        ],
    }

    _prepare_output_path(config.manifest_json, overwrite=config.overwrite)
    config.manifest_json.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest