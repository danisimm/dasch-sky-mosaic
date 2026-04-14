from __future__ import annotations

import json
import logging
import math
import os
import shutil
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from daschlab import open_session
from reproject import reproject_interp
from scipy.ndimage import binary_erosion

from dasch_sky_mosaic.api import StarglassClient, parse_obs_date_jd

LOG = logging.getLogger(__name__)

_EDGE_TRIM_FRACTION = 0.05
_MAD_TO_STDDEV = 1.4826


def _normalize_ra_deg(value: float) -> float:
    return value % 360.0


def jd_to_iso(jd: float) -> str:
    return Time(jd, format="jd", scale="utc").isot # type: ignore


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


def _grid_axis(half_width_deg: float, spacing_deg: float) -> list[float]:
    if half_width_deg == 0:
        return [0.0]

    n_steps = max(1, math.ceil((2.0 * half_width_deg) / spacing_deg))
    axis = np.linspace(-half_width_deg, half_width_deg, num=n_steps + 1)
    return [float(v) for v in axis]


def _should_log_progress(index: int, total: int) -> bool:
    """Log each step for small totals; ~10 checkpoints for larger runs."""
    if total <= 20:
        return True
    stride = max(1, total // 10)
    return index == 1 or index == total or (index % stride == 0)


def _plate_interior_mask(image: np.ndarray, binning: int) -> np.ndarray:
    """Build support mask using geometric trimming of raw mosaic edges.

    We keep all finite data values and remove 5% from each edge of the raw
    plate mosaic (as requested) to avoid boundary artefacts.
    """
    support = np.isfinite(image)

    ny, nx = image.shape
    trim_y = max(1, int(round(_EDGE_TRIM_FRACTION * ny)))
    trim_x = max(1, int(round(_EDGE_TRIM_FRACTION * nx)))
    support[:trim_y, :] = False
    support[-trim_y:, :] = False
    support[:, :trim_x] = False
    support[:, -trim_x:] = False

    # A light erosion avoids one-pixel edge interpolation artefacts.
    erosion_iters = max(1, int(round(6.0 / max(1, binning))))
    return binary_erosion(support, iterations=erosion_iters)


def _robust_clip(arr: np.ndarray, lo_pct: float = 5.0, hi_pct: float = 95.0) -> np.ndarray:
    """Return a boolean mask selecting values within [lo_pct, hi_pct] percentile range."""
    lo, hi = np.nanpercentile(arr, [lo_pct, hi_pct])
    return (arr >= lo) & (arr <= hi)


def _fit_plate_background(image: np.ndarray, degree: int = 2) -> np.ndarray:
    """Fit and return a smooth 2D polynomial background model for a plate.

    Samples the plate on a coarse grid using the 20th-percentile of each block
    so that point sources (stars) are excluded and only the diffuse background
    is modelled.  The resulting surface captures both the large-scale sky
    background and radially symmetric vignetting (bright centre -> dark corners
    in photographic plates), which is a degree-2 radial polynomial.

    This is analogous to the per-image background modelling step in the IPAC
    Montage pipeline (mBackground), which must precede global offset matching so
    that spatially varying illumination does not bias the pairwise statistics.

    Returns an array of the same shape as `image`, dtype float32.
    """
    ny, nx = image.shape
    stride = max(4, min(ny, nx) // 40)
    half = stride // 2

    sample_y: list[float] = []
    sample_x: list[float] = []
    sample_v: list[float] = []

    for y in range(half, ny, stride):
        for x in range(half, nx, stride):
            block = image[max(0, y - half):min(ny, y + half),
                          max(0, x - half):min(nx, x + half)].ravel()
            finite = block[np.isfinite(block)]
            if len(finite) < 4:
                continue
            # Low percentile -> background, not stellar peaks.
            sample_y.append(y / ny - 0.5)   # centred on 0 for numerical stability
            sample_x.append(x / nx - 0.5)
            sample_v.append(float(np.percentile(finite, 20)))

    n_terms = (degree + 1) * (degree + 2) // 2
    if len(sample_v) < n_terms + 1:
        return np.full(image.shape, float(np.nanmedian(image)), dtype=np.float32)

    sy = np.array(sample_y)
    sx = np.array(sample_x)
    sv = np.array(sample_v, dtype=np.float64)

    # Sigma-clip blocks that are dominated by emission (bright galaxies, etc.).
    med = np.median(sv)
    mad = np.median(np.abs(sv - med)) * _MAD_TO_STDDEV + 1e-6
    keep = np.abs(sv - med) <= 5.0 * mad
    if keep.sum() < n_terms + 1:
        keep = np.ones(len(sv), dtype=bool)
    sy, sx, sv = sy[keep], sx[keep], sv[keep]

    def _design(y_arr: np.ndarray, x_arr: np.ndarray) -> np.ndarray:
        cols = []
        for i in range(degree + 1):
            for j in range(degree + 1 - i):
                cols.append(y_arr ** i * x_arr ** j)
        return np.column_stack(cols)

    try:
        coeffs, _, _, _ = np.linalg.lstsq(_design(sy, sx), sv, rcond=None)
    except Exception:
        return np.full(image.shape, float(np.nanmedian(image)), dtype=np.float32)

    Y, X = np.mgrid[0:ny, 0:nx]
    Yn = (Y / ny - 0.5).ravel().astype(np.float64)
    Xn = (X / nx - 0.5).ravel().astype(np.float64)
    bg = (_design(Yn, Xn) @ coeffs).reshape(ny, nx).astype(np.float32)
    return bg


def _solve_global_bg_offsets(
    reprojected_plates: list[np.ndarray],
    good_masks: list[np.ndarray],
    plate_names: list[str],
) -> list[float]:
    """Solve for per-plate additive background corrections to eliminate seam lines.

    Implements the global background-matching algorithm used by IPAC Montage
    (Berriman et al. 2003) and conceptually related to SDSS ubercalibration
    (Padmanabhan et al. 2008).

    Rather than matching each plate sequentially against the accumulated mosaic
    (which drifts with plate order), this collects all pairwise overlap
    measurements simultaneously and solves the overdetermined linear system

        a_i - a_j = d_ij   for all overlapping pairs (i, j)

    by least squares, where d_ij is the robust median of (plate_i - plate_j)
    in the shared footprint.  The most-recently observed plate (last in sorted
    order) is anchored at a = 0, so its background level becomes the
    photometric reference for the whole mosaic.

    Returns a list of additive offsets; add offset[i] to plate i before
    compositing.
    """
    N = len(reprojected_plates)
    if N <= 1:
        return [0.0] * N

    constraint_rows: list[np.ndarray] = []
    constraint_rhs: list[float] = []

    for i in range(N):
        for j in range(i + 1, N):
            overlap = good_masks[i] & good_masks[j]
            n_overlap = int(np.count_nonzero(overlap))
            if n_overlap < 300:
                continue

            diff = (
                reprojected_plates[i][overlap].astype(np.float64)
                - reprojected_plates[j][overlap].astype(np.float64)
            )
            core = diff[_robust_clip(diff)]
            if len(core) < 50:
                continue

            d_ij = float(np.median(core))
            row = np.zeros(N, dtype=np.float64)
            row[i] = 1.0
            row[j] = -1.0
            constraint_rows.append(row)
            constraint_rhs.append(d_ij)
            LOG.info(
                "Overlap constraint (%s, %s): n_pixels=%d  delta=%.2f",
                plate_names[i],
                plate_names[j],
                n_overlap,
                d_ij,
            )

    if not constraint_rows:
        LOG.warning(
            "No overlapping plate pairs with sufficient coverage; "
            "global background matching has no effect."
        )
        return [0.0] * N

    # Anchor the most-recently observed plate (last in date-sorted list).
    anchor = np.zeros(N, dtype=np.float64)
    anchor[N - 1] = 1.0
    constraint_rows.append(anchor)
    constraint_rhs.append(0.0)

    A = np.array(constraint_rows, dtype=np.float64)
    b = np.array(constraint_rhs, dtype=np.float64)
    offsets_raw, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    offsets = [float(x) for x in offsets_raw]
    LOG.info(
        "Global background offsets (reference=%s): %s",
        plate_names[N - 1],
        "  ".join(f"{plate_names[i]}={offsets[i]:+.2f}" for i in range(N)),
    )
    return offsets


def _estimate_overlap_template_refinement(
    reference: np.ndarray,
    candidate: np.ndarray,
    template: np.ndarray,
    overlap: np.ndarray,
) -> tuple[float, float, bool]:
    """Estimate seam-focused correction: candidate += gain*template + offset.

    This uses only overlap pixels, so the adjustment is driven by agreement
    where plates actually meet (not by interior brightness structure).
    """
    n_overlap = int(np.count_nonzero(overlap))
    if n_overlap < 500:
        return 0.0, 0.0, False

    y = (reference[overlap] - candidate[overlap]).astype(np.float64)
    t = template[overlap].astype(np.float64)

    keep = _robust_clip(y) & _robust_clip(t)
    if int(np.count_nonzero(keep)) < 300:
        return 0.0, 0.0, False
    y_fit = y[keep]
    t_fit = t[keep]

    if np.nanstd(t_fit) < 1e-3:
        return 0.0, float(np.nanmedian(y_fit)), True

    X = np.column_stack((t_fit, np.ones_like(t_fit)))
    try:
        gain, offset = np.linalg.lstsq(X, y_fit, rcond=None)[0]
    except Exception:
        return 0.0, 0.0, False

    # Constrain to plausible corrections; reject unstable fits.
    if not np.isfinite(gain) or abs(gain) > 0.6:
        return 0.0, 0.0, False
    if not np.isfinite(offset) or abs(offset) > 1500:
        return 0.0, 0.0, False
    return float(gain), float(offset), True


def iter_query_points(region: Region, query_step_deg: float) -> list[tuple[float, float]]:
    dec_offsets = _grid_axis(region.height_deg / 2.0, query_step_deg)
    points: list[tuple[float, float]] = []

    for dec_offset in dec_offsets:
        dec_deg = max(-90.0, min(90.0, region.dec_deg + dec_offset))
        cos_dec = max(math.cos(math.radians(dec_deg)), 0.05)
        ra_spacing = query_step_deg / cos_dec
        ra_offsets = _grid_axis(region.width_deg / 2.0, ra_spacing)

        for ra_offset in ra_offsets:
            points.append((_normalize_ra_deg(region.ra_deg + ra_offset), dec_deg))

    unique: dict[tuple[int, int], tuple[float, float]] = {}
    for ra_deg, dec_deg in points:
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

    plate_exposures: dict[str, list[tuple[float | None, int | None]]] = {}
    point_winner: dict[tuple[int, int], str] = {}

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
            plate_exposures.setdefault(h.plate_id, []).append((h.obs_date_jd, h.solnum))
        best = max(
            eligible,
            key=lambda h: h.obs_date_jd if h.obs_date_jd is not None else -1.0,
            default=None,
        )
        if best is not None:
            point_winner[(round(ra_deg * 1000), round(dec_deg * 1000))] = best.plate_id

    selection_counts = Counter(point_winner.values())
    LOG.info(
        "Discovered %d candidate plates; %d selected as most-recent at >= 1 query point",
        len(plate_exposures),
        len(selection_counts),
    )

    candidates: list[CandidatePlate] = []
    for plate_id, point_count in selection_counts.items():
        exposures = plate_exposures.get(plate_id, [])
        solnums = {s for _, s in exposures if s is not None}
        n_wcs = len(solnums)
        if not config.allow_multi_solution_plates and n_wcs > 1:
            continue
        valid_jds = [jd for jd, _ in exposures if jd is not None]
        if not valid_jds:
            continue
        candidates.append(CandidatePlate(
            plate_id=plate_id,
            obs_date_jd=max(valid_jds),
            n_wcs_solutions=n_wcs,
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


def _download_mosaic_paths(session_root: Path, plate_ids: list[str], binning: int) -> dict[str, Path]:
    session_root.mkdir(parents=True, exist_ok=True)
    session = open_session(str(session_root))
    paths: dict[str, Path] = {}

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


def _write_fits(path: Path, data: np.ndarray, header: fits.Header, overwrite: bool) -> None:
    _prepare_output_path(path, overwrite)
    fits.PrimaryHDU(data=data.astype(np.float32), header=header).writeto(path, overwrite=overwrite)


def build_mosaic(config: BuildConfig) -> dict[str, Any]:
    candidates = discover_candidate_plates(config)
    if not candidates:
        raise RuntimeError("no DASCH plates matched the requested sky region and date constraints")

    mosaic_paths = _download_mosaic_paths(
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

    # Solve global per-plate additive background offsets (Montage/ubercal approach).
    # This prevents the drift that accumulates when plates are matched sequentially
    # against the growing mosaic: instead every pairwise constraint is used at once.
    LOG.info("Solving global background offsets across %d plates", len(candidates))
    bg_offsets = _solve_global_bg_offsets(all_reprojected, all_good, all_names)

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