# This file was generated with the assistance of GitHub Copilot (AI).
from __future__ import annotations

import logging

import numpy as np
from scipy.ndimage import binary_erosion

LOG = logging.getLogger(__name__)

_EDGE_TRIM_FRACTION = 0.05
_MAD_TO_STDDEV = 1.4826


def _robust_clip(arr: np.ndarray, lo_pct: float = 5.0, hi_pct: float = 95.0) -> np.ndarray:
    """Return a boolean mask selecting values within [lo_pct, hi_pct] percentile range."""
    lo, hi = np.nanpercentile(arr, [lo_pct, hi_pct])
    return (arr >= lo) & (arr <= hi)


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
            y0, y1 = max(0, y - half), min(ny, y + half)
            x0, x1 = max(0, x - half), min(nx, x + half)
            block = image[y0:y1, x0:x1]
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
                "Overlap constraint (%s, %s): n_pixels=%d delta=%.2f",
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


def _fit_overlap_residual_surface(
    reference: np.ndarray,
    candidate: np.ndarray,
    overlap: np.ndarray,
    degree: int = 1,
) -> np.ndarray | None:
    """Fit a low-order polynomial residual surface on overlap pixels."""
    n_overlap = int(np.count_nonzero(overlap))
    if n_overlap < 500:
        return None

    y_idx, x_idx = np.where(overlap)
    delta = (reference[overlap] - candidate[overlap]).astype(np.float64)
    keep = _robust_clip(delta)
    if int(np.count_nonzero(keep)) < 300:
        return None

    y = y_idx[keep].astype(np.float64)
    x = x_idx[keep].astype(np.float64)
    d = delta[keep]

    ny, nx = reference.shape
    yn = y / max(1.0, float(ny - 1)) - 0.5
    xn = x / max(1.0, float(nx - 1)) - 0.5

    if degree == 1:
        X = np.column_stack((np.ones_like(yn), yn, xn))
    elif degree == 2:
        X = np.column_stack((np.ones_like(yn), yn, xn, yn * yn, yn * xn, xn * xn))
    else:
        return None

    try:
        coeffs, _, _, _ = np.linalg.lstsq(X, d, rcond=None)
    except Exception:
        return None
    if not np.all(np.isfinite(coeffs)):
        return None

    if np.max(np.abs(coeffs)) > 2000.0:
        return None
    return coeffs.astype(np.float64)


def _evaluate_residual_surface(shape: tuple[int, int], coeffs: np.ndarray, degree: int = 1) -> np.ndarray:
    """Evaluate residual surface coefficients on a full image grid."""
    ny, nx = shape
    y, x = np.mgrid[0:ny, 0:nx]
    yn = y.astype(np.float64) / max(1.0, float(ny - 1)) - 0.5
    xn = x.astype(np.float64) / max(1.0, float(nx - 1)) - 0.5

    if degree == 1:
        out = coeffs[0] + coeffs[1] * yn + coeffs[2] * xn
    else:
        out = (
            coeffs[0]
            + coeffs[1] * yn
            + coeffs[2] * xn
            + coeffs[3] * yn * yn
            + coeffs[4] * yn * xn
            + coeffs[5] * xn * xn
        )
    return out.astype(np.float32)
