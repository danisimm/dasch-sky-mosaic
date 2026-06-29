"""Per-plate background estimation and global background matching.

Implements the Montage-style additive background matching algorithm:
  1. Fit a 2D polynomial background to each plate (per-plate).
  2. Solve for additive per-plate corrections using pairwise overlap residuals.
  3. Optionally refine on overlap seams using a vignette template.
  4. Optionally fit and apply a low-order residual surface per plate.
"""
from __future__ import annotations

import logging

import numpy as np
from scipy.ndimage import binary_erosion

LOG = logging.getLogger(__name__)

_EDGE_TRIM_FRACTION = 0.05
_MAD_TO_STDDEV = 1.4826


def _robust_clip(arr: np.ndarray, lo_pct: float = 5.0, hi_pct: float = 95.0) -> np.ndarray:
    lo, hi = np.nanpercentile(arr, [lo_pct, hi_pct])
    return (arr >= lo) & (arr <= hi)


def _plate_interior_mask(image: np.ndarray, binning: int) -> np.ndarray:
    """Build support mask: keep finite pixels, trim 5% from each edge, lightly erode."""
    support = np.isfinite(image)
    ny, nx = image.shape
    trim_y = max(1, int(round(_EDGE_TRIM_FRACTION * ny)))
    trim_x = max(1, int(round(_EDGE_TRIM_FRACTION * nx)))
    support[:trim_y, :] = False
    support[-trim_y:, :] = False
    support[:, :trim_x] = False
    support[:, -trim_x:] = False
    erosion_iters = max(1, int(round(6.0 / max(1, binning))))
    return binary_erosion(support, iterations=erosion_iters)


def _fit_plate_background(image: np.ndarray, degree: int = 2) -> np.ndarray:
    """Fit a smooth 2D polynomial background using low-percentile block sampling.

    Analogous to mBackground in the IPAC Montage pipeline. Returns a float32
    array of the same shape as `image`.
    """
    ny, nx = image.shape
    stride = max(4, min(ny, nx) // 40)
    half = stride // 2

    sample_y: list[float] = []
    sample_x: list[float] = []
    sample_v: list[float] = []

    for y in range(half, ny, stride):
        for x in range(half, nx, stride):
            block = image[max(0, y - half):min(ny, y + half), max(0, x - half):min(nx, x + half)]
            finite = block[np.isfinite(block)]
            if len(finite) < 4:
                continue
            sample_y.append(y / ny - 0.5)
            sample_x.append(x / nx - 0.5)
            sample_v.append(float(np.percentile(finite, 20)))

    n_terms = (degree + 1) * (degree + 2) // 2
    if len(sample_v) < n_terms + 1:
        return np.full(image.shape, float(np.nanmedian(image)), dtype=np.float32)

    sy = np.array(sample_y)
    sx = np.array(sample_x)
    sv = np.array(sample_v, dtype=np.float64)

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
    return (_design(Yn, Xn) @ coeffs).reshape(ny, nx).astype(np.float32)


def _solve_global_bg_offsets(
    reprojected_plates: list[np.ndarray],
    good_masks: list[np.ndarray],
    plate_names: list[str],
) -> list[float]:
    """Solve per-plate additive corrections by pairwise overlap least-squares.

    Implements the Montage / SDSS ubercalibration algorithm. The most-recently
    observed plate (last in sorted order) is anchored at offset 0.
    """
    N = len(reprojected_plates)
    if N <= 1:
        return [0.0] * N

    constraint_rows: list[np.ndarray] = []
    constraint_rhs: list[float] = []

    for i in range(N):
        for j in range(i + 1, N):
            overlap = good_masks[i] & good_masks[j]
            if int(np.count_nonzero(overlap)) < 300:
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
                "Overlap constraint (%s, %s): n=%d delta=%.2f",
                plate_names[i], plate_names[j], int(np.count_nonzero(overlap)), d_ij,
            )

    if not constraint_rows:
        LOG.warning("No overlapping pairs with sufficient coverage; background matching has no effect.")
        return [0.0] * N

    anchor = np.zeros(N, dtype=np.float64)
    anchor[N - 1] = 1.0
    constraint_rows.append(anchor)
    constraint_rhs.append(0.0)

    A = np.array(constraint_rows, dtype=np.float64)
    b = np.array(constraint_rhs, dtype=np.float64)
    offsets_raw, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    offsets = [float(x) for x in offsets_raw]
    LOG.info(
        "Global background offsets (ref=%s): %s",
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
    """Fit candidate += gain*template + offset using overlap pixels only."""
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
    if int(np.count_nonzero(overlap)) < 500:
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
    if not np.all(np.isfinite(coeffs)) or np.max(np.abs(coeffs)) > 2000.0:
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
            coeffs[0] + coeffs[1] * yn + coeffs[2] * xn
            + coeffs[3] * yn * yn + coeffs[4] * yn * xn + coeffs[5] * xn * xn
        )
    return out.astype(np.float32)
