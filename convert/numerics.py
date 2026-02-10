from __future__ import annotations

from itertools import permutations

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline


def _sorted_xy(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    x_sorted = x[order]
    y_sorted = y[order]

    # Keep first occurrence for duplicate coordinates.
    keep = np.ones_like(x_sorted, dtype=bool)
    keep[1:] = np.diff(x_sorted) > 0.0
    return x_sorted[keep], y_sorted[keep]


def interp1_linear(x: np.ndarray, y: np.ndarray, xq: np.ndarray | float) -> np.ndarray | float:
    """MATLAB-like interp1(x, y, xq) linear interpolation without extrapolation."""
    x_sorted, y_sorted = _sorted_xy(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    if x_sorted.size < 2:
        raise ValueError("Need at least two points for linear interpolation")

    xq_arr = np.asarray(xq, dtype=float)
    xq_flat = np.atleast_1d(xq_arr)

    out = np.interp(xq_flat, x_sorted, y_sorted)
    outside = (xq_flat < x_sorted[0]) | (xq_flat > x_sorted[-1])
    out[outside] = np.nan

    if xq_arr.ndim == 0:
        return float(out[0])
    return out.reshape(xq_arr.shape)


def interp1_spline(x: np.ndarray, y: np.ndarray, xq: np.ndarray | float) -> np.ndarray | float:
    """MATLAB-like spline(x, y, xq) wrapper."""
    x_sorted, y_sorted = _sorted_xy(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    if x_sorted.size < 2:
        raise ValueError("Need at least two points for spline interpolation")

    if x_sorted.size == 2:
        return interp1_linear(x_sorted, y_sorted, xq)

    spline = CubicSpline(x_sorted, y_sorted, bc_type="not-a-knot", extrapolate=True)
    xq_arr = np.asarray(xq, dtype=float)
    out = np.asarray(spline(xq_arr), dtype=float)
    outside = (xq_arr < x_sorted[0]) | (xq_arr > x_sorted[-1])
    out[outside] = np.nan
    if xq_arr.ndim == 0:
        return float(out)
    return out


def interp2_spline(
    x: np.ndarray,
    y: np.ndarray,
    v_t: np.ndarray,
    xq: float,
    yq: float,
) -> float:
    """
    MATLAB-like interp2(x, y, V, xq, yq, 'spline').

    Expected layout matches MATLAB vector-grid form:
    - len(x) is number of columns in V
    - len(y) is number of rows in V
    - V is passed as v_t with shape (len(y), len(x))
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    v_t = np.asarray(v_t, dtype=float)

    if v_t.shape != (y.size, x.size):
        raise ValueError(
            f"interp2_spline: shape mismatch, got V={v_t.shape}, expected {(y.size, x.size)}"
        )

    kx = min(3, y.size - 1)
    ky = min(3, x.size - 1)
    if xq < x.min() or xq > x.max() or yq < y.min() or yq > y.max():
        return float("nan")

    spline = RectBivariateSpline(y, x, v_t, kx=kx, ky=ky)
    return float(spline.ev(yq, xq))


def coerce_apar_shape(apar: np.ndarray, nx: int, ny: int, nz: int) -> np.ndarray:
    """Try to reorder a loaded MAT array into (nx, ny, nz)."""
    apar = np.asarray(apar)
    if apar.ndim != 3:
        raise ValueError(f"Expected 3D apar array, got ndim={apar.ndim}")

    target = (nx, ny, nz)
    if apar.shape == target:
        return apar.astype(float, copy=False)

    for perm in permutations((0, 1, 2)):
        candidate = np.transpose(apar, perm)
        if candidate.shape == target:
            return candidate.astype(float, copy=False)

    raise ValueError(f"Could not coerce apar shape {apar.shape} into {target}")


def clamp_index_1based(value: float, max_value: int) -> int:
    idx = int(np.rint(value))
    return max(1, min(max_value, idx))
