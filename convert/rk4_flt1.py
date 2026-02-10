from __future__ import annotations

import numpy as np

try:
    from .numerics import interp2_spline
except ImportError:
    from numerics import interp2_spline

TWOPI = 2.0 * np.pi


def rk4_flt1(
    x_start: float,
    y_start: float,
    z_start: float,
    dxdy: np.ndarray,
    dzdy: np.ndarray,
    xarray: np.ndarray,
    zarray: np.ndarray,
    region: int,
    dxdy_pm1: np.ndarray,
    dzdy_pm1: np.ndarray,
    direction: int,
    nypf1: int,
    nypf2: int,
) -> tuple[float, float]:
    """Field-line tracing with RK4 integration translated from RK4_FLT1.m."""
    hh = 0.5
    h6 = 1.0 / 6.0

    y_idx = int(np.rint(y_start)) - 1

    if direction == 1:
        dxdyp = np.asarray(dxdy[:, y_idx, :], dtype=float)
        dzdyp = np.asarray(dzdy[:, y_idx, :], dtype=float)
        if region == 0 and int(np.rint(y_start)) == nypf2:
            dxdyn = np.asarray(dxdy_pm1, dtype=float)
            dzdyn = np.asarray(dzdy_pm1, dtype=float)
            dxdyh = 0.5 * (dxdyp + dxdyn)
            dzdyh = 0.5 * (dzdyp + dzdyn)
        else:
            dxdyn = np.asarray(dxdy[:, y_idx + 1, :], dtype=float)
            dzdyn = np.asarray(dzdy[:, y_idx + 1, :], dtype=float)
            dxdyh = 0.5 * (dxdyp + dxdyn)
            dzdyh = 0.5 * (dzdyp + dzdyn)
    elif direction == -1:
        dxdyp = np.asarray(dxdy[:, y_idx, :], dtype=float)
        dzdyp = np.asarray(dzdy[:, y_idx, :], dtype=float)
        if region == 0 and int(np.rint(y_start)) == nypf1 + 1:
            dxdyn = np.asarray(dxdy_pm1, dtype=float)
            dzdyn = np.asarray(dzdy_pm1, dtype=float)
            dxdyh = 0.5 * (dxdyp + dxdyn)
            dzdyh = 0.5 * (dzdyp + dzdyn)
        else:
            dxdyn = np.asarray(dxdy[:, y_idx - 1, :], dtype=float)
            dzdyn = np.asarray(dzdy[:, y_idx - 1, :], dtype=float)
            dxdyh = 0.5 * (dxdyp + dxdyn)
            dzdyh = 0.5 * (dzdyp + dzdyn)
    else:
        raise ValueError("direction must be +1 or -1")

    # Patch one additional z point (MATLAB: arr(:,end+1)=arr(:,end)).
    dxdyp = np.concatenate([dxdyp, dxdyp[:, -1:]], axis=1)
    dxdyn = np.concatenate([dxdyn, dxdyn[:, -1:]], axis=1)
    dxdyh = np.concatenate([dxdyh, dxdyh[:, -1:]], axis=1)

    dzdyp = np.concatenate([dzdyp, dzdyp[:, -1:]], axis=1)
    dzdyn = np.concatenate([dzdyn, dzdyn[:, -1:]], axis=1)
    dzdyh = np.concatenate([dzdyh, dzdyh[:, -1:]], axis=1)

    dxdy1 = interp2_spline(xarray, zarray, dxdyp.T, x_start, z_start)
    dzdy1 = interp2_spline(xarray, zarray, dzdyp.T, x_start, z_start)
    x1 = x_start + direction * hh * dxdy1
    z1 = z_start + direction * hh * dzdy1

    dxdy2 = interp2_spline(xarray, zarray, dxdyh.T, x1, np.mod(z1, TWOPI))
    dzdy2 = interp2_spline(xarray, zarray, dzdyh.T, x1, np.mod(z1, TWOPI))
    x2 = x_start + direction * hh * dxdy2
    z2 = z_start + direction * hh * dzdy2

    dxdy3 = interp2_spline(xarray, zarray, dxdyh.T, x2, np.mod(z2, TWOPI))
    dzdy3 = interp2_spline(xarray, zarray, dzdyh.T, x2, np.mod(z2, TWOPI))
    x3 = x_start + direction * dxdy3
    z3 = z_start + direction * dzdy3

    dxdy4 = interp2_spline(xarray, zarray, dxdyn.T, x3, np.mod(z3, TWOPI))
    dzdy4 = interp2_spline(xarray, zarray, dzdyn.T, x3, np.mod(z3, TWOPI))

    x_end = x_start + direction * h6 * (dxdy1 + 2.0 * dxdy2 + 2.0 * dxdy3 + dxdy4)
    z_end = z_start + direction * h6 * (dzdy1 + 2.0 * dzdy2 + 2.0 * dzdy3 + dzdy4)
    return float(x_end), float(z_end)
