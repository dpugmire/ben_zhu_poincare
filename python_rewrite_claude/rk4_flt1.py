"""
Field-line tracing with RK4 integration.

RK4 method: for dy/dx = f
    k1 = h*f(xn, yn)
    k2 = h*f(xn + h/2, yn + k1/2)
    k3 = h*f(xn + h/2, yn + k2/2)
    k4 = h*f(xn + h, yn + k3)
    yn+1 = yn + (k1 + 2*k2 + 2*k3 + k4) / 6

Translated from MATLAB (B. Zhu) to Python.

Indexing convention:
    All index variables (yStart, nypf1, nypf2) use 1-based values matching
    the MATLAB original. Coordinate arrays (xarray, zarray) are physical
    values. Array access into Python arrays uses 0-based indexing (with -1
    offset where needed).
"""

import numpy as np
from scipy.interpolate import RectBivariateSpline


def rk4_flt1(xStart, yStart, zStart, dxdy, dzdy, xarray, zarray,
             region, dxdy_pm1, dzdy_pm1, direction, nypf1, nypf2):
    """
    RK4 field-line tracing step.

    Parameters
    ----------
    xStart : float
        Starting x position (psi value).
    yStart : int
        Starting y index (1-based).
    zStart : float
        Starting z position (toroidal angle).
    dxdy : ndarray, shape (nx, ny, nzG)
        dx/dy field on the 3D grid.
    dzdy : ndarray, shape (nx, ny, nzG)
        dz/dy field on the 3D grid.
    xarray : ndarray, shape (nx,)
        Physical x (psi) coordinate values.
    zarray : ndarray, shape (nzG+1,)
        Toroidal z values including wrapped endpoint [0, dz, ..., 2*pi].
    region : int
        Region identifier (0=CFR, 1=SOL, 2=PFR).
    dxdy_pm1 : ndarray, shape (nx, nzG)
        dx/dy at the twist-shift boundary (plus-1 or minus-1).
    dzdy_pm1 : ndarray, shape (nx, nzG)
        dz/dy at the twist-shift boundary (plus-1 or minus-1).
    direction : int
        1 for y-increasing, -1 for y-decreasing.
    nypf1 : int
        Lower private-flux boundary y-index (1-based).
    nypf2 : int
        Upper private-flux boundary y-index (1-based).

    Returns
    -------
    xEnd : float
        Ending x position (psi value).
    zEnd : float
        Ending z position (toroidal angle).
    """
    # Step size
    h = 1.0
    hh = h / 2.0
    h6 = h / 6.0

    # Need half-step and full-step info
    # Note: yStart is 1-based; Python array index = yStart - 1
    if direction == 1:
        dxdyp = dxdy[:, yStart - 1, :].copy()
        dzdyp = dzdy[:, yStart - 1, :].copy()

        if region == 0 and yStart == nypf2:
            dxdyn = dxdy_pm1.copy()
            dxdyh = 0.5 * (dxdyp + dxdyn)
            dzdyn = dzdy_pm1.copy()
            dzdyh = 0.5 * (dzdyp + dzdyn)
        else:
            # yStart+1 in MATLAB (1-based) = index yStart in Python (0-based)
            dxdyn = dxdy[:, yStart, :].copy()
            dxdyh = 0.5 * (dxdy[:, yStart - 1, :] + dxdy[:, yStart, :])
            dzdyn = dzdy[:, yStart, :].copy()
            dzdyh = 0.5 * (dzdy[:, yStart - 1, :] + dzdy[:, yStart, :])

    elif direction == -1:
        dxdyp = dxdy[:, yStart - 1, :].copy()
        dzdyp = dzdy[:, yStart - 1, :].copy()

        if region == 0 and yStart == nypf1 + 1:
            dxdyn = dxdy_pm1.copy()
            dxdyh = 0.5 * (dxdyp + dxdyn)
            dzdyn = dzdy_pm1.copy()
            dzdyh = 0.5 * (dzdyp + dzdyn)
        else:
            # yStart-1 in MATLAB (1-based) = index yStart-2 in Python (0-based)
            dxdyn = dxdy[:, yStart - 2, :].copy()
            dxdyh = 0.5 * (dxdy[:, yStart - 1, :] + dxdy[:, yStart - 2, :])
            dzdyn = dzdy[:, yStart - 2, :].copy()
            dzdyh = 0.5 * (dzdy[:, yStart - 1, :] + dzdy[:, yStart - 2, :])
    else:
        raise ValueError("Check direction parameter setting!")

    # Patch last z point (periodicity padding)
    dxdyp = np.column_stack([dxdyp, dxdyp[:, -1]])
    dxdyn = np.column_stack([dxdyn, dxdyn[:, -1]])
    dxdyh = np.column_stack([dxdyh, dxdyh[:, -1]])
    dzdyp = np.column_stack([dzdyp, dzdyp[:, -1]])
    dzdyn = np.column_stack([dzdyn, dzdyn[:, -1]])
    dzdyh = np.column_stack([dzdyh, dzdyh[:, -1]])

    # Build 2D spline interpolators: axes are (xarray, zarray)
    # Each data array has shape (nx, nzG+1)
    spl_dxdyp = RectBivariateSpline(xarray, zarray, dxdyp)
    spl_dzdyp = RectBivariateSpline(xarray, zarray, dzdyp)
    spl_dxdyh = RectBivariateSpline(xarray, zarray, dxdyh)
    spl_dzdyh = RectBivariateSpline(xarray, zarray, dzdyh)
    spl_dxdyn = RectBivariateSpline(xarray, zarray, dxdyn)
    spl_dzdyn = RectBivariateSpline(xarray, zarray, dzdyn)

    two_pi = 2.0 * np.pi

    # First RK4 step
    dxdy1 = float(spl_dxdyp.ev(xStart, zStart))
    dzdy1 = float(spl_dzdyp.ev(xStart, zStart))
    x1 = xStart + direction * hh * dxdy1
    z1 = zStart + direction * hh * dzdy1

    # Second RK4 step
    dxdy2 = float(spl_dxdyh.ev(x1, z1 % two_pi))
    dzdy2 = float(spl_dzdyh.ev(x1, z1 % two_pi))
    x2 = xStart + direction * hh * dxdy2
    z2 = zStart + direction * hh * dzdy2

    # Third RK4 step
    dxdy3 = float(spl_dxdyh.ev(x2, z2 % two_pi))
    dzdy3 = float(spl_dzdyh.ev(x2, z2 % two_pi))
    x3 = xStart + direction * dxdy3
    z3 = zStart + direction * dzdy3

    # Fourth RK4 step
    dxdy4 = float(spl_dxdyn.ev(x3, z3 % two_pi))
    dzdy4 = float(spl_dzdyn.ev(x3, z3 % two_pi))

    # Accumulate increments with proper weights
    xEnd = xStart + direction * h6 * (dxdy1 + 2 * dxdy2 + 2 * dxdy3 + dxdy4)
    zEnd = zStart + direction * h6 * (dzdy1 + 2 * dzdy2 + 2 * dzdy3 + dzdy4)

    return xEnd, zEnd
