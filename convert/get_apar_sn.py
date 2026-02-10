from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

try:
    from .numerics import interp1_spline
except ImportError:
    from numerics import interp1_spline


def _estimate_dpsi(psixy: np.ndarray) -> float:
    nx, ny = psixy.shape
    i1 = min(100, nx - 1)
    i0 = max(0, i1 - 1)
    j = min(38, ny - 1)
    dpsi = float(psixy[i1, j] - psixy[i0, j])
    if dpsi == 0.0:
        col = psixy[:, j]
        diffs = np.diff(col)
        nz = diffs[np.nonzero(diffs)]
        if nz.size == 0:
            return 1.0
        dpsi = float(np.mean(np.abs(nz)))
    return abs(dpsi)


def _kz_modes(nz: int, zperiod: int) -> np.ndarray:
    return np.concatenate((np.arange(0, nz // 2 + 1), np.arange(-nz // 2 + 1, 0))) * zperiod


def get_apar_sn(
    psi: np.ndarray,
    bxy: np.ndarray,
    psixy: np.ndarray,
    zs: np.ndarray,
    sa: np.ndarray,
    sinty: np.ndarray,
    dy0: float,
    dz: float,
    ixsep: int,
    nypf1: int,
    nypf2: int,
    zperiod: int,
    interp_opt: int,
    deriv_opt: int,
    true_apar: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Prepare and calculate apar and derivatives for single-null grids."""
    psi = np.asarray(psi, dtype=float).copy()
    bxy = np.asarray(bxy, dtype=float)
    psixy = np.asarray(psixy, dtype=float)
    zs = np.asarray(zs, dtype=float)
    sa = np.asarray(sa, dtype=float)
    sinty = np.asarray(sinty, dtype=float)

    nx, ny, nz = psi.shape

    apar = np.zeros((nx, ny, nz), dtype=float)
    dapardx = np.zeros((nx, ny, nz), dtype=float)
    dapardy = np.zeros((nx, ny, nz), dtype=float)
    dapardz = np.zeros((nx, ny, nz), dtype=float)

    # MATLAB indexing reminder:
    # iy_cfr = [nypf1+1:nypf2] (1-based) -> zero-based [nypf1:nypf2-1]
    iy_cfr = np.arange(nypf1, nypf2)
    iy_pfr = np.concatenate([np.arange(0, nypf1), np.arange(nypf2, ny)])
    ny_cfr = nypf2 - nypf1 + 1

    zarray = np.arange(nz, dtype=float) * dz
    zarrayp = np.arange(nz + 1, dtype=float) * dz
    zmax = dz * nz

    if not true_apar:
        psis = psi.copy()
        psisp1 = np.zeros((ixsep, nz), dtype=float)

        for i in range(ixsep):
            zarray_shift = np.mod(zarray + sa[i], zmax)
            psis_tmp = np.concatenate([psi[i, nypf1, :], psi[i, nypf1, 0:1]])
            psisp1[i, :] = interp1_spline(zarrayp, psis_tmp, zarray_shift)

        if interp_opt == 0:
            for i in range(nx):
                if i < ixsep:
                    psi[i, iy_cfr[:-1], :] = 0.5 * (
                        psis[i, iy_cfr[:-1], :] + psis[i, iy_cfr[1:], :]
                    )
                    psi[i, iy_cfr[-1], :] = 0.5 * (psisp1[i, :] + psis[i, iy_cfr[-1], :])

                    psi[i, iy_pfr[:-1], :] = 0.5 * (
                        psis[i, iy_pfr[:-1], :] + psis[i, iy_pfr[1:], :]
                    )
                    psi[i, iy_pfr[-1], :] = 1.5 * psis[i, iy_pfr[-1], :] - 0.5 * psis[i, iy_pfr[-2], :]
                else:
                    psi[i, 0 : ny - 1, :] = 0.5 * (psis[i, 0 : ny - 1, :] + psis[i, 1:ny, :])
                    psi[i, ny - 1, :] = 1.5 * psis[i, ny - 1, :] - 0.5 * psis[i, ny - 2, :]

        elif interp_opt == 1:
            x_cfr = np.arange(1, ny_cfr + 1, dtype=float)
            xq_cfr = np.arange(1, ny_cfr, dtype=float) + 0.5
            x_pfr = np.arange(1, 2 * nypf1 + 2, dtype=float)
            xq_pfr = np.arange(1, 2 * nypf1 + 1, dtype=float) + 0.5
            x_sol = np.arange(1, ny + 2, dtype=float)
            xq_sol = np.arange(1, ny + 1, dtype=float) + 0.5

            for i in range(nx):
                if i < ixsep:
                    for k in range(nz):
                        psis_cfr = np.concatenate([psis[i, iy_cfr, k], [psisp1[i, k]]])
                        psi[i, iy_cfr, k] = interp1_spline(x_cfr, psis_cfr, xq_cfr)

                        psis_pfr = psis[i, iy_pfr, k]
                        pfr_ext = 2.0 * psis_pfr[-1] - psis_pfr[-2]
                        psis_pfr = np.concatenate([psis_pfr, [pfr_ext]])
                        psi[i, iy_pfr, k] = interp1_spline(x_pfr, psis_pfr, xq_pfr)
                else:
                    for k in range(nz):
                        psis_tmp = np.concatenate([psis[i, :, k], [psis[i, -1, k]]])
                        psi[i, :, k] = interp1_spline(x_sol, psis_tmp, xq_sol)
        else:
            raise ValueError("Unknown interpolation method for psi")

        apar = psi * bxy[:, :, None]
    else:
        apar = psi.copy()

    apars = np.zeros((nx, ny, nz), dtype=float)
    dapardpsi = np.zeros((nx, ny, nz), dtype=float)
    kz = _kz_modes(nz, zperiod)

    for i in range(nx):
        for j in range(ny):
            apars[i, j, :] = np.real(
                np.fft.ifft(np.fft.fft(apar[i, j, :]) * np.exp(-1j * zs[i, j] * kz))
            )

    if deriv_opt == 0:
        dpsi = _estimate_dpsi(psixy)
        dapardpsi[1 : nx - 1, :, :] = 0.5 * (apars[2:nx, :, :] - apars[0 : nx - 2, :, :]) / dpsi

        dapardz[:, :, 0] = 0.5 * (apar[:, :, 1] - apar[:, :, nz - 1]) / dz
        dapardz[:, :, 1 : nz - 1] = 0.5 * (apar[:, :, 2:nz] - apar[:, :, 0 : nz - 2]) / dz
        dapardz[:, :, nz - 1] = 0.5 * (apar[:, :, 0] - apar[:, :, nz - 2]) / dz

        for k in range(nz):
            dapardx[:, :, k] = dapardpsi[:, :, k] + sinty * dapardz[:, :, k]

        for i in range(nx):
            if i < ixsep:
                dapardy[i, iy_cfr[1:-1], :] = 0.5 * (
                    apar[i, iy_cfr[2:], :] - apar[i, iy_cfr[:-2], :]
                ) / dy0

                zarray_shift = np.mod(zarray - sa[i], zmax)
                apar_tmp = np.concatenate([apar[i, nypf2 - 1, :], apar[i, nypf2 - 1, 0:1]])
                aparm1 = interp1_spline(zarrayp, apar_tmp, zarray_shift)
                if nypf1 + 1 < ny:
                    dapardy[i, nypf1, :] = 0.5 * (apar[i, nypf1 + 1, :] - aparm1) / dy0

                zarray_shift = np.mod(zarray + sa[i], zmax)
                apar_tmp = np.concatenate([apar[i, nypf1, :], apar[i, nypf1, 0:1]])
                aparp1 = interp1_spline(zarrayp, apar_tmp, zarray_shift)
                dapardy[i, nypf2 - 1, :] = 0.5 * (aparp1 - apar[i, nypf2 - 2, :]) / dy0

                dapardy[i, iy_pfr[1:-1], :] = 0.5 * (
                    apar[i, iy_pfr[2:], :] - apar[i, iy_pfr[:-2], :]
                ) / dy0
            else:
                dapardy[i, 1 : ny - 1, :] = 0.5 * (apar[i, 2:ny, :] - apar[i, 0 : ny - 2, :]) / dy0

    elif deriv_opt == 1:
        for k in range(nz):
            for j in range(ny):
                apar_tmp = np.concatenate([apars[:, j, k], [apars[-1, j, k]]])
                xarray = np.concatenate([psixy[:, j], [2.0 * psixy[-1, j] - psixy[-2, j]]])
                cs = CubicSpline(xarray, apar_tmp, bc_type="not-a-knot", extrapolate=True)
                dapardpsi[:, j, k] = cs(psixy[:, j], 1)

        yi_sol = np.arange(1, ny + 2, dtype=float)
        yi_sol_eval = np.arange(1, ny + 1, dtype=float)
        yi_cfr = np.arange(1, ny_cfr + 1, dtype=float)
        yi_cfr_eval = np.arange(1, ny_cfr, dtype=float)
        yi_pfr = np.arange(1, 2 * nypf1 + 2, dtype=float)
        yi_pfr_eval = np.arange(1, 2 * nypf1 + 1, dtype=float)

        for k in range(nz):
            for i in range(nx):
                if i < ixsep:
                    zarray_shift = np.mod(zarray[k] + sa[i], zmax)
                    apar_tmp = np.concatenate([apar[i, nypf1, :], apar[i, nypf1, 0:1]])
                    aparp1 = float(interp1_spline(zarrayp, apar_tmp, zarray_shift))

                    apar_cfr = np.concatenate([apar[i, iy_cfr, k], [aparp1]])
                    cs_cfr = CubicSpline(yi_cfr, apar_cfr, bc_type="not-a-knot", extrapolate=True)
                    dapardy[i, iy_cfr, k] = cs_cfr(yi_cfr_eval, 1)

                    apar_pfr = np.concatenate([apar[i, iy_pfr, k], [apar[i, iy_pfr[-1], k]]])
                    cs_pfr = CubicSpline(yi_pfr, apar_pfr, bc_type="not-a-knot", extrapolate=True)
                    dapardy[i, iy_pfr, k] = cs_pfr(yi_pfr_eval, 1)
                else:
                    apar_sol = np.concatenate([apar[i, :, k], [apar[i, -1, k]]])
                    cs_sol = CubicSpline(yi_sol, apar_sol, bc_type="not-a-knot", extrapolate=True)
                    dapardy[i, :, k] = cs_sol(yi_sol_eval, 1)

        dapardy /= dy0

        for j in range(ny):
            for i in range(nx):
                apar_tmp = np.concatenate([apar[i, j, :], apar[i, j, 0:1]])
                cs = CubicSpline(zarrayp, apar_tmp, bc_type="not-a-knot", extrapolate=True)
                dapardz[i, j, :] = cs(zarray, 1)

        for k in range(nz):
            dapardx[:, :, k] = -dapardpsi[:, :, k] + sinty * dapardz[:, :, k]
    else:
        raise ValueError("Unknown derivative method for apar")

    return apar, dapardx, dapardy, dapardz
