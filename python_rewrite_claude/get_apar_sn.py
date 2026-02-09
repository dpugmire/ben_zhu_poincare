"""
Prepare and calculate Apar and its derivatives for field-line tracing
for single-null configuration.

Note that psi is defined at CELL_YLOW, while B field is at CELL_CENTER,
so interpolation/shifting along y is needed UNLESS true_apar is True.

Translated from MATLAB (B. Zhu) to Python.

Indexing convention:
    1-based index variables for interpolation coordinates.
    ixsep, nypf1, nypf2 use 1-based MATLAB convention (matching the values
    in the original code). Python array access uses 0-based (with -1 offset).
"""

import numpy as np
from scipy.interpolate import CubicSpline


def get_apar_sn(psi, bxy, psixy, zs, sa, sinty, dy0, dz,
                ixsep, nypf1, nypf2, zperiod, interp_opt, deriv_opt,
                true_apar):
    """
    Compute Apar and its derivatives for single-null configuration.

    Parameters
    ----------
    psi : ndarray, shape (nx, ny, nz)
        BOUT++ psi output (or true Apar if true_apar=True).
    bxy : ndarray, shape (nx, ny)
        Total magnetic field strength.
    psixy : ndarray, shape (nx, ny)
        Equilibrium poloidal flux.
    zs : ndarray, shape (nx, ny)
        zShift field.
    sa : ndarray, shape (nx,)
        ShiftAngle.
    sinty : ndarray, shape (nx, ny)
        Integrated shear (sinty field).
    dy0 : float
        Poloidal grid spacing.
    dz : float
        Toroidal grid spacing.
    ixsep : int
        Separatrix x-index (1-based MATLAB convention).
    nypf1 : int
        Lower PF boundary index (1-based: jyseps1_1 + 1 from file).
    nypf2 : int
        Upper PF boundary index (1-based: jyseps2_2 + 1 from file).
    zperiod : int
        Toroidal periodicity factor.
    interp_opt : int
        0 = linear interpolation, 1 = cubic spline interpolation.
    deriv_opt : int
        0 = central finite difference, 1 = cubic spline differentiation.
    true_apar : bool
        If True, input is already true Apar at cell center.

    Returns
    -------
    apar : ndarray, shape (nx, ny, nz)
    dapardx : ndarray, shape (nx, ny, nz)
    dapardy : ndarray, shape (nx, ny, nz)
    dapardz : ndarray, shape (nx, ny, nz)
    """
    nx, ny, nz = psi.shape

    apar = np.zeros((nx, ny, nz))
    dapardx = np.zeros((nx, ny, nz))
    dapardy = np.zeros((nx, ny, nz))
    dapardz = np.zeros((nx, ny, nz))

    ny_cfr = nypf2 - nypf1 + 1  # number of closed flux region y-points + 1

    # 0-based index arrays for the closed flux region and private flux region
    # MATLAB iy_cfr = [nypf1+1:nypf2] (1-based) → Python [nypf1:nypf2] (0-based)
    iy_cfr = np.arange(nypf1, nypf2)  # 0-based indices into arrays
    # MATLAB iy_pfr = [1:nypf1, nypf2+1:ny] → Python [0:nypf1, nypf2:ny]
    iy_pfr = np.concatenate([np.arange(0, nypf1), np.arange(nypf2, ny)])

    zarray = np.arange(nz) * dz          # [0, dz, ..., (nz-1)*dz]
    zarrayp = np.arange(nz + 1) * dz     # [0, dz, ..., nz*dz]
    zmax = dz * nz

    psi = psi.copy()

    if not true_apar:
        # Step 1: interpolation of psi back to CELL_CENTER
        psis = psi.copy()

        # Shift first CFR point across branch-cut → becomes ny+1 cell
        psisp1 = np.zeros((ixsep, 1, nz))
        for i in range(ixsep):  # MATLAB 1:ixsep → Python 0:ixsep
            zarray_shift = np.mod(zarray + sa[i], zmax)
            # MATLAB: psi(i, nypf1+1, :) → Python: psi[i, nypf1, :]
            psis_tmp = np.append(psis[i, nypf1, :], psis[i, nypf1, 0])
            cs = CubicSpline(zarrayp, psis_tmp)
            psisp1[i, 0, :] = cs(zarray_shift)

        if interp_opt == 0:
            # Linear interpolation
            print("Linear interpolation of psi back to cell center.")
            for i in range(nx):
                if i < ixsep:
                    # Closed flux surface region
                    psi[i, iy_cfr[:-1], :] = 0.5 * (
                        psis[i, iy_cfr[:-1], :] + psis[i, iy_cfr[1:], :])
                    # MATLAB: psi(i, nypf2, :) → Python: psi[i, nypf2-1, :]
                    psi[i, nypf2 - 1, :] = 0.5 * (
                        psisp1[i, 0, :] + psis[i, nypf2 - 1, :])
                    # Private flux region
                    psi[i, iy_pfr[:-1], :] = 0.5 * (
                        psis[i, iy_pfr[:-1], :] + psis[i, iy_pfr[1:], :])
                    # Simple linear extrapolation for last PFR point
                    psi[i, iy_pfr[-1], :] = (
                        1.5 * psis[i, iy_pfr[-1], :] -
                        0.5 * psis[i, iy_pfr[-2], :])
                else:
                    # SOL
                    psi[i, :-1, :] = 0.5 * (psis[i, :-1, :] + psis[i, 1:, :])
                    # Linear extrapolation for last point
                    psi[i, -1, :] = 1.5 * psis[i, -1, :] - 0.5 * psis[i, -2, :]

        elif interp_opt == 1:
            # Cubic spline interpolation
            print("Cubic spline interpolation of psi back to cell center.")
            for i in range(nx):
                if i < ixsep:
                    for k in range(nz):
                        # CFR
                        psis_cfr = np.append(psis[i, iy_cfr, k], psisp1[i, 0, k])
                        ycoords = np.arange(1, ny_cfr + 1, dtype=float)
                        cs = CubicSpline(ycoords, psis_cfr)
                        psi[i, nypf1:nypf2, k] = cs(ycoords[:-1] + 0.5)
                        # PFR
                        psis_pfr = psis[i, iy_pfr, k].copy()
                        psis_pfr = np.append(
                            psis_pfr,
                            2.0 * psis_pfr[-1] - psis_pfr[-2])
                        n_pfr = len(iy_pfr)
                        ycoords_pfr = np.arange(1, n_pfr + 2, dtype=float)
                        cs_pfr = CubicSpline(ycoords_pfr, psis_pfr)
                        psi[i, iy_pfr, k] = cs_pfr(
                            np.arange(1, n_pfr + 1, dtype=float) + 0.5)
                else:
                    for k in range(nz):
                        psis_tmp = np.append(psis[i, :, k], psis[i, -1, k])
                        ycoords = np.arange(1, ny + 2, dtype=float)
                        cs = CubicSpline(ycoords, psis_tmp)
                        psi[i, :, k] = cs(np.arange(1, ny + 1, dtype=float) + 0.5)
        else:
            print("Unknown interpolation method for psi!")

        # Step 2: apar = psi * bxy
        for k in range(nz):
            apar[:, :, k] = psi[:, :, k] * bxy

    else:
        # Input is true Apar at cell center
        apar = psi.copy()

    # Step 3: get Apar derivatives

    # Shift to flux coordinate for d/dpsi
    apars = np.zeros((nx, ny, nz))
    dapardpsi = np.zeros((nx, ny, nz))
    kz = np.concatenate([np.arange(0, nz // 2 + 1),
                         np.arange(-nz // 2 + 1, 0)]) * zperiod
    for i in range(nx):
        for j in range(ny):
            apars[i, j, :] = np.real(
                np.fft.ifft(np.fft.fft(apar[i, j, :]) *
                            np.exp(-1j * zs[i, j] * kz))
            )

    if deriv_opt == 0:
        # Central finite difference
        print("Central finite difference for Apar.")
        dpsi = abs(psixy[100, 38] - psixy[99, 38])
        dapardpsi[1:-1, :, :] = (
            0.5 * (apars[2:, :, :] - apars[:-2, :, :]) / dpsi)

        dapardz[:, :, 0] = 0.5 * (apar[:, :, 1] - apar[:, :, -1]) / dz
        dapardz[:, :, 1:-1] = 0.5 * (apar[:, :, 2:] - apar[:, :, :-2]) / dz
        dapardz[:, :, -1] = 0.5 * (apar[:, :, 0] - apar[:, :, -2]) / dz

        for k in range(nz):
            dapardx[:, :, k] = dapardpsi[:, :, k] + sinty * dapardz[:, :, k]

        # d/dy (tricky at branch cuts)
        for i in range(nx):
            if i < ixsep:
                # CFR interior
                dapardy[i, iy_cfr[1:-1], :] = 0.5 * (
                    apar[i, iy_cfr[2:], :] - apar[i, iy_cfr[:-2], :]) / dy0

                # Reverse shift last CFR point → becomes -1 cell
                zarray_shift = np.mod(zarray - sa[i], zmax)
                # MATLAB: apar(i, nypf2, :) → Python: apar[i, nypf2-1, :]
                apar_tmp = np.append(apar[i, nypf2 - 1, :],
                                     apar[i, nypf2 - 1, 0])
                cs = CubicSpline(zarrayp, apar_tmp)
                aparm1 = cs(zarray_shift)
                # NOTE: MATLAB original uses hardcoded index 5: dapardy(i,5,:)
                # This appears to be a config-specific bug (should be nypf1+1
                # in 1-based, i.e., the first CFR boundary point).
                # Python uses the correct dynamic index: nypf1 (0-based).
                dapardy[i, nypf1, :] = 0.5 * (
                    apar[i, nypf1 + 1, :] - aparm1) / dy0

                # Forward shift first CFR point → becomes ny+1 cell
                zarray_shift = np.mod(zarray + sa[i], zmax)
                apar_tmp = np.append(apar[i, nypf1, :], apar[i, nypf1, 0])
                cs = CubicSpline(zarrayp, apar_tmp)
                aparp1 = cs(zarray_shift)
                dapardy[i, nypf2 - 1, :] = 0.5 * (
                    aparp1 - apar[i, nypf2 - 2, :]) / dy0

                # PFR (parallel gradient vanishes at divertor target)
                dapardy[i, iy_pfr[1:-1], :] = 0.5 * (
                    apar[i, iy_pfr[2:], :] - apar[i, iy_pfr[:-2], :]) / dy0
            else:
                # SOL
                dapardy[i, 1:-1, :] = 0.5 * (
                    apar[i, 2:, :] - apar[i, :-2, :]) / dy0

    elif deriv_opt == 1:
        # Cubic spline fit then differentiation
        print("Cubic spline fit then differentiation of Apar.")

        for k in range(nz):
            for j in range(ny):
                apar_tmp = np.append(apars[:, j, k], apars[-1, j, k])
                xcoords = np.append(psixy[:, j],
                                    2.0 * psixy[-1, j] - psixy[-2, j])
                cs = CubicSpline(xcoords, apar_tmp)
                dapardpsi[:, j, k] = cs.c[2, :]

        # d/dy
        yiarray_sol = np.arange(1, ny + 2, dtype=float)
        yiarray_cfr = np.arange(1, ny_cfr + 1, dtype=float)

        for k in range(nz):
            for i in range(nx):
                if i < ixsep:
                    # CFR
                    zarray_shift = np.mod(zarray[k] + sa[i], zmax)
                    apar_tmp_z = np.append(apar[i, nypf1, :],
                                           apar[i, nypf1, 0])
                    cs_z = CubicSpline(zarrayp, apar_tmp_z)
                    aparp1 = cs_z(zarray_shift)

                    apar_tmp_y = np.append(apar[i, iy_cfr, k], aparp1)
                    cs_y = CubicSpline(yiarray_cfr, apar_tmp_y)
                    dapardy[i, iy_cfr, k] = cs_y.c[2, :]

                    # PFR
                    apar_tmp_pfr = np.append(apar[i, iy_pfr, k],
                                             apar[i, iy_pfr[-1], k])
                    n_pfr = len(iy_pfr)
                    ycoords_pfr = np.arange(1, n_pfr + 2, dtype=float)
                    cs_pfr = CubicSpline(ycoords_pfr, apar_tmp_pfr)
                    dapardy[i, iy_pfr, k] = cs_pfr.c[2, :]

                else:
                    # SOL
                    apar_tmp = np.append(apar[i, :, k], apar[i, -1, k])
                    cs = CubicSpline(yiarray_sol, apar_tmp)
                    dapardy[i, :, k] = cs.c[2, :]

        dapardy = dapardy / dy0

        for j in range(ny):
            for i in range(nx):
                apar_tmp = np.append(apar[i, j, :], apar[i, j, 0])
                cs = CubicSpline(zarrayp, apar_tmp)
                dapardz[i, j, :] = cs.c[2, :]

        for k in range(nz):
            dapardx[:, :, k] = -dapardpsi[:, :, k] + sinty * dapardz[:, :, k]

    return apar, dapardx, dapardy, dapardz
