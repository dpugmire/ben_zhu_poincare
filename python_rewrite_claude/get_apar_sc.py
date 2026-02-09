"""
Prepare and calculate Apar and its derivatives for field-line tracing
for shifted-circular configuration.

Note that psi is defined at CELL_YLOW, while B field is at CELL_CENTER,
so interpolation/shifting along y is needed UNLESS true_apar is True.

Translated from MATLAB (B. Zhu) to Python.

Indexing convention:
    1-based index variables for interpolation coordinates, 0-based for
    Python array access. See fieldline_tracing.py header for details.
"""

import numpy as np
from scipy.interpolate import CubicSpline


def get_apar_sc(psi, bxy, psixy, zs, sa, sinty, dy0, dz,
                zperiod, interp_opt, deriv_opt, true_apar):
    """
    Compute Apar and its derivatives for shifted-circular configuration.

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
    zperiod : int
        Toroidal periodicity factor.
    interp_opt : int
        0 = linear interpolation, 1 = cubic spline interpolation.
    deriv_opt : int
        0 = central finite difference, 1 = cubic spline differentiation.
    true_apar : bool
        If True, input psi is already true Apar at cell center.

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

    # For circular case, the entire y range is the closed flux region
    ny_cfr = ny + 1
    iy_cfr = np.arange(ny)  # 0-based indices [0, 1, ..., ny-1]

    zarray = np.arange(nz) * dz          # [0, dz, ..., (nz-1)*dz]
    zarrayp = np.arange(nz + 1) * dz     # [0, dz, ..., nz*dz]
    zmax = dz * nz

    psi = psi.copy()  # avoid modifying input

    if not true_apar:
        # Step 1: interpolation of psi back to CELL_CENTER
        psis = psi.copy()

        # Shift first point across the branch-cut → becomes ny+1 cell
        psisp1 = np.zeros((nx, 1, nz))
        for i in range(nx):
            zarray_shift = np.mod(zarray + sa[i], zmax)
            psis_tmp = np.append(psis[i, 0, :], psis[i, 0, 0])  # wrap-around
            cs = CubicSpline(zarrayp, psis_tmp)
            psisp1[i, 0, :] = cs(zarray_shift)

        if interp_opt == 0:
            # Linear interpolation
            print("Linear interpolation of psi back to cell center.")
            for i in range(nx):
                # Closed flux surface region
                psi[i, :-1, :] = 0.5 * (psis[i, :-1, :] + psis[i, 1:, :])
                psi[i, -1, :] = 0.5 * (psisp1[i, 0, :] + psis[i, -1, :])

        elif interp_opt == 1:
            # Cubic spline interpolation
            print("Cubic spline interpolation of psi back to cell center.")
            for i in range(nx):
                for k in range(nz):
                    psis_cfr = np.append(psis[i, :, k], psisp1[i, 0, k])
                    ycoords = np.arange(1, ny_cfr + 1, dtype=float)
                    cs = CubicSpline(ycoords, psis_cfr)
                    psi[i, :, k] = cs(ycoords[:-1] + 0.5)
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
        # Use a representative dpsi (same indices as MATLAB: 101,39 → 100,38)
        dpsi = psixy[100, 38] - psixy[99, 38]
        dapardpsi[1:-1, :, :] = 0.5 * (apars[2:, :, :] - apars[:-2, :, :]) / dpsi

        dapardz[:, :, 0] = 0.5 * (apar[:, :, 1] - apar[:, :, -1]) / dz
        dapardz[:, :, 1:-1] = 0.5 * (apar[:, :, 2:] - apar[:, :, :-2]) / dz
        dapardz[:, :, -1] = 0.5 * (apar[:, :, 0] - apar[:, :, -2]) / dz

        for k in range(nz):
            dapardx[:, :, k] = dapardpsi[:, :, k] + sinty * dapardz[:, :, k]

        # d/dy (tricky at branch cut)
        for i in range(nx):
            dapardy[i, 1:-1, :] = 0.5 * (apar[i, 2:, :] - apar[i, :-2, :]) / dy0

            # Reverse shift last point → becomes -1 cell
            zarray_shift = np.mod(zarray - sa[i], zmax)
            apar_tmp = np.append(apar[i, -1, :], apar[i, -1, 0])
            cs = CubicSpline(zarrayp, apar_tmp)
            aparm1 = cs(zarray_shift)
            # NOTE: MATLAB original uses hardcoded index 5: dapardy(i,5,:)
            # This appears to be a config-specific bug (should be index 1
            # in 1-based, i.e., the first boundary point at the branch cut).
            # Python uses the correct dynamic index: 0 (0-based).
            dapardy[i, 0, :] = 0.5 * (apar[i, 1, :] - aparm1) / dy0

            # Forward shift first point → becomes ny+1 cell
            zarray_shift = np.mod(zarray + sa[i], zmax)
            apar_tmp = np.append(apar[i, 0, :], apar[i, 0, 0])
            cs = CubicSpline(zarrayp, apar_tmp)
            aparp1 = cs(zarray_shift)
            dapardy[i, -1, :] = 0.5 * (aparp1 - apar[i, -2, :]) / dy0

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
        yiarray_cfr = np.arange(1, ny_cfr + 1, dtype=float)
        for k in range(nz):
            for i in range(nx):
                # Forward shift first point → ny+1 cell
                zarray_shift = np.mod(zarray[k] + sa[i], zmax)
                apar_tmp_z = np.append(apar[i, 0, :], apar[i, 0, 0])
                cs_z = CubicSpline(zarrayp, apar_tmp_z)
                aparp1 = cs_z(zarray_shift)

                apar_tmp_y = np.append(apar[i, :, k], aparp1)
                cs_y = CubicSpline(yiarray_cfr, apar_tmp_y)
                dapardy[i, :, k] = cs_y.c[2, :]
        dapardy = dapardy / dy0

        for j in range(ny):
            for i in range(nx):
                apar_tmp = np.append(apar[i, j, :], apar[i, j, 0])
                cs = CubicSpline(zarrayp, apar_tmp)
                dapardz[i, j, :] = cs.c[2, :]

        for k in range(nz):
            dapardx[:, :, k] = dapardpsi[:, :, k] + sinty * dapardz[:, :, k]

    return apar, dapardx, dapardy, dapardz
