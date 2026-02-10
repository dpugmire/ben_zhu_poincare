"""
Field-line tracing for poloidal Poincare plot using precomputed apar data.

This script reads precomputed apar and derivatives from a NetCDF file,
loads grid information, and performs field-line tracing to generate
Poincare plots and trajectory data.

Usage:
    python trace_fieldlines.py --apar-file apar_data.nc --grid-file grid.nc
"""

import os
import sys
import argparse
import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline
import netCDF4

# Local module imports
from rk4_flt1 import rk4_flt1


def main():
    # --------------------------------------------------------------------------
    # Parse command-line arguments
    # --------------------------------------------------------------------------

    parser = argparse.ArgumentParser(
        description='Field-line tracing for Poincare plots')
    parser.add_argument('--apar-file', type=str, required=True,
                        help='Input NetCDF file with apar data')
    parser.add_argument('--grid-file', type=str, required=True,
                        help='BOUT++ grid file')
    parser.add_argument('--direction', type=int, default=1, choices=[1, -1],
                        help='Tracing direction: 1 (forward) or -1 (backward)')
    parser.add_argument('--nlines', type=int, default=256,
                        help='Number of field lines to trace')
    parser.add_argument('--range', dest='line_range', nargs=3, metavar=('X0', 'X1', 'N'),
                        default=None,
                        help='Trace N line values between X0 and X1 (supports floating point X0/X1)')
    parser.add_argument('--nturns', type=int, default=25,
                        help='Approximate number of poloidal turns')
    parser.add_argument('--lines', type=str, default=None,
                        help='Comma-separated list of specific line values to trace (e.g., "50,75,100.5")')
    args = parser.parse_args()

    # --------------------------------------------------------------------------
    # Load precomputed apar data
    # --------------------------------------------------------------------------

    print(f"Loading apar data from {args.apar_file} ...")

    ds_apar = netCDF4.Dataset(args.apar_file, 'r')

    # Read dimensions
    nx = int(ds_apar.nx)
    ny = int(ds_apar.ny)
    nz = int(ds_apar.nz)
    zperiod = int(ds_apar.zperiod)
    divertor = int(ds_apar.divertor)

    # Read apar arrays (shape: nx, ny, nz)
    apar0 = np.array(ds_apar.variables['apar'][:], dtype=float)
    dapardx0 = np.array(ds_apar.variables['dapardx'][:], dtype=float)
    dapardy0 = np.array(ds_apar.variables['dapardy'][:], dtype=float)
    dapardz0 = np.array(ds_apar.variables['dapardz'][:], dtype=float)

    if divertor == 1:
        ixsep = int(ds_apar.ixsep)
        nypf1 = int(ds_apar.nypf1)
        nypf2 = int(ds_apar.nypf2)

    # Store gridfile attribute for error messages (if available)
    try:
        apar_gridfile = ds_apar.gridfile
    except:
        apar_gridfile = None

    ds_apar.close()

    print(f"Loaded apar data: shape = {apar0.shape}")

    # --------------------------------------------------------------------------
    # Setup
    # --------------------------------------------------------------------------

    direction = args.direction
    nlines = args.nlines
    nturns = args.nturns
    nsteps = nturns * ny
    np_max = 1250  # maximum puncture points

    # Output options
    save_traj = True
    save_pp = True

    # Fill the torus if needed (zperiod > 1)
    nzG = nz * zperiod  # now nz extends full torus
    zmin = 0.0
    zmax = 2.0 * np.pi
    dz_torus = (zmax - zmin) / nzG

    # 1-based coordinate arrays for interpolation
    ziarray = np.arange(1, nzG + 2, dtype=float)       # [1, 2, ..., nzG+1]
    zarray  = (ziarray - 1) * dz_torus                   # [0, dz, ..., 2*pi]

    xiarray = np.arange(1, nx + 1, dtype=float)          # [1, 2, ..., nx]
    yiarray = np.arange(1, ny + 1, dtype=float)          # [1, 2, ..., ny]

    # --------------------------------------------------------------------------
    # Load grid information
    # --------------------------------------------------------------------------

    print("Loading grid information ...")

    ds = netCDF4.Dataset(args.grid_file, 'r')

    zShift = np.array(ds.variables['zShift'][:], dtype=float)
    rxy    = np.array(ds.variables['Rxy'][:], dtype=float)
    zxy    = np.array(ds.variables['Zxy'][:], dtype=float)
    psixy  = np.array(ds.variables['psixy'][:], dtype=float)
    rmag   = float(ds.variables['rmag'][:])

    ixsep1 = int(ds.variables['ixseps1'][:])
    ixsep2 = int(ds.variables['ixseps2'][:])

    if ixsep2 < nx:
        divertor = 2  # double null
        print("\tDouble null configuration")
        nypf11 = int(ds.variables['jyseps1_1'][:]) + 1
        nypf21 = int(ds.variables['jyseps2_1'][:]) + 1
        nypf12 = int(ds.variables['jyseps1_2'][:]) + 1
        nypf22 = int(ds.variables['jyseps2_2'][:]) + 1

    elif ixsep1 < nx:
        divertor = 1
        print("\tSingle null configuration")
        ixsep = ixsep1
        nypf1 = int(ds.variables['jyseps1_1'][:]) + 1
        nypf2 = int(ds.variables['jyseps2_2'][:]) + 1

    else:
        divertor = 0
        ixsep = nx
        nypf1 = 0
        nypf2 = ny
        print("\tCircular configuration")

    # jyomp: y-index of the outboard midplane (1-based)
    jyomp = int(np.argmax(rxy[-1, :]) + 1)

    xarray = psixy[:, jyomp - 1]  # psi values at the outboard midplane
    xMin = np.min(xarray)
    xMax = np.max(xarray)

    # Magnetic field and geometry variables
    bxy   = np.array(ds.variables['Bxy'][:], dtype=float)
    btxy  = np.array(ds.variables['Btxy'][:], dtype=float)
    bpxy  = np.array(ds.variables['Bpxy'][:], dtype=float)
    hthe  = np.array(ds.variables['hthe'][:], dtype=float)
    sinty = np.array(ds.variables['sinty'][:], dtype=float)

    bxcvx = np.array(ds.variables['bxcvx'][:], dtype=float)
    bxcvy = np.array(ds.variables['bxcvy'][:], dtype=float)
    bxcvz = np.array(ds.variables['bxcvz'][:], dtype=float)

    jpar0 = np.array(ds.variables['Jpar0'][:], dtype=float)
    dy_raw = np.array(ds.variables['dy'][:], dtype=float)
    dy0 = float(dy_raw.flat[0])

    sa = np.array(ds.variables['ShiftAngle'][:], dtype=float).flatten()

    ds.close()

    # --------------------------------------------------------------------------
    # Validate dimensions match between apar file and grid file
    # --------------------------------------------------------------------------

    grid_nx, grid_ny = rxy.shape

    if nx != grid_nx or ny != grid_ny:
        print("\n" + "="*70)
        print("ERROR: Dimension mismatch between apar file and grid file!")
        print("="*70)
        print(f"\nApar file dimensions (from {args.apar_file}):")
        print(f"  nx = {nx}, ny = {ny}, nz = {nz}")
        print(f"\nGrid file dimensions (from {args.grid_file}):")
        print(f"  nx = {grid_nx}, ny = {grid_ny}")
        print(f"\nThe apar data was generated from a different grid file.")
        print(f"Grid file used to generate apar data:")
        if apar_gridfile is not None:
            print(f"  {apar_gridfile}")
        else:
            print(f"  (not recorded in apar file)")
        print("\nSolutions:")
        print("  1. Use the correct grid file that matches the apar data")
        print("  2. Regenerate apar data using the current grid file:")
        print(f"     python generate_apar_data.py --output new_apar.nc")
        print("="*70 + "\n")
        sys.exit(1)

    dz = 2.0 * np.pi / zperiod / nz
    nu = btxy * hthe / bpxy / rxy

    # --------------------------------------------------------------------------
    # Calculate SOL, x-point locations and poloidal angle theta
    # --------------------------------------------------------------------------

    theta = np.zeros(ny)

    if divertor == 1:  # single null
        # Geometric center of the core
        core_rxy = rxy[0, nypf1:ny - nypf1]
        core_zxy = zxy[0, nypf1:ny - nypf1]
        center_x = 0.5 * (np.max(core_rxy) + np.min(core_rxy))
        center_y = 0.5 * (np.max(core_zxy) + np.min(core_zxy))

        # Separatrix indices (0-based Python) — rearranged y ordering
        tmp = np.concatenate([
            np.arange(0, nypf1),
            np.arange(ny - nypf1 - 1, nypf1 - 1, -1),
            np.arange(ny - nypf1, ny)
        ])
        sepx = 0.5 * (rxy[ixsep - 1, tmp] + rxy[ixsep, tmp])
        sepy = 0.5 * (zxy[ixsep - 1, tmp] + zxy[ixsep, tmp])

        # X-point location
        xpoint_x = 0.25 * (
            sepx[nypf1 - 1] + sepx[nypf1] +
            sepx[ny - nypf1 - 1] + sepx[ny - nypf1])
        xpoint_y = 0.25 * (
            sepy[nypf1 - 1] + sepy[nypf1] +
            sepy[ny - nypf1 - 1] + sepy[ny - nypf1])

        u = np.array([center_x - xpoint_x, center_y - xpoint_y, 0.0])

        for iy in range(ny):
            v = np.array([center_x - rxy[0, iy], center_y - zxy[0, iy], 0.0])
            theta[iy] = np.arctan2(np.linalg.norm(np.cross(u, v)), np.dot(u, v))

        theta = theta / np.pi
        itheta = np.argmax(theta)
        theta[itheta:ny] = 2.0 - theta[itheta:ny]
        itheta = np.argmax(theta)
        if itheta != ny - 1:
            theta[itheta:ny] = 4.0 - theta[itheta:ny]

        # Shift reference to (ix=1, iy=nypf1+1) in MATLAB = (ix=0, iy=nypf1) in Python
        theta = theta - theta[nypf1]

    elif divertor == 0:  # shifted-circular
        center_x = 0.5 * (np.max(rxy[0, :]) + np.min(rxy[0, :]))
        center_y = 0.5 * (np.max(zxy[0, :]) + np.min(zxy[0, :]))
        u = np.array([center_x - rxy[0, 0], center_y - zxy[0, 0], 0.0])

        for iy in range(ny):
            v = np.array([center_x - rxy[0, iy], center_y - zxy[0, iy], 0.0])
            theta[iy] = np.arctan2(np.linalg.norm(np.cross(u, v)), np.dot(u, v))

        theta = theta / np.pi
        itheta = np.argmax(theta)
        theta[itheta:ny] = 2.0 - theta[itheta:ny]
        itheta = np.argmax(theta)
        if itheta != ny - 1:
            theta[itheta:ny] = 4.0 - theta[itheta:ny]

    else:
        print("\t\tConfiguration to be implemented!")

    # --------------------------------------------------------------------------
    # Construct/patch closed flux surface region for better Poincare plot
    # --------------------------------------------------------------------------

    # 1-based coordinate array for CFR x-indices
    xiarray_cfr = np.arange(1, ixsep + 1, dtype=float)

    if divertor == 0:
        # 1-based y-coords for CFR: [1, 2, ..., ny+1] (one extra patched point)
        yiarray_cfr = np.arange(1, ny + 2, dtype=float)
        theta_cfr = np.zeros(ny + 1)
        theta_cfr[:ny] = theta[:ny]
        theta_cfr[-1] = 2.0  # theta is pi-based

    else:
        # 1-based y-coords: [nypf1+1, ..., nypf2+1]
        yiarray_cfr = np.arange(nypf1 + 1, nypf2 + 2, dtype=float)
        theta_cfr = theta[nypf1:nypf2 + 1].copy()
        theta_cfr[-1] = 2.0

    # rxy_cfr: closed flux region with one extra wrapped column
    rxy_cfr = rxy[:ixsep, nypf1:nypf2].copy()
    rxy_cfr = np.column_stack([rxy_cfr, rxy_cfr[:, 0]])

    zxy_cfr = zxy[:ixsep, nypf1:nypf2].copy()
    zxy_cfr = np.column_stack([zxy_cfr, zxy_cfr[:, 0]])

    zs_cfr = zShift[:ixsep, nypf1:nypf2].copy()
    # Additional zShift across the branch cut
    zs_last = (0.5 * (nu[:ixsep, nypf1] + nu[:ixsep, nypf2 - 1]) * dy0
                + zs_cfr[:, -1])
    zs_cfr = np.column_stack([zs_cfr, zs_last])

    # --------------------------------------------------------------------------
    # Geometric coefficients
    # --------------------------------------------------------------------------
    A1 = rxy * bpxy * btxy / hthe
    A2 = bxy ** 2
    A3 = sinty * A1
    JJ = 4.0 * np.pi * 1.0e-7 * bpxy / hthe / (bxy ** 2) * jpar0

    # --------------------------------------------------------------------------
    # Fill the full torus
    # --------------------------------------------------------------------------

    apar    = np.zeros((nx, ny, nzG))
    dapardx = np.zeros((nx, ny, nzG))
    dapardy = np.zeros((nx, ny, nzG))
    dapardz = np.zeros((nx, ny, nzG))

    if apar0.shape[2] == nzG:
        apar = apar0.copy()
        dapardx = dapardx0.copy()
        dapardy = dapardy0.copy()
        dapardz = dapardz0.copy()
    elif apar0.shape[2] == nz:
        for zp in range(zperiod):
            s = zp * nz
            e = (zp + 1) * nz
            apar[:, :, s:e]    = apar0
            dapardx[:, :, s:e] = dapardx0
            dapardy[:, :, s:e] = dapardy0
            dapardz[:, :, s:e] = dapardz0
    else:
        raise ValueError(
            f"Unexpected Apar z-dimension: {apar0.shape[2]} "
            f"(expected {nz} or {nzG}).")

    # --------------------------------------------------------------------------
    # Calculate perturbed field
    # --------------------------------------------------------------------------

    b0dgy = bpxy / hthe

    bdgx = np.zeros((nx, ny, nzG))
    bdgy = np.zeros((nx, ny, nzG))
    bdgz = np.zeros((nx, ny, nzG))
    dxdy = np.zeros((nx, ny, nzG))
    dzdy = np.zeros((nx, ny, nzG))

    for k in range(nzG):
        bdgx[:, :, k] = ((1.0 / bxy) *
                          (-A1 * dapardy[:, :, k] - A2 * dapardz[:, :, k]) +
                          apar[:, :, k] * bxcvx)
        bdgy[:, :, k] = ((1.0 / bxy) *
                          (A1 * dapardx[:, :, k] - A3 * dapardz[:, :, k]) +
                          apar[:, :, k] * (bxcvy + JJ))
        bdgz[:, :, k] = ((1.0 / bxy) *
                          (A2 * dapardx[:, :, k] + A3 * dapardz[:, :, k]) +
                          apar[:, :, k] * bxcvz)

        dxdy[:, :, k] = bdgx[:, :, k] / (b0dgy + bdgy[:, :, k])
        dzdy[:, :, k] = bdgz[:, :, k] / (b0dgy + bdgy[:, :, k])

    # --------------------------------------------------------------------------
    # Compute twist-shift boundary data for the closed flux surface region
    # --------------------------------------------------------------------------

    # p1: from y=nypf1+1 (first CFR point) forward-shifted across the branch cut
    dxdyt = dxdy[:, nypf1, :].copy()              # (nx, nzG)
    dxdyt = np.column_stack([dxdyt, dxdyt[:, 0]])  # (nx, nzG+1) — wrap-around
    dzdyt = dzdy[:, nypf1, :].copy()
    dzdyt = np.column_stack([dzdyt, dzdyt[:, 0]])

    dxdy_p1 = np.zeros((nx, nzG))
    dzdy_p1 = np.zeros((nx, nzG))
    for ix in range(ixsep):
        zarray_shift = np.mod(zarray[:nzG] + sa[ix], zmax)
        cs_dx = CubicSpline(zarray, dxdyt[ix, :])
        dxdy_p1[ix, :] = cs_dx(zarray_shift)
        cs_dz = CubicSpline(zarray, dzdyt[ix, :])
        dzdy_p1[ix, :] = cs_dz(zarray_shift)

    # m1: from y=nypf2 (last CFR point) reverse-shifted across the branch cut
    dxdyt = dxdy[:, nypf2 - 1, :].copy()
    dxdyt = np.column_stack([dxdyt, dxdyt[:, 0]])
    dzdyt = dzdy[:, nypf2 - 1, :].copy()
    dzdyt = np.column_stack([dzdyt, dzdyt[:, 0]])

    dxdy_m1 = np.zeros((nx, nzG))
    dzdy_m1 = np.zeros((nx, nzG))
    for ix in range(ixsep):
        zarray_rshift = np.mod(zarray[:nzG] - sa[ix], zmax)
        cs_dx = CubicSpline(zarray, dxdyt[ix, :])
        dxdy_m1[ix, :] = cs_dx(zarray_rshift)
        cs_dz = CubicSpline(zarray, dzdyt[ix, :])
        dzdy_m1[ix, :] = cs_dz(zarray_rshift)

    # --------------------------------------------------------------------------
    # Main field-line tracing loop
    # --------------------------------------------------------------------------

    print("Starting field-line tracing ...\n")

    # Create output directories
    os.makedirs('npz_traj', exist_ok=True)
    os.makedirs('npz_pp', exist_ok=True)

    ip_fid = open('ip_xyz.txt', 'w')
    ip_fid.write('iline it ipx ipy ipz\n')
    ip_thetapsi_fid = open('ip_thetapsi.txt', 'w')
    ip_thetapsi_fid.write('iline theta psi\n')
    traj_fid = open('traj_xyz.txt', 'w')
    traj_fid.write('iline it x y z\n')

    # Which lines to trace
    if args.line_range is not None:
        x0_str, x1_str, n_str = args.line_range
        try:
            x0 = float(x0_str)
            x1 = float(x1_str)
            n_samples = int(n_str)
        except ValueError:
            parser.error('--range expects X0 X1 N where X0/X1 are float and N is int')
        if n_samples <= 0:
            parser.error('--range requires N > 0')
        LINES = np.linspace(x0, x1, n_samples, dtype=float).tolist()
    elif args.lines is not None:
        try:
            LINES = [float(x.strip()) for x in args.lines.split(',') if x.strip()]
        except ValueError:
            parser.error('--lines must be a comma-separated list of numeric values')
        if len(LINES) == 0:
            parser.error('--lines list is empty')
    else:
        LINES = list(range(1, nlines + 1))

    for iline in LINES:
        if iline < 1 or iline > nx:
            print(f"\tSkipping line {iline:.16g}: outside valid index range [1, {nx}]")
            continue

        # Pick starting point (1-based indices)
        xind = float(iline)
        xStart = float(np.interp(xind, xiarray, xarray))  # psi value (supports float xind)
        yyy = jyomp  # 1-based
        yStart = jyomp  # 1-based
        zzz = 1  # 1-based (for filename)
        zStart = zarray[0]  # = 0.0

        # Trajectory and puncture point arrays
        traj = np.zeros((7, nsteps))
        fl_x3d = np.zeros(nsteps)
        fl_y3d = np.zeros(nsteps)
        fl_z3d = np.zeros(nsteps)
        px = np.zeros(np_max)
        py = np.zeros(np_max)
        pz = np.zeros(np_max)
        ptheta = np.zeros(np_max)
        ppsi = np.zeros(np_max)

        it = 0  # 0-based step counter (traj column index)
        iturn = 1

        # Determine initial region
        if xind < float(ixsep) + 0.5:
            region = 0  # closed flux surface
            if yStart < nypf1 + 1 or yStart > nypf2:
                region = 2  # private flux region
        else:
            region = 1  # SOL

        yind = yStart
        zind = float(np.interp(zStart, zarray, ziarray))

        print(f"\tline {iline} started at indices ({xind:.1f},{yind},{zind:.1f}),")

        # Check if starting on divertor
        if divertor == 1:
            if yStart == ny and direction == 1:
                print(f"\tline {iline} starts on the divertor.")
                region = 14
                traj[5, it] = 0.0
            elif yStart == 1 and direction == -1:
                print(f"\tline {iline} starts on the divertor.")
                region = 13
                traj[5, it] = 0.0

        # Record initial position
        traj[0, it] = 1       # turn number
        traj[1, it] = xind    # fractional x-index (1-based)
        traj[2, it] = yStart  # integer y-index (1-based)
        traj[3, it] = zind    # fractional z-index (1-based)
        traj[4, it] = region  # region flag
        traj[6, it] = zStart  # z angle value

        while region < 10 and iturn < nturns:

            if iturn % 50 == 1:
                print(f"\t\t line{iline}, turn {iturn}/{nturns} ...")

            for iy in range(ny - 1):

                # ----------------------------------------------------------
                # CFR stepping
                # ----------------------------------------------------------
                if region == 0 and yStart > nypf1 and yStart < nypf2 + 1:

                    if direction == 1:
                        xEnd, zEnd = rk4_flt1(
                            xStart, yStart, zStart, dxdy, dzdy,
                            xarray, zarray, region,
                            dxdy_p1, dzdy_p1, 1, nypf1, nypf2)
                        yEnd = yStart + 1
                    elif direction == -1:
                        xEnd, zEnd = rk4_flt1(
                            xStart, yStart, zStart, dxdy, dzdy,
                            xarray, zarray, region,
                            dxdy_m1, dzdy_m1, -1, nypf1, nypf2)
                        yEnd = yStart - 1

                    traj[6, it + 1] = zEnd

                    # Check field-line endpoint
                    if xEnd > xMax:
                        print(f"\tstarting xind={xind:.1f}, "
                              f"line {iline} reaches outer bndry")
                        region = 12
                    elif xEnd < xMin:
                        print(f"\tstarting xind={xind:.1f}, "
                              f"line {iline} reaches inner bndry")
                        region = 11
                    else:
                        xind = float(np.interp(xEnd, xarray, xiarray))
                        if xind > float(ixsep) + 0.5:
                            region = 1
                            print(f"\tending xind={xind:.1f}, "
                                  f"line {iline} enters the SOL")

                    # Twist-shift at the branch cut
                    if (direction == 1 and yStart == nypf2 and region == 0):
                        shiftangle = float(np.interp(xind, xiarray, sa))
                        zEnd = zEnd + shiftangle
                        yEnd = nypf1 + 1

                    if (direction == -1 and yStart == nypf1 + 1 and region == 0):
                        shiftangle = float(np.interp(xind, xiarray, sa))
                        zEnd = zEnd - shiftangle
                        yEnd = nypf2

                    # Re-label toroidal location if necessary
                    if zEnd < zmin or zEnd > zmax:
                        zEnd = zEnd % zmax
                    zind = float(np.interp(zEnd, zarray, ziarray))

                    it += 1
                    traj[0, it] = iturn
                    traj[1, it] = xind
                    traj[2, it] = yEnd
                    traj[3, it] = zind
                    traj[4, it] = region
                    traj[5, it] = hthe[int(np.round(xind)) - 1, yEnd - 1]

                    # End-point becomes new start-point
                    xStart = xEnd
                    yStart = yEnd
                    zStart = zEnd

                # ----------------------------------------------------------
                # SOL / PFR stepping
                # ----------------------------------------------------------
                elif region == 1 or region == 2:

                    if direction == 1:
                        xEnd, zEnd = rk4_flt1(
                            xStart, yStart, zStart, dxdy, dzdy,
                            xarray, zarray, region,
                            dxdy_p1, dzdy_p1, 1, nypf1, nypf2)
                        yEnd = yStart + 1
                    elif direction == -1:
                        xEnd, zEnd = rk4_flt1(
                            xStart, yStart, zStart, dxdy, dzdy,
                            xarray, zarray, region,
                            dxdy_m1, dzdy_m1, -1, nypf1, nypf2)
                        yEnd = yStart - 1

                    traj[6, it + 1] = zEnd

                    # Correct yEnd for PFR wraparound
                    if (direction == 1 and yStart == nypf1 and region == 2):
                        yEnd = nypf2 + 1
                    elif (direction == -1 and yStart == nypf2 + 1 and region == 2):
                        yEnd = nypf1

                    # Check field-line endpoint
                    if xEnd > xMax:
                        print(f"\tstarting xind={xind:.1f}, "
                              f"line {iline} reaches outer bndry")
                        region = 12
                    elif xEnd < xMin:
                        print(f"\tstarting xind={xind:.1f}, "
                              f"line {iline} reaches inner bndry")
                        region = 11
                    else:
                        xind = float(np.interp(xEnd, xarray, xiarray))
                        if (xind < float(ixsep) + 0.5 and
                                yEnd > nypf1 and yEnd < nypf2 + 1):
                            if region != 0:
                                print(f"\tending xind={xind:.1f}, "
                                      f"line {iline} enters the CFR")
                            region = 0
                        elif (xind < float(ixsep) + 0.5 and
                              (yEnd > nypf2 - 1 or yEnd < nypf1)):
                            if region != 2:
                                print(f"\tending xind={xind:.1f}, "
                                      f"line {iline} enters the PFR")
                            region = 2

                    # Check divertor endpoints
                    if direction == 1 and yEnd == ny:
                        print(f"\tstarting xind={xind:.1f}, "
                              f"line {iline} reaches divertor")
                        region = 14
                    elif direction == -1 and yEnd == 1:
                        print(f"\tstarting xind={xind:.1f}, "
                              f"line {iline} reaches divertor")
                        region = 13

                    # Re-label toroidal location if necessary
                    if zEnd < zmin or zEnd > zmax:
                        zEnd = zEnd % zmax
                    zind = float(np.interp(zEnd, zarray, ziarray))

                    it += 1
                    traj[0, it] = iturn
                    traj[1, it] = xind
                    traj[2, it] = yEnd
                    traj[3, it] = zind
                    traj[4, it] = region
                    traj[5, it] = hthe[int(np.round(xind)) - 1, yStart - 1]

                    # End-point becomes new start-point
                    xStart = xEnd
                    yStart = yEnd
                    zStart = zEnd

                # End of single step (CFR/SOL/PFR)

                # If region changed to >= 10, break out of the y loop
                if region >= 10:
                    break

            # End of approximate one turn (ny steps)
            iturn += 1

        # End of nturns

        # Record maximum valid steps
        itmax = it + 1  # number of valid entries (0-based it → itmax entries)

        traj = traj[:, :itmax]

        if save_traj:
            suffix = 'p' if direction == 1 else 'm'
            fname = (f"npz_traj/x{iline}y{yyy}z{zzz}"
                     f"_v3lc-01-250{suffix}.npz")
            np.savez(fname, v=traj)

        # ------------------------------------------------------------------
        # Calculate puncture points and generate Poincare plot data
        # ------------------------------------------------------------------

        # Compute 3D Cartesian coordinates of the field-line trajectory
        for istep in range(itmax):
            yidx = int(traj[2, istep]) - 1  # 0-based y-index for array access
            rxyvalue = float(np.interp(traj[1, istep], xiarray, rxy[:, yidx]))
            zsvalue  = float(np.interp(traj[1, istep], xiarray, zShift[:, yidx]))
            zvalue   = float(np.interp(traj[3, istep], ziarray, zarray))
            x3d_tmp = rxyvalue * np.cos(zsvalue)
            y3d_tmp = rxyvalue * np.sin(zsvalue)
            fl_x3d[istep] = x3d_tmp * np.cos(zvalue) - y3d_tmp * np.sin(zvalue)
            fl_y3d[istep] = x3d_tmp * np.sin(zvalue) + y3d_tmp * np.cos(zvalue)
            fl_z3d[istep] = float(np.interp(traj[1, istep], xiarray, zxy[:, yidx]))

        for istep in range(itmax):
            traj_fid.write(
                f"{iline} {istep + 1} "
                f"{fl_x3d[istep]:.16g} {fl_y3d[istep]:.16g} "
                f"{fl_z3d[istep]:.16g}\n")

        # ------------------------------------------------------------------
        # Get puncture points for Poincare plot
        # ------------------------------------------------------------------

        if itmax > 1:
            # Find zero-crossings in fl_x3d (x = 0 plane in Cartesian)
            # Use fine spline interpolation
            itarray = np.arange(1, itmax + 1, dtype=float)  # 1-based
            fit = np.arange(1, itmax + 0.0001, 0.0001)
            cs_x3d = CubicSpline(itarray, fl_x3d[:itmax])
            ffl_x3d = cs_x3d(fit)

            # Zero-crossing detection
            v_prod = ffl_x3d[:-1] * ffl_x3d[1:]
            iit = np.where(v_prod <= 0)[0]

            nc = len(iit)
            ip = 0
            id_count = 0

            if nc > 0:
                for i in range(nc):
                    it_val = fit[iit[i]]
                    it_idx = int(np.floor(it_val))  # 1-based traj column index
                    a = it_val - it_idx
                    b = 1.0 - a

                    # Traj column indices: it_idx-1 and it_idx (0-based Python)
                    tc0 = it_idx - 1  # 0-based
                    tc1 = it_idx      # 0-based

                    # Ensure we don't exceed array bounds
                    if tc1 >= itmax:
                        continue

                    # Linear interpolation along field-line
                    xind_tmp = b * traj[1, tc0] + a * traj[1, tc1]
                    yind_tmp = b * traj[2, tc0] + a * traj[2, tc1]

                    # Raw zEnd interpolation
                    zvalue = b * traj[6, tc0] + a * traj[6, tc1]
                    if abs(traj[6, tc0] - traj[6, tc1]) > 1.0:
                        zvalue = (b * (traj[6, tc0] % zmax) +
                                  a * (traj[6, tc1] % zmax))

                    # Edge case adjustments at branch cut
                    if (traj[2, tc0] == float(nypf2) and direction == 1 and
                            xind_tmp < float(ixsep) + 0.5):
                        yind_tmp = b * traj[2, tc0] + a * float(nypf2 + 1)
                    elif (traj[2, tc0] == float(nypf1 + 1) and direction == -1 and
                          xind_tmp < float(ixsep) + 0.5):
                        yind_tmp = b * float(nypf2 + 1) + a * traj[2, tc1]
                        shiftangle = float(np.interp(xind_tmp, xiarray, sa))
                        zvalue = (zvalue - shiftangle) % zmax
                    elif (tc0 > 0 and
                          (traj[2, tc0 - 1] == float(nypf2) or
                           traj[2, tc0 - 1] == float(nypf1 + 1))):
                        zvalue = (
                            b * float(np.interp(traj[3, tc0], ziarray, zarray)) +
                            a * float(np.interp(traj[3, tc1], ziarray, zarray)))

                    # 2D spline interpolation for R, Z, zShift at puncture point
                    if xind_tmp < float(ixsep) + 0.5:
                        spl_r = RectBivariateSpline(
                            xiarray_cfr, yiarray_cfr, rxy_cfr)
                        spl_z = RectBivariateSpline(
                            xiarray_cfr, yiarray_cfr, zxy_cfr)
                        spl_zs = RectBivariateSpline(
                            xiarray_cfr, yiarray_cfr, zs_cfr)
                        rxyvalue = float(spl_r.ev(xind_tmp, yind_tmp))
                        zxyvalue = float(spl_z.ev(xind_tmp, yind_tmp))
                        zsvalue  = float(spl_zs.ev(xind_tmp, yind_tmp))
                    else:
                        spl_r = RectBivariateSpline(xiarray, yiarray, rxy)
                        spl_z = RectBivariateSpline(xiarray, yiarray, zxy)
                        spl_zs = RectBivariateSpline(xiarray, yiarray, zShift)
                        rxyvalue = float(spl_r.ev(xind_tmp, yind_tmp))
                        zxyvalue = float(spl_z.ev(xind_tmp, yind_tmp))
                        zsvalue  = float(spl_zs.ev(xind_tmp, yind_tmp))

                    # Convert to 3D Cartesian
                    ipx3d_tmp = rxyvalue * np.cos(zsvalue)
                    ipy3d_tmp = rxyvalue * np.sin(zsvalue)
                    ipx = ipx3d_tmp * np.cos(zvalue) - ipy3d_tmp * np.sin(zvalue)
                    ipy = ipx3d_tmp * np.sin(zvalue) + ipy3d_tmp * np.cos(zvalue)
                    ipz = zxyvalue

                    if ipy > 0:
                        px[ip] = ipx
                        py[ip] = ipy
                        pz[ip] = ipz
                        ip_fid.write(
                            f"{iline} {it_idx} "
                            f"{ipx:.16g} {ipy:.16g} {ipz:.16g}\n")

                        ptheta[ip] = float(
                            np.interp(yind_tmp, yiarray_cfr, theta_cfr))
                        ppsi[ip] = float(
                            np.interp(xind_tmp, xiarray, xarray))

                        # Write theta (converted to radians 0-2π) and psi
                        theta_rad = ptheta[ip] * np.pi  # Convert from π units to radians
                        ip_thetapsi_fid.write(
                            f"{iline} {theta_rad:.16g} {ppsi[ip]:.16g}\n")

                        ip += 1

                print(f"\t\tline {iline} has {ip}(+{id_count}) "
                      f"interception points.")

            # Pack puncture point data
            if ip > 0:
                pxp = px[:ip]
                pyp = py[:ip]
                pzp = pz[:ip]
                ptp = ptheta[:ip]
                ppp = ppsi[:ip]
            else:
                pxp = np.array([0.0])
                pyp = np.array([0.0])
                pzp = np.array([0.0])
                ptp = np.array([0.0])
                ppp = np.array([0.0])

            if nc == 0:
                pxp = np.array([0.0])
                pyp = np.array([0.0])
                pzp = np.array([0.0])
                ptp = np.array([0.0])
                ppp = np.array([0.0])

            # 3D trajectory for output
            traj0 = np.zeros((3, itmax))
            traj0[0, :itmax] = fl_x3d[:itmax]
            traj0[1, :itmax] = fl_y3d[:itmax]
            traj0[2, :itmax] = fl_z3d[:itmax]

            # Connection length
            lc = np.sum(traj[5, :itmax])
            print(f"\t\tline {iline}: connection length is {lc:.6f} "
                  f"and ends at region={region}.")

            # Save puncture point info
            if save_pp:
                suffix = 'p' if direction == 1 else 'm'
                fname = (f"npz_pp/x{iline}y{yyy}z{zzz}"
                         f"_v3lc-00-250{suffix}.npz")
                np.savez(fname,
                         v1=pxp, v2=pyp, v3=pzp,
                         v4=ptp, v5=ppp, v6=traj0,
                         v7=np.array([lc]), v8=np.array([region]))

        # End itmax > 1

    # End iline loop

    ip_fid.close()
    ip_thetapsi_fid.close()
    traj_fid.close()

    print("\nField-line tracing complete.")


if __name__ == '__main__':
    main()
