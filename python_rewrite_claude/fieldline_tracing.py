"""
Field-line tracing for poloidal Poincare plot and other analysis.

Translated from MATLAB (B. Zhu, last updated 12/2023) to Python.

Indexing convention:
    All "index variables" (ixsep, nypf1, nypf2, xiarray, yiarray, ziarray,
    xind, yind, zind, yStart, yEnd, jyomp, etc.) use **1-based** values,
    matching the MATLAB original.  They serve as coordinate values for
    interpolation and are stored in the traj array as-is.

    When accessing Python (numpy) arrays, we subtract 1 to convert from
    1-based to 0-based.  For example:
        MATLAB:  rxy(ixsep, jy)       ->  Python:  rxy[ixsep-1, jy-1]
        MATLAB:  rxy(1:ixsep, nypf1+1:nypf2)
                                       ->  Python:  rxy[:ixsep, nypf1:nypf2]
    The second pattern works because Python slicing is exclusive at the
    upper bound, so  rxy[:ixsep]  gives the first ixsep rows (indices 0
    through ixsep-1), same as MATLAB's  rxy(1:ixsep).
"""

import os
import sys
import argparse
import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline
import netCDF4

# Local module imports (same directory)
from rk4_flt1 import rk4_flt1
from get_apar_sc import get_apar_sc
from get_apar_sn import get_apar_sn


# --------------------------------------------------------------------------
# STEP 0: Parse command-line arguments
# --------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Field-line tracing for Poincare plots')
parser.add_argument('--export-apar', type=str, default=None,
                    help='Export apar0, dapardx0, dapardy0 to NetCDF file and exit')
args = parser.parse_args()

# --------------------------------------------------------------------------
# STEP 1: User setup
# --------------------------------------------------------------------------

divertorCase = 1

if divertorCase == 0:
    gridfile = '/Users/dpn/proj/bout++/nersc_data/circle-zonal/cbm18_dens3_0.5BS_516nx64ny.grid.nc'
    aparfile = '/Users/dpn/proj/bout++/nersc_data/circle-zonal/apar_cbm18_dens3_0.5BS_516nx64ny64nz_t500.npy'
    apar_raw = np.load(aparfile)
    print(f"Apar shape={apar_raw.shape}")
    divertor = 0
    nx = 516; ny = 64; nz = 64
    zperiod = 5

elif divertorCase == 1:
    gridfile = '/Users/dpn/proj/bout++/poincare/boutpp_poincare/data/kstar_30306_7850_psi085105_nx260ny128_f2_v0.nc'
    aparfile = '/Users/dpn/proj/bout++/nersc_data/xpoint_singlenull/apar_kstar_30306_7850_psi085105_nx260ny128_f2_nz256.mat'
    print(f"aparfile = {aparfile}")

    # Load .mat file (try scipy first, fall back to h5py for v7.3 files)
    try:
        import scipy.io as sio
        apar_raw = sio.loadmat(aparfile)['apar']
    except NotImplementedError:
        import h5py
        with h5py.File(aparfile, 'r') as f:
            apar_raw = np.array(f['apar']).T  # HDF5 stores in transposed order

    divertor = 1
    nx = 260; ny = 128; nz = 256
    zperiod = 1

else:
    print("To be implemented!")
    sys.exit(1)

# Field-line tracing direction: 1 (y index increasing); -1 (y index decreasing)
direction = 1
# Individual field-lines to be traced radially
nlines = 256
deltaix = 1; ixoffset = 1

# (Roughly) total poloidal turns
nturns = 25
nsteps = nturns * ny
np_max = 1250  # maximum puncture points

# Output options
save_traj = True
save_pp = True

# --------------------------------------------------------------------------
# No user setup in this section
# --------------------------------------------------------------------------

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
# STEP 1: Load grid info and simulation output
# --------------------------------------------------------------------------

print("Loading grid information ...")

ds = netCDF4.Dataset(gridfile, 'r')

# Note: Python netCDF4 reads in C order matching file dimension order.
# MATLAB reads in Fortran order then transposes.  The result is the same
# (nx, ny) shape in both languages, with no transpose needed in Python.

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
    # MATLAB: tmp = [1:nypf1  ny-nypf1:-1:nypf1+1  ny-nypf1+1:ny]
    tmp = np.concatenate([
        np.arange(0, nypf1),
        np.arange(ny - nypf1 - 1, nypf1 - 1, -1),
        np.arange(ny - nypf1, ny)
    ])
    sepx = 0.5 * (rxy[ixsep - 1, tmp] + rxy[ixsep, tmp])
    sepy = 0.5 * (zxy[ixsep - 1, tmp] + zxy[ixsep, tmp])

    # X-point location (average of four adjacent separatrix points)
    # In 1-based: sepx at indices nypf1, nypf1+1, ny-nypf1, ny-nypf1+1
    # tmp was built with rearranged ordering, so we need the indices into tmp:
    #   nypf1 elements in part1, then the reversed core, then nypf1 elements in part3
    # The x-point is near the boundary between part1/part2 and part2/part3
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

# rxy_cfr: closed flux region (ixsep x ny_cfr_pts) with one extra wrapped column
# MATLAB: rxy(1:ixsep, nypf1+1:nypf2)  →  Python: rxy[:ixsep, nypf1:nypf2]
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
# Load perturbed field information
# --------------------------------------------------------------------------

print("Loading perturbed field information ...")

apar    = np.zeros((nx, ny, nzG))
dapardx = np.zeros((nx, ny, nzG))
dapardy = np.zeros((nx, ny, nzG))
dapardz = np.zeros((nx, ny, nzG))

if divertor == 0:
    apar0, dapardx0, dapardy0, dapardz0 = get_apar_sc(
        apar_raw, bxy, psixy, zShift, sa, sinty, dy0, dz,
        zperiod, 0, 1, False)

elif divertor == 1:
    apar0, dapardx0, dapardy0, dapardz0 = get_apar_sn(
        apar_raw, bxy, psixy, zShift, sa, sinty, dy0, dz,
        ixsep, nypf1, nypf2, zperiod, 0, 1, False)

else:
    print("\tConfiguration to be implemented!")
    sys.exit(1)

# --------------------------------------------------------------------------
# Optional: Export apar data and exit
# --------------------------------------------------------------------------

if args.export_apar is not None:
    print(f"\nExporting apar data to {args.export_apar} ...")

    # Create NetCDF file
    ncout = netCDF4.Dataset(args.export_apar, 'w', format='NETCDF4')

    # Create dimensions
    ncout.createDimension('nx', nx)
    ncout.createDimension('ny', ny)
    ncout.createDimension('nz', nz)

    # Create variables
    var_apar = ncout.createVariable('apar', 'f8', ('nx', 'ny', 'nz'))
    var_dapardx = ncout.createVariable('dapardx', 'f8', ('nx', 'ny', 'nz'))
    var_dapardy = ncout.createVariable('dapardy', 'f8', ('nx', 'ny', 'nz'))

    # Add attributes
    var_apar.long_name = 'Parallel vector potential'
    var_dapardx.long_name = 'x-derivative of parallel vector potential'
    var_dapardy.long_name = 'y-derivative of parallel vector potential'

    # Write data
    var_apar[:, :, :] = apar0
    var_dapardx[:, :, :] = dapardx0
    var_dapardy[:, :, :] = dapardy0

    # Add global attributes
    ncout.description = 'Precomputed parallel vector potential and derivatives'
    ncout.divertor_case = divertorCase
    ncout.divertor = divertor
    if divertor == 1:
        ncout.ixsep = ixsep
        ncout.nypf1 = nypf1
        ncout.nypf2 = nypf2
    ncout.zperiod = zperiod

    ncout.close()

    print(f"Successfully exported apar data to {args.export_apar}")
    print("Shape: apar0 =", apar0.shape)
    print("Exiting as requested.")
    sys.exit(0)

# --------------------------------------------------------------------------
# Fill the full torus
# --------------------------------------------------------------------------
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
# STEP 2: Field-line tracing
# --------------------------------------------------------------------------

# Calculate perturbed field
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
# MATLAB: dxdy(:, nypf1+1, :)  →  Python: dxdy[:, nypf1, :]
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
# MATLAB: dxdy(:, nypf2, :)  →  Python: dxdy[:, nypf2-1, :]
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
traj_fid = open('traj_xyz.txt', 'w')
traj_fid.write('iline it x y z\n')

# Which lines to trace (override with specific lines for testing)
LINES = list(range(1, nlines + 1))
LINES = [50, 75, 100, 125, 150]

for iline in LINES:
    # Pick starting point (1-based indices)
    xind = float(iline)
    xStart = psixy[iline - 1, jyomp - 1]  # psi value
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
    # STEP 3: Calculate puncture points and generate Poincare plot data
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
traj_fid.close()

print("\nField-line tracing complete.")
