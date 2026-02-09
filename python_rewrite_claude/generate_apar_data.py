"""
Generate and export apar and derivatives to NetCDF file.

This script loads the grid file and perturbed field data, computes
apar, dapardx, dapardy, and dapardz using the appropriate geometry
(shifted-circular or single-null), and saves them to a NetCDF file.

Usage:
    python generate_apar_data.py \\
        --grid-file grid.nc \\
        --apar-file apar_raw.npy \\
        --zperiod 5 \\
        --output apar_data.nc
"""

import os
import sys
import argparse
import numpy as np
import netCDF4

# Local module imports
from get_apar_sc import get_apar_sc
from get_apar_sn import get_apar_sn


def main():
    # --------------------------------------------------------------------------
    # Parse command-line arguments
    # --------------------------------------------------------------------------

    parser = argparse.ArgumentParser(
        description='Generate apar data and export to NetCDF')
    parser.add_argument('--grid-file', type=str, required=True,
                        help='BOUT++ grid file (.nc)')
    parser.add_argument('--apar-file', type=str, required=True,
                        help='Raw apar data file (.npy or .mat)')
    parser.add_argument('--zperiod', type=int, required=True,
                        help='Toroidal periodicity (e.g., 1 or 5)')
    parser.add_argument('--output', type=str, required=True,
                        help='Output NetCDF file path')
    args = parser.parse_args()

    gridfile = args.grid_file
    aparfile = args.apar_file
    zperiod = args.zperiod

    # --------------------------------------------------------------------------
    # Load raw apar data
    # --------------------------------------------------------------------------

    print(f"Loading apar data from {aparfile} ...")

    # Detect file format and load accordingly
    if aparfile.endswith('.npy'):
        apar_raw = np.load(aparfile)
    elif aparfile.endswith('.mat'):
        # Try scipy first, fall back to h5py for v7.3 files
        try:
            import scipy.io as sio
            apar_raw = sio.loadmat(aparfile)['apar']
        except NotImplementedError:
            import h5py
            with h5py.File(aparfile, 'r') as f:
                apar_raw = np.array(f['apar']).T  # HDF5 stores in transposed order
    else:
        print(f"Error: Unsupported apar file format: {aparfile}")
        print("Supported formats: .npy, .mat")
        sys.exit(1)

    print(f"  Apar shape = {apar_raw.shape}")

    # Get dimensions from apar_raw
    nx, ny, nz = apar_raw.shape
    print(f"  nx = {nx}, ny = {ny}, nz = {nz}")

    # --------------------------------------------------------------------------
    # Load grid information
    # --------------------------------------------------------------------------

    print(f"Loading grid information from {gridfile} ...")

    ds = netCDF4.Dataset(gridfile, 'r')

    # Validate grid dimensions match apar dimensions
    grid_vars = ds.variables
    grid_shape = grid_vars['Rxy'][:].shape
    if grid_shape != (nx, ny):
        print(f"\nError: Grid dimensions {grid_shape} don't match apar dimensions ({nx}, {ny})")
        print("Make sure the grid file and apar file are compatible.")
        ds.close()
        sys.exit(1)

    zShift = np.array(ds.variables['zShift'][:], dtype=float)
    psixy  = np.array(ds.variables['psixy'][:], dtype=float)
    bxy    = np.array(ds.variables['Bxy'][:], dtype=float)
    sinty  = np.array(ds.variables['sinty'][:], dtype=float)

    ixsep1 = int(ds.variables['ixseps1'][:])
    ixsep2 = int(ds.variables['ixseps2'][:])

    # Determine divertor configuration from grid file
    if ixsep2 < nx:
        divertor = 2  # double null
        print("  Double null configuration")
        nypf11 = int(ds.variables['jyseps1_1'][:]) + 1
        nypf21 = int(ds.variables['jyseps2_1'][:]) + 1
        nypf12 = int(ds.variables['jyseps1_2'][:]) + 1
        nypf22 = int(ds.variables['jyseps2_2'][:]) + 1

    elif ixsep1 < nx:
        divertor = 1
        print("  Single null configuration")
        ixsep = ixsep1
        nypf1 = int(ds.variables['jyseps1_1'][:]) + 1
        nypf2 = int(ds.variables['jyseps2_2'][:]) + 1

    else:
        divertor = 0
        ixsep = nx
        nypf1 = 0
        nypf2 = ny
        print("  Circular configuration")

    dy_raw = np.array(ds.variables['dy'][:], dtype=float)
    dy0 = float(dy_raw.flat[0])

    sa = np.array(ds.variables['ShiftAngle'][:], dtype=float).flatten()

    ds.close()

    dz = 2.0 * np.pi / zperiod / nz

    # --------------------------------------------------------------------------
    # Compute apar and derivatives
    # --------------------------------------------------------------------------

    print("Computing apar and derivatives ...")

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
    # Export to NetCDF file
    # --------------------------------------------------------------------------

    print(f"\nExporting apar data to {args.output} ...")

    # Create NetCDF file
    ncout = netCDF4.Dataset(args.output, 'w', format='NETCDF4')

    # Create dimensions
    ncout.createDimension('nx', nx)
    ncout.createDimension('ny', ny)
    ncout.createDimension('nz', nz)

    # Create variables
    var_apar = ncout.createVariable('apar', 'f8', ('nx', 'ny', 'nz'))
    var_dapardx = ncout.createVariable('dapardx', 'f8', ('nx', 'ny', 'nz'))
    var_dapardy = ncout.createVariable('dapardy', 'f8', ('nx', 'ny', 'nz'))
    var_dapardz = ncout.createVariable('dapardz', 'f8', ('nx', 'ny', 'nz'))

    # Add attributes
    var_apar.long_name = 'Parallel vector potential'
    var_dapardx.long_name = 'x-derivative of parallel vector potential'
    var_dapardy.long_name = 'y-derivative of parallel vector potential'
    var_dapardz.long_name = 'z-derivative of parallel vector potential'

    # Write data
    var_apar[:, :, :] = apar0
    var_dapardx[:, :, :] = dapardx0
    var_dapardy[:, :, :] = dapardy0
    var_dapardz[:, :, :] = dapardz0

    # Add global attributes
    ncout.description = 'Precomputed parallel vector potential and derivatives'
    ncout.divertor = divertor
    ncout.gridfile = gridfile
    ncout.aparfile = aparfile
    if divertor == 1:
        ncout.ixsep = ixsep
        ncout.nypf1 = nypf1
        ncout.nypf2 = nypf2
    elif divertor == 2:
        ncout.nypf11 = nypf11
        ncout.nypf21 = nypf21
        ncout.nypf12 = nypf12
        ncout.nypf22 = nypf22
    ncout.zperiod = zperiod
    ncout.nx = nx
    ncout.ny = ny
    ncout.nz = nz
    ncout.dz = dz
    ncout.dy0 = dy0

    ncout.close()

    print(f"\nSuccessfully exported apar data to {args.output}")
    print(f"  Shape: apar0 = {apar0.shape}")
    print(f"  Configuration: divertor = {divertor}")
    print("Done.")


if __name__ == '__main__':
    main()
