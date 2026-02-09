"""
Quick test to demonstrate dimension validation.
"""
import netCDF4
import numpy as np

# Create a fake apar file with dimensions (100, 50, 64)
print("Creating test apar file with dimensions (100, 50, 64)...")
with netCDF4.Dataset('test_apar.nc', 'w') as nc:
    nc.createDimension('nx', 100)
    nc.createDimension('ny', 50)
    nc.createDimension('nz', 64)

    var = nc.createVariable('apar', 'f8', ('nx', 'ny', 'nz'))
    var[:] = np.zeros((100, 50, 64))

    var = nc.createVariable('dapardx', 'f8', ('nx', 'ny', 'nz'))
    var[:] = np.zeros((100, 50, 64))

    var = nc.createVariable('dapardy', 'f8', ('nx', 'ny', 'nz'))
    var[:] = np.zeros((100, 50, 64))

    var = nc.createVariable('dapardz', 'f8', ('nx', 'ny', 'nz'))
    var[:] = np.zeros((100, 50, 64))

    nc.nx = 100
    nc.ny = 50
    nc.nz = 64
    nc.zperiod = 1
    nc.divertor = 0
    nc.gridfile = "/path/to/original/grid_100x50.nc"

# Create a fake grid file with dimensions (200, 80)
print("Creating test grid file with dimensions (200, 80)...")
with netCDF4.Dataset('test_grid.nc', 'w') as nc:
    nc.createDimension('x', 200)
    nc.createDimension('y', 80)

    var = nc.createVariable('Rxy', 'f8', ('x', 'y'))
    var[:] = np.ones((200, 80))

    var = nc.createVariable('Zxy', 'f8', ('x', 'y'))
    var[:] = np.ones((200, 80))

    var = nc.createVariable('psixy', 'f8', ('x', 'y'))
    var[:] = np.ones((200, 80))

    var = nc.createVariable('zShift', 'f8', ('x', 'y'))
    var[:] = np.zeros((200, 80))

    var = nc.createVariable('Bxy', 'f8', ('x', 'y'))
    var[:] = np.ones((200, 80))

    var = nc.createVariable('Btxy', 'f8', ('x', 'y'))
    var[:] = np.ones((200, 80))

    var = nc.createVariable('Bpxy', 'f8', ('x', 'y'))
    var[:] = np.ones((200, 80))

    var = nc.createVariable('hthe', 'f8', ('x', 'y'))
    var[:] = np.ones((200, 80))

    var = nc.createVariable('sinty', 'f8', ('x', 'y'))
    var[:] = np.zeros((200, 80))

    var = nc.createVariable('bxcvx', 'f8', ('x', 'y'))
    var[:] = np.zeros((200, 80))

    var = nc.createVariable('bxcvy', 'f8', ('x', 'y'))
    var[:] = np.zeros((200, 80))

    var = nc.createVariable('bxcvz', 'f8', ('x', 'y'))
    var[:] = np.zeros((200, 80))

    var = nc.createVariable('Jpar0', 'f8', ('x', 'y'))
    var[:] = np.zeros((200, 80))

    var = nc.createVariable('dy', 'f8', ('x', 'y'))
    var[:] = np.ones((200, 80)) * 0.1

    var = nc.createVariable('ShiftAngle', 'f8', ('x',))
    var[:] = np.zeros(200)

    var = nc.createVariable('rmag', 'f8', ())
    var[:] = 1.0

    var = nc.createVariable('ixseps1', 'i', ())
    var[:] = 300  # > nx means circular

    var = nc.createVariable('ixseps2', 'i', ())
    var[:] = 300

print("\nNow testing dimension validation...")
print("Running: python trace_fieldlines.py --apar-file test_apar.nc --grid-file test_grid.nc --lines \"1\"\n")
