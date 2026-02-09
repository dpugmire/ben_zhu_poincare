"""
Collect puncture point data and generate Poincare plot.

Translated from MATLAB poincare.m (B. Zhu, 12/2023) to Python.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4


# --------------------------------------------------------------------------
# Step 0: Equilibrium and grid resolution
# --------------------------------------------------------------------------

gridfile = ("/Users/dpn/proj/bout++/poincare/boutpp_poincare/data/"
            "kstar_30306_7850_psi085105_nx260ny128_f2_v0.nc")

nx = 260; ny = 128; nz = 256; zperiod = 1
nlines = 13
direction = 1

# --------------------------------------------------------------------------
# Read grid info
# --------------------------------------------------------------------------

ds = netCDF4.Dataset(gridfile, 'r')
rxy   = np.array(ds.variables['Rxy'][:], dtype=float)
zxy   = np.array(ds.variables['Zxy'][:], dtype=float)
psixy = np.array(ds.variables['psixy'][:], dtype=float)
hthe  = np.array(ds.variables['hthe'][:], dtype=float)
nxsep = int(ds.variables['ixseps1'][:]) + 1  # 1-based for compatibility
psi_bndry = float(ds.variables['psi_bndry'][:])
psi_axis  = float(ds.variables['psi_axis'][:])
# psin at a representative y-column (MATLAB index 55 → Python 54)
psin = (psixy[:, 54] - psi_axis) / (psi_bndry - psi_axis)
ds.close()

nypf = 16

# Core boundary
corex = rxy[0, nypf:ny - nypf]
corey = zxy[0, nypf:ny - nypf]

# SOL boundary (concatenated segments)
solx_parts = [
    rxy[0, :nypf],
    rxy[0, ny - nypf:ny],
    rxy[:, -1],
    rxy[-1, ::-1],
    rxy[::-1, 0],
]
soly_parts = [
    zxy[0, :nypf],
    zxy[0, ny - nypf:ny],
    zxy[:, -1],
    zxy[-1, ::-1],
    zxy[::-1, 0],
]

# Separatrix
# MATLAB: tmp = [1:nypf ny-nypf:-1:nypf+1 ny-nypf+1:ny]
tmp = np.concatenate([
    np.arange(0, nypf),
    np.arange(ny - nypf - 1, nypf - 1, -1),
    np.arange(ny - nypf, ny)
])
# nxsep is 1-based: MATLAB rxy(nxsep-1, :) → Python rxy[nxsep-2, :]
sepx = 0.5 * (rxy[nxsep - 2, tmp] + rxy[nxsep - 1, tmp])
sepy = 0.5 * (zxy[nxsep - 2, tmp] + zxy[nxsep - 1, tmp])

# --------------------------------------------------------------------------
# Step 1: Generate Poincare plot
# --------------------------------------------------------------------------

cm = plt.cm.jet(np.linspace(0, 1, nlines))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

for iline in range(1, nlines + 1):
    print(f"\tradial index: {iline}")

    if direction == 1:
        filename = f"npz_pp/x{iline}y55z1_v3lc-01-250p.npz"
    elif direction == -1:
        filename = f"npz_pp/x{iline}y55z1_v3lc-01-250m.npz"

    if os.path.isfile(filename):
        data = np.load(filename)
        v2 = data['v2']  # py (R values)
        v3 = data['v3']  # pz (Z values)
        v4 = data['v4']  # ptheta
        v5 = data['v5']  # ppsi

        ax1.plot(v2, v3, '.', color=cm[iline - 1], markersize=2)
        ax2.plot(v5, v4, '.', color=cm[iline - 1], markersize=2)

ax1.plot(sepx, sepy, '--k')
ax1.set_aspect('equal')
ax1.set_xlim([1.2, 2.3])
ax1.set_xlabel('R (m)')
ax1.set_ylabel('Z (m)')

# Reference psi line at MATLAB index 195 → Python 194
ax2.axvline(x=psixy[194, 54], color='k', linestyle='--')
ax2.set_xlabel(r'$\psi$')
ax2.set_ylabel(r'$\theta/\pi$')

plt.tight_layout()
plt.savefig('poincare_plot.png', dpi=150, bbox_inches='tight')
plt.show()

print("Poincare plot generated.")
