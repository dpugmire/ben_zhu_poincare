# C++ field-line tracer for `fieldline_tracing.m`

This directory contains a C++ implementation of the MATLAB tracing workflow in
`/Users/dpn/proj/bout++/ben_zhu_poincare/codex_cxx/matlab/fieldline_tracing.m`,
using precomputed `apar` derivatives from NetCDF.

## MATLAB-to-C++ scope

- `get_apar_sc.m` / `get_apar_sn.m` are **not** recomputed in C++.
- Instead, C++ reads the precomputed arrays from:
  - `apar_circular.nc` for `divertorCase = 0`
  - `apar_single.nc` for `divertorCase = 1`
- Required variables in NetCDF are:
  - `apar(nx,ny,nz)`
  - `dapardx(nx,ny,nz)`
  - `dapardy(nx,ny,nz)`
  - `dapardz(nx,ny,nz)`

The tracer reproduces the RK4 stepping, region logic (CFR/SOL/PFR/divertor),
twist-shift handling, and puncture extraction to `ip_xyz.txt`.

## Build

```bash
cd /Users/dpn/proj/bout++/ben_zhu_poincare/codex_cxx
make
```

## Run

### Single-null (`divertorCase = 1`)

```bash
./trace_fieldlines --divertor-case 1 --lines 151 --nturns 100
```

### Circular (`divertorCase = 0`)

```bash
./trace_fieldlines --divertor-case 0 --lines 151 --nturns 100
```

You can also pass explicit paths:

```bash
./trace_fieldlines --apar-file /path/to/apar.nc --grid-file /path/to/grid.nc --lines 151 --nturns 100
```

## Outputs

- `ip_xyz.txt`
- `ip_thetapsi.txt`
- `traj_xyz.txt`

Reference MATLAB puncture outputs are in:

- `/Users/dpn/proj/bout++/ben_zhu_poincare/codex_cxx/matlab/matlab_ip_xyz.circular.txt`
- `/Users/dpn/proj/bout++/ben_zhu_poincare/codex_cxx/matlab/matlab_ip_xyz.single.txt`
