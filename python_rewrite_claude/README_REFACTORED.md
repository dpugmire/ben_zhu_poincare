# Refactored Field-Line Tracing Code

This directory contains a refactored version of the field-line tracing code, split into two separate scripts for better modularity and performance.

## Overview

The code has been split into two stages:

1. **Data Generation** (`generate_apar_data.py`): Computes apar and derivatives (expensive, run once)
2. **Field-Line Tracing** (`trace_fieldlines.py`): Performs the actual particle tracing (can be run multiple times)

## Workflow

### Step 1: Generate apar data (run once)

```bash
# For single-null configuration
python generate_apar_data.py \
    --grid-file /path/to/grid.nc \
    --apar-file /path/to/apar_raw.mat \
    --zperiod 1 \
    --output apar_data.nc

# For circular configuration
python generate_apar_data.py \
    --grid-file /path/to/circular_grid.nc \
    --apar-file /path/to/apar_raw.npy \
    --zperiod 5 \
    --output apar_data_circular.nc
```

This will:

- Load the grid file and raw perturbed field data (.npy or .mat)
- Automatically detect the configuration type (circular/single-null/double-null) from the grid file
- Compute `apar`, `dapardx`, `dapardy`, `dapardz` using `get_apar_sc` or `get_apar_sn`
- Save the results to a NetCDF file
- Exit

**Note**: This step is computationally expensive but only needs to be run once. The resulting NetCDF file can be reused for multiple tracing runs with different parameters.

### Step 2: Trace field lines (run multiple times)

```bash
# Basic usage - trace default lines
python trace_fieldlines.py \
    --apar-file apar_data.nc \
    --grid-file /path/to/grid.nc

# Trace specific lines
python trace_fieldlines.py \
    --apar-file apar_data.nc \
    --grid-file /path/to/grid.nc \
    --lines "50,75,100,125,150"

# Customize tracing parameters
python trace_fieldlines.py \
    --apar-file apar_data.nc \
    --grid-file /path/to/grid.nc \
    --direction 1 \
    --nlines 256 \
    --nturns 25
```

This will:

- Read the precomputed apar data from the NetCDF file
- Load the grid information
- Perform field-line tracing
- Save puncture points to `ip_xyz.txt`
- Save trajectories to `traj_xyz.txt`
- Save detailed data to `npz_pp/` and `npz_traj/` directories

## Output Files

### Text Files

- **`ip_xyz.txt`**: Puncture points (intersections with the poloidal plane)
  - Format: `iline it ipx ipy ipz`
- **`traj_xyz.txt`**: Complete 3D trajectories
  - Format: `iline it x y z`

### NPZ Files

- **`npz_pp/`**: Puncture point data for each field line
  - Naming: `x{iline}y{yyy}z{zzz}_v3lc-00-250{p/m}.npz`
  - Contains: positions, theta, psi, connection length, region
- **`npz_traj/`**: Trajectory data for each field line
  - Naming: `x{iline}y{yyy}z{zzz}_v3lc-01-250{p/m}.npz`
  - Contains: trajectory in index coordinates

## Command-Line Arguments

### `generate_apar_data.py`

| Argument      | Required | Default | Description                         |
| ------------- | -------- | ------- | ----------------------------------- |
| `--grid-file` | Yes      | -       | BOUT++ grid file (.nc)              |
| `--apar-file` | Yes      | -       | Raw apar data file (.npy or .mat)   |
| `--zperiod`   | Yes      | -       | Toroidal periodicity (e.g., 1 or 5) |
| `--output`    | Yes      | -       | Output NetCDF file path             |

**Note**: The configuration type (circular/single-null/double-null) is automatically detected from the grid file.

### `trace_fieldlines.py`

| Argument      | Required | Default | Description                                                |
| ------------- | -------- | ------- | ---------------------------------------------------------- |
| `--apar-file` | Yes      | -       | Input NetCDF file with apar data                           |
| `--grid-file` | Yes      | -       | BOUT++ grid file                                           |
| `--direction` | No       | 1       | Tracing direction: 1 (forward) or -1 (backward)            |
| `--nlines`    | No       | 256     | Number of field lines to trace                             |
| `--nturns`    | No       | 25      | Approximate number of poloidal turns                       |
| `--lines`     | No       | -       | Comma-separated list of specific lines (e.g., "50,75,100") |

## Benefits of This Refactoring

1. **Performance**: The expensive apar computation is done once and cached
2. **Flexibility**: Can run multiple tracing jobs with different parameters without recomputing apar
3. **Modularity**: Clear separation between data generation and analysis
4. **Debugging**: Easier to debug each stage independently
5. **C++ Integration**: The NetCDF format makes it easy to read precomputed data from C++

## NetCDF File Contents

The `apar_data.nc` file contains:

### Dimensions

- `nx`: Radial grid points
- `ny`: Poloidal grid points
- `nz`: Toroidal grid points

### Variables

- `apar(nx, ny, nz)`: Parallel vector potential
- `dapardx(nx, ny, nz)`: x-derivative of apar
- `dapardy(nx, ny, nz)`: y-derivative of apar
- `dapardz(nx, ny, nz)`: z-derivative of apar

### Global Attributes

- `divertor_case`: Configuration type
- `divertor`: Divertor flag (0=circular, 1=single-null, 2=double-null)
- `gridfile`: Path to original grid file
- `aparfile`: Path to original apar file
- `ixsep`, `nypf1`, `nypf2`: Grid parameters (for single-null)
- `zperiod`: Toroidal periodicity
- `nx`, `ny`, `nz`: Grid dimensions
- `dz`, `dy0`: Grid spacings

## Example Workflow

```bash
# 1. Generate apar data (expensive, run once)
python generate_apar_data.py \
    --grid-file /path/to/grid.nc \
    --apar-file /path/to/apar_raw.mat \
    --zperiod 1 \
    --output apar_data.nc

# 2. Trace a few test lines
python trace_fieldlines.py \
    --apar-file apar_data.nc \
    --grid-file /path/to/grid.nc \
    --lines "50,100,150"

# 3. Trace all lines for production run
python trace_fieldlines.py \
    --apar-file apar_data.nc \
    --grid-file /path/to/grid.nc \
    --nlines 256 \
    --nturns 25

# 4. Run backward tracing
python trace_fieldlines.py \
    --apar-file apar_data.nc \
    --grid-file /path/to/grid.nc \
    --direction -1 \
    --lines "50,100,150"
```

## Dependencies

- Python 3.6+
- numpy
- scipy
- netCDF4
- h5py (for MATLAB v7.3 files)

###

# circular case

python generate_apar_data.py --grid-file /Users/dpn/proj/bout++/nersc_data/circle-zonal/cbm18_dens3_0.5BS_516nx64ny.grid.nc --apar-file /Users/dpn/proj/bout++/nersc_data/circle-zonal/apar_cbm18_dens3_0.5BS_516nx64ny64nz_t500.npy --zperiod 5 --output apar_circular.nc

python trace_fieldlines.py --apar-file apar_circular.nc --grid-file /Users/dpn/proj/bout++/nersc_data/circle-zonal/cbm18_dens3_0.5BS_516nx64ny.grid.nc --lines "50,75,100,125,150"

# single null case

python generate_apar_data.py --grid-file /Users/dpn/proj/bout++/poincare/boutpp_poincare/data/kstar_30306_7850_psi085105_nx260ny128_f2_v0.nc --apar-file /Users/dpn/proj/bout++/nersc_data/xpoint_singlenull/apar_kstar_30306_7850_psi085105_nx260ny128_f2_nz256.mat --zperiod 5 --output apar_single.nc

python trace_fieldlines.py --apar-file apar_single.nc --grid-file /Users/dpn/proj/bout++/poincare/boutpp_poincare/data/kstar_30306_7850_psi085105_nx260ny128_f2_v0.nc --lines "50,75,100,125,150"
