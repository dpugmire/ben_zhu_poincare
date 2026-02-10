# Python Conversion of MATLAB Field-Line Tools

This directory contains a Python conversion of:

- `fieldline_tracing.m`
- `RK4_FLT1.m`
- `get_apar_sc.m`
- `get_apar_sn.m`
- `parsave.m`
- `parsave8.m`
- `poincare.m`

## Files

- `fieldline_tracing.py`: main tracer (`run_fieldline_tracing`) and CLI
- `rk4_flt1.py`: RK4 stepping routine
- `get_apar_sc.py`: circular (`divertor=0`) Apar and derivative prep
- `get_apar_sn.py`: single-null (`divertor=1`) Apar and derivative prep
- `numerics.py`: interpolation helpers used to mimic MATLAB behavior
- `io_helpers.py`: MAT save helpers (`parsave`, `parsave8`)
- `poincare.py`: plotting helper equivalent to `poincare.m`

## Dependencies

- Python 3.9+
- `numpy`
- `scipy`
- `netCDF4`
- `matplotlib` (for `poincare.py`)

## Run Field-Line Tracing

Single-null (`divertor=1`, KSTAR defaults):

```bash
cd /Users/dpn/proj/bout++/ben_zhu_poincare/convert
python3 fieldline_tracing.py --case 1 --output-dir .
```

Circular (`divertor=0`, circle-zonal defaults):

```bash
cd /Users/dpn/proj/bout++/ben_zhu_poincare/convert
python3 fieldline_tracing.py --case 0 --output-dir .
```

Useful overrides:

- `--gridfile ...`
- `--aparfile ...` (`.mat` and `.npy` are both supported)
- `--apar-variable apar` (MAT variable name)
- `--time-index -1` (for 4D Apar arrays; default is last slice)
- `--nlines`, `--nturns`, `--direction`
- `--lines "50, 100"` (trace only specified radial indices; overrides `--nlines`)
- `--equilibrium` (set `dxdy=dzdy=0`, matching MATLAB equilibrium quick-test behavior)

## Run Poincare Plot Script

```bash
cd /Users/dpn/proj/bout++/ben_zhu_poincare/convert
python3 poincare.py --pp-dir ./mat_pp
```

## Notes

- The script checks that the requested `divertor` mode matches what is inferred from the grid file.
- Double-null (`divertor=2`) is still intentionally not implemented, matching the MATLAB script status.
