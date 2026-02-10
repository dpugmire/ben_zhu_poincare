# C++ Field-Line Tracing Prototype

This directory contains the first-pass C++ conversion of `trace_fieldlines.py`.
This was created with openAI Codex.

## Design choices

- All multidimensional arrays are flat `std::vector<float>`.
- `FieldData` (`field_data.h/.cpp`) owns all grid/apar/field arrays and interpolation helpers.
- RK4 stepping is in `rk4.h/.cpp`.
- The main particle tracing and puncture loop is a standalone function:
  `trace_field_lines(...)` in `trace.cpp`.
- `generate_apar_data.py` remains Python and is used to produce the NetCDF apar input.

## Build

```bash
cd /Users/dpn/proj/bout++/ben_zhu_poincare/python_rewrite_claude/cxx
make
```

## Run

```bash
./trace_fieldlines \
  --apar-file /path/to/apar_data.nc \
  --grid-file /path/to/grid.nc \
  --lines "50,75,100"
```

Outputs are written to:

- `ip_xyz.txt`
- `ip_thetapsi.txt`
- `traj_xyz.txt`
