from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import scipy.io as sio

try:
    from .fieldline_tracing import DEFAULT_GRID
except ImportError:
    from fieldline_tracing import DEFAULT_GRID


def _read_scalar(ds: netCDF4.Dataset, name: str) -> float:
    return float(np.asarray(ds.variables[name][:]).squeeze())


def _read_2d(ds: netCDF4.Dataset, name: str, nx: int, ny: int) -> np.ndarray:
    arr = np.asarray(ds.variables[name][:], dtype=float)
    if arr.shape == (ny, nx):
        return arr.T.copy()
    if arr.shape == (nx, ny):
        return arr.copy()
    raise ValueError(f"Unexpected shape for {name}: {arr.shape}")


def build_poincare_plot(
    gridfile: str,
    pp_dir: str,
    output_png: str | None,
    direction: int,
    nlines: int,
    y_start: int,
    z_start: int,
    nypf: int,
    nx: int,
    ny: int,
) -> None:
    with netCDF4.Dataset(gridfile, "r") as ds:
        rxy = _read_2d(ds, "Rxy", nx, ny)
        zxy = _read_2d(ds, "Zxy", nx, ny)
        psixy = _read_2d(ds, "psixy", nx, ny)
        ixsep = int(_read_scalar(ds, "ixseps1")) + 1

    tmp = np.concatenate(
        [
            np.arange(1, nypf + 1),
            np.arange(ny - nypf, nypf, -1),
            np.arange(ny - nypf + 1, ny + 1),
        ]
    )
    tmp0 = tmp - 1
    sepx = 0.5 * (rxy[ixsep - 2, tmp0] + rxy[ixsep - 1, tmp0])
    sepy = 0.5 * (zxy[ixsep - 2, tmp0] + zxy[ixsep - 1, tmp0])

    color_map = plt.cm.jet(np.linspace(0.0, 1.0, nlines))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    tail = "p" if direction == 1 else "m"
    pp_path = Path(pp_dir)
    for iline in range(1, nlines + 1):
        mat_file = pp_path / f"x{iline}y{y_start}z{z_start}_v3lc-01-250{tail}.mat"
        if not mat_file.exists():
            continue

        data = sio.loadmat(mat_file)
        v2 = np.ravel(data.get("v2", np.array([])))
        v3 = np.ravel(data.get("v3", np.array([])))
        v4 = np.ravel(data.get("v4", np.array([])))
        v5 = np.ravel(data.get("v5", np.array([])))
        if v2.size and v3.size:
            ax1.plot(v2, v3, ".", color=color_map[iline - 1], markersize=2)
        if v5.size and v4.size:
            ax2.plot(v5, v4, ".", color=color_map[iline - 1], markersize=2)

    ax1.plot(sepx, sepy, "--k")
    ax1.set_aspect("equal", adjustable="box")
    ax1.set_xlim([1.2, 2.3])
    ax1.set_xlabel("R(m)")
    ax1.set_ylabel("Z(m)")

    sep_psi = float(psixy[194, 54])
    ax2.plot([sep_psi, sep_psi], [0.0, 2.0], "--k")
    ax2.set_xlabel("psi")
    ax2.set_ylabel("theta/pi")

    fig.tight_layout()
    if output_png:
        fig.savefig(output_png, dpi=150)
    else:
        plt.show()


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Python port of poincare.m")
    parser.add_argument("--gridfile", default=DEFAULT_GRID)
    parser.add_argument("--pp-dir", default="./mat_pp")
    parser.add_argument("--output-png")
    parser.add_argument("--direction", type=int, default=1)
    parser.add_argument("--nlines", type=int, default=13)
    parser.add_argument("--y-start", type=int, default=55)
    parser.add_argument("--z-start", type=int, default=1)
    parser.add_argument("--nypf", type=int, default=16)
    parser.add_argument("--nx", type=int, default=260)
    parser.add_argument("--ny", type=int, default=128)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    build_poincare_plot(
        gridfile=args.gridfile,
        pp_dir=args.pp_dir,
        output_png=args.output_png,
        direction=args.direction,
        nlines=args.nlines,
        y_start=args.y_start,
        z_start=args.z_start,
        nypf=args.nypf,
        nx=args.nx,
        ny=args.ny,
    )


if __name__ == "__main__":
    main()
