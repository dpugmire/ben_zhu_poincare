from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import netCDF4
import numpy as np
import scipy.io as sio
from scipy.interpolate import CubicSpline

try:
    from .get_apar_sc import get_apar_sc
    from .get_apar_sn import get_apar_sn
    from .io_helpers import parsave, parsave8
    from .numerics import (
        clamp_index_1based,
        coerce_apar_shape,
        interp1_linear,
        interp1_spline,
        interp2_spline,
    )
    from .rk4_flt1 import rk4_flt1
except ImportError:
    from get_apar_sc import get_apar_sc
    from get_apar_sn import get_apar_sn
    from io_helpers import parsave, parsave8
    from numerics import (
        clamp_index_1based,
        coerce_apar_shape,
        interp1_linear,
        interp1_spline,
        interp2_spline,
    )
    from rk4_flt1 import rk4_flt1

TWOPI = 2.0 * np.pi


@dataclass
class FieldlineTracingConfig:
    gridfile: str
    aparfile: str
    nx: int = 260
    ny: int = 128
    nz: int = 256
    zperiod: int = 1
    divertor: int = 1
    direction: int = 1
    nlines: int = 256
    nturns: int = 250
    np_max: int = 1250
    save_traj: bool = True
    save_pp: bool = True
    output_dir: str = "."
    apar_variable: str = "apar"
    time_index: int = -1
    lines: tuple[int, ...] | None = None
    equilibrium: bool = False


DEFAULT_GRID = "/Users/dpn/proj/bout++/poincare/boutpp_poincare/data/kstar_30306_7850_psi085105_nx260ny128_f2_v0.nc"
DEFAULT_APAR = "/Users/dpn/proj/bout++/poincare/boutpp_poincare/data/apar_kstar_30306_7850_psi085105_nx260ny128_f2_nz256.mat"
CASE0_GRID = "/Users/dpn/proj/bout++/nersc_data/circle-zonal/cbm18_dens3_0.5BS_516nx64ny.grid.nc"
CASE0_APAR = "/Users/dpn/proj/bout++/nersc_data/circle-zonal/apar_cbm18_dens3_0.5BS_516nx64ny64nz_t500.npy"


def _case_defaults(case: int) -> dict[str, Any]:
    if case == 0:
        return {
            "gridfile": CASE0_GRID,
            "aparfile": CASE0_APAR,
            "nx": 516,
            "ny": 64,
            "nz": 64,
            "zperiod": 5,
            "divertor": 0,
        }
    if case == 1:
        return {
            "gridfile": DEFAULT_GRID,
            "aparfile": DEFAULT_APAR,
            "nx": 260,
            "ny": 128,
            "nz": 256,
            "zperiod": 1,
            "divertor": 1,
        }
    raise ValueError(f"Unsupported case={case}, expected 0 or 1")


def _read_scalar(ds: netCDF4.Dataset, name: str) -> float:
    return float(np.asarray(ds.variables[name][:]).squeeze())


def _read_2d(ds: netCDF4.Dataset, name: str, nx: int, ny: int) -> np.ndarray:
    arr = np.asarray(ds.variables[name][:], dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be 2D, got shape={arr.shape}")
    if arr.shape == (ny, nx):
        return arr.T.copy()
    if arr.shape == (nx, ny):
        return arr.copy()
    raise ValueError(f"Unexpected shape for {name}: {arr.shape}, expected {(ny, nx)} or {(nx, ny)}")


def _load_grid(config: FieldlineTracingConfig):
    nx = config.nx
    ny = config.ny

    with netCDF4.Dataset(config.gridfile, "r") as ds:
        z_shift = _read_2d(ds, "zShift", nx, ny)
        rxy = _read_2d(ds, "Rxy", nx, ny)
        zxy = _read_2d(ds, "Zxy", nx, ny)
        psixy = _read_2d(ds, "psixy", nx, ny)
        rmag = _read_scalar(ds, "rmag")

        ixsep1 = int(_read_scalar(ds, "ixseps1"))
        ixsep2 = int(_read_scalar(ds, "ixseps2"))

        if ixsep2 < nx:
            divertor = 2
            nypf11 = int(_read_scalar(ds, "jyseps1_1")) + 1
            nypf21 = int(_read_scalar(ds, "jyseps2_1")) + 1
            nypf12 = int(_read_scalar(ds, "jyseps1_2")) + 1
            nypf22 = int(_read_scalar(ds, "jyseps2_2")) + 1
            _ = (nypf11, nypf21, nypf12, nypf22)
            ixsep = ixsep1
            nypf1 = 0
            nypf2 = ny
        elif ixsep1 < nx:
            divertor = 1
            ixsep = ixsep1
            nypf1 = int(_read_scalar(ds, "jyseps1_1")) + 1
            nypf2 = int(_read_scalar(ds, "jyseps2_2")) + 1
        else:
            divertor = 0
            ixsep = nx
            nypf1 = 0
            nypf2 = ny

        bxy = _read_2d(ds, "Bxy", nx, ny)
        btxy = _read_2d(ds, "Btxy", nx, ny)
        bpxy = _read_2d(ds, "Bpxy", nx, ny)
        hthe = _read_2d(ds, "hthe", nx, ny)
        sinty = _read_2d(ds, "sinty", nx, ny)

        bxcvx = _read_2d(ds, "bxcvx", nx, ny)
        bxcvy = _read_2d(ds, "bxcvy", nx, ny)
        bxcvz = _read_2d(ds, "bxcvz", nx, ny)

        jpar0 = _read_2d(ds, "Jpar0", nx, ny)
        dy = np.asarray(ds.variables["dy"][:], dtype=float)
        dy0 = float(dy.flat[0])

        sa = np.asarray(ds.variables["ShiftAngle"][:], dtype=float).squeeze()

    return {
        "z_shift": z_shift,
        "rxy": rxy,
        "zxy": zxy,
        "psixy": psixy,
        "rmag": rmag,
        "ixsep1": ixsep1,
        "ixsep2": ixsep2,
        "ixsep": ixsep,
        "nypf1": nypf1,
        "nypf2": nypf2,
        "divertor": divertor,
        "bxy": bxy,
        "btxy": btxy,
        "bpxy": bpxy,
        "hthe": hthe,
        "sinty": sinty,
        "bxcvx": bxcvx,
        "bxcvy": bxcvy,
        "bxcvz": bxcvz,
        "jpar0": jpar0,
        "dy0": dy0,
        "sa": sa,
    }


def _theta_profile(divertor: int, rxy: np.ndarray, zxy: np.ndarray, nypf1: int, nypf2: int, ixsep: int, ny: int):
    theta = np.zeros(ny, dtype=float)

    if divertor == 1:
        center_x = 0.5 * (np.max(rxy[0, nypf1 : ny - nypf1]) + np.min(rxy[0, nypf1 : ny - nypf1]))
        center_y = 0.5 * (np.max(zxy[0, nypf1 : ny - nypf1]) + np.min(zxy[0, nypf1 : ny - nypf1]))

        tmp = np.concatenate(
            [
                np.arange(1, nypf1 + 1),
                np.arange(ny - nypf1, nypf1, -1),
                np.arange(ny - nypf1 + 1, ny + 1),
            ]
        )
        tmp0 = tmp - 1
        sepx = 0.5 * (rxy[ixsep - 1, tmp0] + rxy[ixsep, tmp0])
        sepy = 0.5 * (zxy[ixsep - 1, tmp0] + zxy[ixsep, tmp0])

        xpoint_x = 0.25 * (
            sepx[nypf1 - 1] + sepx[nypf1] + sepx[ny - nypf1 - 1] + sepx[ny - nypf1]
        )
        xpoint_y = 0.25 * (
            sepy[nypf1 - 1] + sepy[nypf1] + sepy[ny - nypf1 - 1] + sepy[ny - nypf1]
        )
        u = np.array([center_x - xpoint_x, center_y - xpoint_y, 0.0])

        for iy in range(ny):
            v = np.array([center_x - rxy[0, iy], center_y - zxy[0, iy], 0.0])
            theta[iy] = np.arctan2(np.linalg.norm(np.cross(u, v)), np.dot(u, v))

        theta = theta / np.pi
        itheta = int(np.argmax(theta))
        theta[itheta:] = 2.0 - theta[itheta:]
        itheta = int(np.argmax(theta))
        if itheta != ny - 1:
            theta[itheta:] = 4.0 - theta[itheta:]

        theta = theta - theta[nypf1]

    elif divertor == 0:
        center_x = 0.5 * (np.max(rxy[0, :]) + np.min(rxy[0, :]))
        center_y = 0.5 * (np.max(zxy[0, :]) + np.min(zxy[0, :]))
        u = np.array([center_x - rxy[0, 0], center_y - zxy[0, 0], 0.0])

        for iy in range(ny):
            v = np.array([center_x - rxy[0, iy], center_y - zxy[0, iy], 0.0])
            theta[iy] = np.arctan2(np.linalg.norm(np.cross(u, v)), np.dot(u, v))

        theta = theta / np.pi
        itheta = int(np.argmax(theta))
        theta[itheta:] = 2.0 - theta[itheta:]
        itheta = int(np.argmax(theta))
        if itheta != ny - 1:
            theta[itheta:] = 4.0 - theta[itheta:]

    else:
        raise NotImplementedError("Double-null theta profile is not implemented")

    return theta


def _find_crossings(fl_x3d: np.ndarray, quant_step: float = 1.0e-4) -> list[float]:
    """Find x=0 crossings following MATLAB spline-based puncture logic."""
    n = fl_x3d.size
    if n < 2:
        return []

    xk = np.arange(1, n + 1, dtype=float)
    cs = CubicSpline(xk, fl_x3d, bc_type="not-a-knot", extrapolate=True)

    roots: list[float] = []
    for i in range(n - 1):
        # Polynomial on interval i: c0*t^3 + c1*t^2 + c2*t + c3, t in [0, 1]
        c0, c1, c2, c3 = (float(cs.c[0, i]), float(cs.c[1, i]), float(cs.c[2, i]), float(cs.c[3, i]))
        poly_roots = np.roots([c0, c1, c2, c3])
        x0 = float(i + 1)
        for rr in poly_roots:
            if abs(rr.imag) > 1.0e-10:
                continue
            t = float(rr.real)
            # Include right endpoint only for the last interval to avoid duplicates.
            if (0.0 <= t < 1.0) or (i == n - 2 and np.isclose(t, 1.0, atol=1.0e-12)):
                roots.append(x0 + t)

    if not roots:
        return []

    roots = sorted(roots)
    dedup: list[float] = []
    for r in roots:
        if not dedup or abs(r - dedup[-1]) > 1.0e-10:
            dedup.append(r)

    if quant_step > 0.0:
        quantized = [np.floor(r / quant_step) * quant_step for r in dedup]
        # Deduplicate again after quantization.
        out: list[float] = []
        for r in quantized:
            if not out or abs(r - out[-1]) > 0.5 * quant_step:
                out.append(float(r))
        return out

    return [float(r) for r in dedup]


def _infer_time_axis(shape: tuple[int, ...], nx: int, ny: int, nz: int) -> int:
    target = sorted((nx, ny, nz))
    for axis in range(len(shape)):
        rem = shape[:axis] + shape[axis + 1 :]
        if sorted(rem) == target:
            return axis
    return 0


def _load_apar_array(
    aparfile: str,
    nx: int,
    ny: int,
    nz: int,
    apar_variable: str,
    time_index: int,
) -> np.ndarray:
    path = Path(aparfile)
    suffix = path.suffix.lower()

    if suffix == ".npy":
        raw = np.load(path, mmap_mode="r")
    elif suffix == ".mat":
        data = sio.loadmat(path)
        if apar_variable in data and isinstance(data[apar_variable], np.ndarray):
            raw = data[apar_variable]
        else:
            candidates = [
                key
                for key, value in data.items()
                if not key.startswith("__") and isinstance(value, np.ndarray) and value.ndim in (3, 4)
            ]
            if not candidates:
                raise KeyError(
                    f"No 3D/4D array found in {aparfile}. Tried variable '{apar_variable}' first."
                )
            raw = data[candidates[0]]
    else:
        raise ValueError(f"Unsupported Apar file extension: {suffix}")

    raw = np.asarray(raw)
    if raw.ndim == 4:
        t_axis = _infer_time_axis(raw.shape, nx, ny, nz)
        tidx = time_index if time_index >= 0 else raw.shape[t_axis] + time_index
        if tidx < 0 or tidx >= raw.shape[t_axis]:
            raise IndexError(
                f"time_index={time_index} is out of bounds for axis size {raw.shape[t_axis]}"
            )
        raw = np.take(raw, tidx, axis=t_axis)

    if raw.ndim != 3:
        raise ValueError(f"Expected a 3D Apar array after loading, got shape={raw.shape}")

    return coerce_apar_shape(raw, nx, ny, nz)


def _parse_line_selection(lines_arg: str | None) -> tuple[int, ...] | None:
    if lines_arg is None:
        return None

    raw_tokens = [tok.strip() for tok in lines_arg.split(",")]
    values: list[int] = []
    for tok in raw_tokens:
        if not tok:
            continue
        try:
            values.append(int(tok))
        except ValueError as exc:
            raise ValueError(f"Invalid line index '{tok}' in --lines") from exc

    if not values:
        raise ValueError("--lines provided but no valid indices were found")

    seen: set[int] = set()
    ordered_unique: list[int] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            ordered_unique.append(value)
    return tuple(ordered_unique)


def run_fieldline_tracing(config: FieldlineTracingConfig) -> None:
    nx = config.nx
    ny = config.ny
    nz = config.nz
    zperiod = config.zperiod
    direction = config.direction
    if direction not in (-1, 1):
        raise ValueError(f"direction must be +/-1, got {direction}")
    nlines = min(config.nlines, nx)
    nturns = config.nturns

    nzg = nz * zperiod
    zmin = 0.0
    zmax = TWOPI
    dz = (zmax - zmin) / nzg
    ziarray = np.arange(1, nzg + 2, dtype=float)
    zarray = (ziarray - 1.0) * dz

    xiarray = np.arange(1, nx + 1, dtype=float)
    yiarray = np.arange(1, ny + 1, dtype=float)

    print("Loading grid information ...")
    grid = _load_grid(config)

    z_shift = grid["z_shift"]
    rxy = grid["rxy"]
    zxy = grid["zxy"]
    psixy = grid["psixy"]

    ixsep1 = int(grid["ixsep1"])
    ixsep = int(grid["ixsep"])
    nypf1 = int(grid["nypf1"])
    nypf2 = int(grid["nypf2"])
    divertor = int(grid["divertor"])
    if divertor != int(config.divertor):
        raise ValueError(
            f"Grid-derived divertor={divertor} does not match requested divertor={config.divertor}."
        )

    bxy = grid["bxy"]
    btxy = grid["btxy"]
    bpxy = grid["bpxy"]
    hthe = grid["hthe"]
    sinty = grid["sinty"]
    bxcvx = grid["bxcvx"]
    bxcvy = grid["bxcvy"]
    bxcvz = grid["bxcvz"]
    jpar0 = grid["jpar0"]
    dy0 = float(grid["dy0"])
    sa = np.asarray(grid["sa"], dtype=float).squeeze()

    jyomp = int(np.argmax(rxy[-1, :])) + 1
    xarray = psixy[:, jyomp - 1]
    x_min = float(np.min(xarray))
    x_max = float(np.max(xarray))

    dz = TWOPI / zperiod / nz
    nu = btxy * hthe / bpxy / rxy

    theta = _theta_profile(divertor, rxy, zxy, nypf1, nypf2, ixsep, ny)

    if divertor == 2:
        raise NotImplementedError("Double-null configuration is not implemented in this Python port")

    xiarray_cfr = np.arange(1, ixsep + 1, dtype=float)
    if divertor == 0:
        yiarray_cfr = np.arange(1, ny + 2, dtype=float)
        theta_cfr = np.zeros(ny + 1, dtype=float)
        theta_cfr[:ny] = theta
        theta_cfr[-1] = 2.0
    else:
        yiarray_cfr = np.arange(nypf1 + 1, nypf2 + 2, dtype=float)
        theta_cfr = theta[yiarray_cfr.astype(int) - 1].copy()
        theta_cfr[-1] = 2.0

    rxy_cfr = rxy[:ixsep, nypf1:nypf2].copy()
    rxy_cfr = np.concatenate([rxy_cfr, rxy_cfr[:, 0:1]], axis=1)

    zxy_cfr = zxy[:ixsep, nypf1:nypf2].copy()
    zxy_cfr = np.concatenate([zxy_cfr, zxy_cfr[:, 0:1]], axis=1)

    zs_cfr = z_shift[:ixsep, nypf1:nypf2].copy()
    zs_patch = 0.5 * (nu[:ixsep, nypf1] + nu[:ixsep, nypf2 - 1]) * dy0 + zs_cfr[:, -1]
    zs_cfr = np.concatenate([zs_cfr, zs_patch[:, None]], axis=1)

    a1 = rxy * bpxy * btxy / hthe
    a2 = bxy**2
    a3 = sinty * a1
    jj = 4.0 * np.pi * 1.0e-7 * bpxy / hthe / (bxy**2) * jpar0

    print("Loading perturbed field information ...")
    apar_raw = _load_apar_array(
        config.aparfile,
        nx,
        ny,
        nz,
        config.apar_variable,
        config.time_index,
    )

    if divertor == 0:
        apar0, dapardx0, dapardy0, dapardz0 = get_apar_sc(
            apar_raw,
            bxy,
            psixy,
            z_shift,
            sa,
            sinty,
            dy0,
            dz,
            zperiod,
            0,
            1,
            False,
        )
    elif divertor == 1:
        apar0, dapardx0, dapardy0, dapardz0 = get_apar_sn(
            apar_raw,
            bxy,
            psixy,
            z_shift,
            sa,
            sinty,
            dy0,
            dz,
            ixsep,
            nypf1,
            nypf2,
            zperiod,
            0,
            1,
            False,
        )
    else:
        raise NotImplementedError("Double-null configuration is not implemented")

    if apar0.shape[2] == nzg:
        apar = apar0.copy()
        dapardx = dapardx0.copy()
        dapardy = dapardy0.copy()
        dapardz = dapardz0.copy()
    elif apar0.shape[2] == nz:
        apar = np.tile(apar0, (1, 1, zperiod))
        dapardx = np.tile(dapardx0, (1, 1, zperiod))
        dapardy = np.tile(dapardy0, (1, 1, zperiod))
        dapardz = np.tile(dapardz0, (1, 1, zperiod))
    else:
        raise ValueError(
            f"Unexpected Apar z-size={apar0.shape[2]} (expected {nz} or {nzg})"
        )

    b0dgy = bpxy / hthe
    inv_bxy = 1.0 / bxy[:, :, None]

    bdgx = inv_bxy * (-a1[:, :, None] * dapardy - a2[:, :, None] * dapardz) + apar * bxcvx[:, :, None]
    bdgy = inv_bxy * (a1[:, :, None] * dapardx - a3[:, :, None] * dapardz) + apar * (
        bxcvy[:, :, None] + jj[:, :, None]
    )
    bdgz = inv_bxy * (a2[:, :, None] * dapardx + a3[:, :, None] * dapardz) + apar * bxcvz[:, :, None]

    dxdy = bdgx / (b0dgy[:, :, None] + bdgy)
    dzdy = bdgz / (b0dgy[:, :, None] + bdgy)

    if config.equilibrium:
        dxdy.fill(0.0)
        dzdy.fill(0.0)

    dxdyt = dxdy[:, nypf1, :].copy()
    dxdyt = np.concatenate([dxdyt, dxdyt[:, 0:1]], axis=1)
    dzdyt = dzdy[:, nypf1, :].copy()
    dzdyt = np.concatenate([dzdyt, dzdyt[:, 0:1]], axis=1)

    dxdy_p1 = np.zeros((nx, nzg), dtype=float)
    dzdy_p1 = np.zeros((nx, nzg), dtype=float)
    for ix in range(ixsep):
        zarray_shift = np.mod(zarray[:nzg] + sa[ix], zmax)
        dxdy_p1[ix, :] = interp1_spline(zarray, dxdyt[ix, :], zarray_shift)
        dzdy_p1[ix, :] = interp1_spline(zarray, dzdyt[ix, :], zarray_shift)

    dxdyt = dxdy[:, nypf2 - 1, :].copy()
    dxdyt = np.concatenate([dxdyt, dxdyt[:, 0:1]], axis=1)
    dzdyt = dzdy[:, nypf2 - 1, :].copy()
    dzdyt = np.concatenate([dzdyt, dzdyt[:, 0:1]], axis=1)

    dxdy_m1 = np.zeros((nx, nzg), dtype=float)
    dzdy_m1 = np.zeros((nx, nzg), dtype=float)
    for ix in range(ixsep):
        zarray_rshift = np.mod(zarray[:nzg] - sa[ix], zmax)
        dxdy_m1[ix, :] = interp1_spline(zarray, dxdyt[ix, :], zarray_rshift)
        dzdy_m1[ix, :] = interp1_spline(zarray, dzdyt[ix, :], zarray_rshift)

    print("Starting field-line tracing ...")

    output_dir = Path(config.output_dir)
    if config.save_traj:
        (output_dir / "mat_traj").mkdir(parents=True, exist_ok=True)
    if config.save_pp:
        (output_dir / "mat_pp").mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    ip_fid = (output_dir / "ip_xyz.txt").open("w", encoding="utf-8")
    traj_fid = (output_dir / "traj_xyz.txt").open("w", encoding="utf-8")
    ip_fid.write("iline it ipx ipy ipz\n")
    traj_fid.write("iline it x y z\n")

    if config.lines is None:
        lines_to_trace = list(range(1, nlines + 1))
    else:
        lines_to_trace = list(config.lines)
        for iline in lines_to_trace:
            if iline < 1 or iline > nx:
                raise ValueError(f"Line index {iline} is out of bounds [1, {nx}]")

    for iline in lines_to_trace:
        xind = float(iline)
        x_start = float(psixy[iline - 1, jyomp - 1])
        yyy = jyomp
        y_start = float(jyomp)
        zzz = 1
        z_start = float(zarray[0])

        px_list: list[float] = []
        py_list: list[float] = []
        pz_list: list[float] = []
        pt_list: list[float] = []
        pp_list: list[float] = []

        traj_records: list[list[float]] = []
        iturn = 1

        if xind < float(ixsep) + 0.5:
            region = 0
            if y_start < nypf1 + 1 or y_start > nypf2:
                region = 2
        else:
            region = 1

        zind = float(interp1_linear(zarray, ziarray, z_start))

        print(f"\tline {iline} started at indeices ({xind:.6f},{y_start:.6f},{zind:.6f}),")

        if divertor == 1:
            if int(np.rint(y_start)) == ny and direction == 1:
                print(f"\tline {iline} starts on the divertor.")
                region = 14
            elif int(np.rint(y_start)) == 1 and direction == -1:
                print(f"\tline {iline} starts on the divertor.")
                region = 13

        traj_records.append([1.0, xind, y_start, zind, float(region), 0.0, z_start])

        while region < 10 and iturn < nturns:
            if iturn % 50 == 1:
                print(f"\t\t line{iline}, turn {iturn}/{nturns} ...")

            for _ in range(1, ny):
                if region == 0 and y_start > nypf1 and y_start < nypf2 + 1:
                    if config.equilibrium:
                        x_end = x_start
                        z_end = z_start
                        y_end = y_start + 1.0 if direction == 1 else y_start - 1.0
                    elif direction == 1:
                        x_end, z_end = rk4_flt1(
                            x_start,
                            y_start,
                            z_start,
                            dxdy,
                            dzdy,
                            xarray,
                            zarray,
                            region,
                            dxdy_p1,
                            dzdy_p1,
                            1,
                            nypf1,
                            nypf2,
                        )
                        y_end = y_start + 1.0
                    else:
                        x_end, z_end = rk4_flt1(
                            x_start,
                            y_start,
                            z_start,
                            dxdy,
                            dzdy,
                            xarray,
                            zarray,
                            region,
                            dxdy_m1,
                            dzdy_m1,
                            -1,
                            nypf1,
                            nypf2,
                        )
                        y_end = y_start - 1.0

                    z_raw = z_end

                    if x_end > x_max:
                        print(f"\tstarting xind={xind:.6f}, line {iline} reaches outer bndry")
                        region = 12
                    elif x_end < x_min:
                        print(f"\tstarting xind={xind:.6f}, line {iline} reaches inner bndry")
                        region = 11
                    else:
                        xind = float(interp1_linear(xarray, xiarray, x_end))
                        if xind > float(ixsep1) + 0.5:
                            region = 1
                            print(f"\tending xind={xind:.6f}, line {iline} enters the SOL")

                    if direction == 1 and int(np.rint(y_start)) == nypf2 and region == 0:
                        shiftangle = float(interp1_linear(xiarray, sa, xind))
                        z_end = z_end + shiftangle
                        y_end = float(nypf1 + 1)

                    if direction == -1 and int(np.rint(y_start)) == nypf1 + 1 and region == 0:
                        shiftangle = float(interp1_linear(xiarray, sa, xind))
                        z_end = z_end - shiftangle
                        y_end = float(nypf2)

                    if z_end < zmin or z_end > zmax:
                        z_end = float(np.mod(z_end, zmax))

                    zind = float(interp1_linear(zarray, ziarray, z_end))
                    hx = clamp_index_1based(xind, nx) - 1
                    hy = clamp_index_1based(y_end, ny) - 1
                    seg_len = float(hthe[hx, hy])

                    traj_records.append([
                        float(iturn),
                        xind,
                        y_end,
                        zind,
                        float(region),
                        seg_len,
                        float(z_raw),
                    ])

                    x_start = x_end
                    y_start = y_end
                    z_start = z_end

                elif region in (1, 2):
                    if config.equilibrium:
                        x_end = x_start
                        z_end = z_start
                        y_end = y_start + 1.0 if direction == 1 else y_start - 1.0
                    elif direction == 1:
                        x_end, z_end = rk4_flt1(
                            x_start,
                            y_start,
                            z_start,
                            dxdy,
                            dzdy,
                            xarray,
                            zarray,
                            region,
                            dxdy_p1,
                            dzdy_p1,
                            1,
                            nypf1,
                            nypf2,
                        )
                        y_end = y_start + 1.0
                    else:
                        x_end, z_end = rk4_flt1(
                            x_start,
                            y_start,
                            z_start,
                            dxdy,
                            dzdy,
                            xarray,
                            zarray,
                            region,
                            dxdy_m1,
                            dzdy_m1,
                            -1,
                            nypf1,
                            nypf2,
                        )
                        y_end = y_start - 1.0

                    z_raw = z_end

                    if direction == 1 and int(np.rint(y_start)) == nypf1 and region == 2:
                        y_end = float(nypf2 + 1)
                    elif direction == -1 and int(np.rint(y_start)) == nypf2 + 1 and region == 2:
                        y_end = float(nypf1)

                    if x_end > x_max:
                        print(f"\tstarting xind={xind:.6f}, line {iline} reaches outer bndry")
                        region = 12
                    elif x_end < x_min:
                        print(f"\tstarting xind={xind:.6f}, line {iline} reaches inner bndry")
                        region = 11
                    else:
                        xind = float(interp1_linear(xarray, xiarray, x_end))
                        if xind < float(ixsep1) + 0.5 and y_end > nypf1 and y_end < nypf2 + 1:
                            if region != 0:
                                print(f"\tending xind={xind:.6f}, line {iline} enters the CFR")
                            region = 0
                        elif xind < float(ixsep1) + 0.5 and (y_end > nypf2 - 1 or y_end < nypf1):
                            if region != 2:
                                print(f"\tending xind={xind:.6f}, line {iline} enters the PFR")
                            region = 2

                    if direction == 1 and int(np.rint(y_end)) == ny:
                        print(f"\tstarting xind={xind:.6f}, line {iline} reaches divertor")
                        region = 14
                    elif direction == -1 and int(np.rint(y_end)) == 1:
                        print(f"\tstarting xind={xind:.6f}, line {iline} reaches divertor")
                        region = 13

                    if z_end < zmin or z_end > zmax:
                        z_end = float(np.mod(z_end, zmax))

                    zind = float(interp1_linear(zarray, ziarray, z_end))
                    hx = clamp_index_1based(xind, nx) - 1
                    hy = clamp_index_1based(y_start, ny) - 1
                    seg_len = float(hthe[hx, hy])

                    traj_records.append([
                        float(iturn),
                        xind,
                        y_end,
                        zind,
                        float(region),
                        seg_len,
                        float(z_raw),
                    ])

                    x_start = x_end
                    y_start = y_end
                    z_start = z_end
                else:
                    break

            iturn += 1

        traj = np.asarray(traj_records, dtype=float).T
        itmax = traj.shape[1]

        if config.save_traj:
            suffix = "p.mat" if direction == 1 else "m.mat"
            traj_file = output_dir / "mat_traj" / f"x{iline}y{yyy}z{zzz}_v3lc-01-250{suffix}"
            parsave(str(traj_file), traj)

        fl_x3d = np.zeros(itmax, dtype=float)
        fl_y3d = np.zeros(itmax, dtype=float)
        fl_z3d = np.zeros(itmax, dtype=float)

        for istep in range(itmax):
            ystep = clamp_index_1based(traj[2, istep], ny)
            xstep = float(traj[1, istep])

            rxyvalue = float(interp1_linear(xiarray, rxy[:, ystep - 1], xstep))
            zsvalue = float(interp1_linear(xiarray, z_shift[:, ystep - 1], xstep))
            zvalue = float(interp1_linear(ziarray, zarray, traj[3, istep]))

            x3d_tmp = rxyvalue * np.cos(zsvalue)
            y3d_tmp = rxyvalue * np.sin(zsvalue)
            fl_x3d[istep] = x3d_tmp * np.cos(zvalue) - y3d_tmp * np.sin(zvalue)
            fl_y3d[istep] = x3d_tmp * np.sin(zvalue) + y3d_tmp * np.cos(zvalue)
            fl_z3d[istep] = float(interp1_linear(xiarray, zxy[:, ystep - 1], xstep))
            traj_fid.write(
                f"{iline} {istep + 1} {fl_x3d[istep]:.16g} {fl_y3d[istep]:.16g} {fl_z3d[istep]:.16g}\n"
            )

        if itmax > 1:
            crossings = _find_crossings(fl_x3d)
            id_count = 0

            if crossings:
                for root in crossings:
                    it = int(np.floor(root))
                    if it < 1 or it >= itmax:
                        continue

                    a = root - it
                    b = 1.0 - a

                    xind_tmp = b * traj[1, it - 1] + a * traj[1, it]
                    yind_tmp = b * traj[2, it - 1] + a * traj[2, it]

                    zvalue = b * traj[6, it - 1] + a * traj[6, it]
                    if abs(traj[6, it - 1] - traj[6, it]) > 1.0:
                        zvalue = b * np.mod(traj[6, it - 1], zmax) + a * np.mod(traj[6, it], zmax)

                    if (
                        int(np.rint(traj[2, it - 1])) == nypf2
                        and direction == 1
                        and xind_tmp < float(ixsep) + 0.5
                    ):
                        yind_tmp = b * traj[2, it - 1] + a * float(nypf2 + 1)
                    elif (
                        int(np.rint(traj[2, it - 1])) == nypf1 + 1
                        and direction == -1
                        and xind_tmp < float(ixsep) + 0.5
                    ):
                        yind_tmp = b * float(nypf2 + 1) + a * traj[2, it]
                        shiftangle = float(interp1_linear(xiarray, sa, xind_tmp))
                        zvalue = float(np.mod(zvalue - shiftangle, zmax))
                    elif it > 1 and (
                        int(np.rint(traj[2, it - 2])) == nypf2
                        or int(np.rint(traj[2, it - 2])) == nypf1 + 1
                    ):
                        zvalue = b * float(interp1_linear(ziarray, zarray, traj[3, it - 1])) + a * float(
                            interp1_linear(ziarray, zarray, traj[3, it])
                        )

                    if xind_tmp < float(ixsep) + 0.5:
                        rxyvalue = interp2_spline(xiarray_cfr, yiarray_cfr, rxy_cfr.T, xind_tmp, yind_tmp)
                        zxyvalue = interp2_spline(xiarray_cfr, yiarray_cfr, zxy_cfr.T, xind_tmp, yind_tmp)
                        zsvalue = interp2_spline(xiarray_cfr, yiarray_cfr, zs_cfr.T, xind_tmp, yind_tmp)
                    else:
                        rxyvalue = interp2_spline(xiarray, yiarray, rxy.T, xind_tmp, yind_tmp)
                        zxyvalue = interp2_spline(xiarray, yiarray, zxy.T, xind_tmp, yind_tmp)
                        zsvalue = interp2_spline(xiarray, yiarray, z_shift.T, xind_tmp, yind_tmp)

                    ipx3d_tmp = rxyvalue * np.cos(zsvalue)
                    ipy3d_tmp = rxyvalue * np.sin(zsvalue)
                    ipx = ipx3d_tmp * np.cos(zvalue) - ipy3d_tmp * np.sin(zvalue)
                    ipy = ipx3d_tmp * np.sin(zvalue) + ipy3d_tmp * np.cos(zvalue)
                    ipz = zxyvalue

                    if ipy > 0.0:
                        px_list.append(float(ipx))
                        py_list.append(float(ipy))
                        pz_list.append(float(ipz))
                        ip_fid.write(f"{iline} {it} {ipx:.16g} {ipy:.16g} {ipz:.16g}\n")
                        pt_list.append(float(interp1_linear(yiarray_cfr, theta_cfr, yind_tmp)))
                        pp_list.append(float(interp1_linear(xiarray, xarray, xind_tmp)))

                print(f"\t\tline {iline} has {len(px_list)}(+{id_count}) interception points.")

            if px_list:
                pxp = np.asarray(px_list, dtype=float)
                pyp = np.asarray(py_list, dtype=float)
                pzp = np.asarray(pz_list, dtype=float)
                ptp = np.asarray(pt_list, dtype=float)
                ppp = np.asarray(pp_list, dtype=float)
            else:
                pxp = np.array([0.0], dtype=float)
                pyp = np.array([0.0], dtype=float)
                pzp = np.array([0.0], dtype=float)
                ptp = np.array([0.0], dtype=float)
                ppp = np.array([0.0], dtype=float)

            traj0 = np.vstack([fl_x3d, fl_y3d, fl_z3d])
            lc = float(np.sum(traj[5, :itmax]))
            print(f"\t\tline {iline}: connection length is {lc:.6f} and ends at region={int(region)}.")

            if config.save_pp:
                suffix = "p.mat" if direction == 1 else "m.mat"
                pp_file = output_dir / "mat_pp" / f"x{iline}y{yyy}z{zzz}_v3lc-00-250{suffix}"
                parsave8(str(pp_file), pxp, pyp, pzp, ptp, ppp, traj0, lc, region)

    ip_fid.close()
    traj_fid.close()


def _parse_args() -> FieldlineTracingConfig:
    parser = argparse.ArgumentParser(description="Python port of fieldline_tracing.m")
    parser.add_argument("--case", type=int, choices=[0, 1], default=1)
    parser.add_argument("--gridfile")
    parser.add_argument("--aparfile")
    parser.add_argument("--nx", type=int)
    parser.add_argument("--ny", type=int)
    parser.add_argument("--nz", type=int)
    parser.add_argument("--zperiod", type=int)
    parser.add_argument("--divertor", type=int, choices=[0, 1])
    parser.add_argument("--direction", type=int, default=1)
    parser.add_argument("--nlines", type=int, default=256)
    parser.add_argument("--nturns", type=int, default=250)
    parser.add_argument("--np-max", type=int, default=1250)
    parser.add_argument("--save-traj", action="store_true", default=True)
    parser.add_argument("--no-save-traj", action="store_false", dest="save_traj")
    parser.add_argument("--save-pp", action="store_true", default=True)
    parser.add_argument("--no-save-pp", action="store_false", dest="save_pp")
    parser.add_argument("--output-dir", default=".")
    parser.add_argument("--apar-variable", default="apar")
    parser.add_argument("--time-index", type=int, default=-1)
    parser.add_argument("--lines", help='Comma-separated line indices, e.g. "50, 100"')
    parser.add_argument("--equilibrium", action="store_true", help="Trace equilibrium field (dxdy=dzdy=0)")

    args = parser.parse_args()
    defaults = _case_defaults(args.case)

    return FieldlineTracingConfig(
        gridfile=args.gridfile or defaults["gridfile"],
        aparfile=args.aparfile or defaults["aparfile"],
        nx=args.nx if args.nx is not None else defaults["nx"],
        ny=args.ny if args.ny is not None else defaults["ny"],
        nz=args.nz if args.nz is not None else defaults["nz"],
        zperiod=args.zperiod if args.zperiod is not None else defaults["zperiod"],
        divertor=args.divertor if args.divertor is not None else defaults["divertor"],
        direction=args.direction,
        nlines=args.nlines,
        nturns=args.nturns,
        np_max=args.np_max,
        save_traj=args.save_traj,
        save_pp=args.save_pp,
        output_dir=args.output_dir,
        apar_variable=args.apar_variable,
        time_index=args.time_index,
        lines=_parse_line_selection(args.lines),
        equilibrium=args.equilibrium,
    )


def main() -> None:
    cfg = _parse_args()
    run_fieldline_tracing(cfg)


if __name__ == "__main__":
    main()
