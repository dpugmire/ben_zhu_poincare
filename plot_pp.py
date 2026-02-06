#!/usr/bin/env python3
import argparse
import pathlib
from typing import Any

import numpy as np
import matplotlib.pyplot as plt

## plot v2, v3 shows punctures.


def loadOctaveAsciiMulti(path: str) -> dict[str, Any]:
    """
    Load an Octave ASCII save file containing blocks like:

      # name: v2
      # type: matrix
      # rows: ...
      # columns: ...
      <data...>

    and/or:

      # name: v7
      # type: scalar
      <value>

    Returns dict: name -> numpy array (matrix) or float (scalar)
    """
    out: dict[str, Any] = {}

    p = pathlib.Path(path)
    with p.open("r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    n = len(lines)
    i = 0

    def stripCommentTail(s: str) -> str:
        # Some of your lines have trailing '%' comment tails.
        if "%" in s:
            s = s.split("%", 1)[0]
        return s.strip()

    def advanceTo(prefix: str) -> None:
        nonlocal i
        while i < n and not lines[i].lstrip().startswith(prefix):
            i += 1

    while i < n:
        line = lines[i].strip()
        if not line.startswith("# name:"):
            i += 1
            continue

        name = line.split(":", 1)[1].strip()
        i += 1

        advanceTo("# type:")
        if i >= n:
            break
        vtype = lines[i].split(":", 1)[1].strip()
        i += 1

        if vtype == "matrix":
            advanceTo("# rows:")
            if i >= n:
                raise ValueError(f"Variable {name!r}: missing rows")
            rows = int(lines[i].split(":", 1)[1].strip())
            i += 1

            advanceTo("# columns:")
            if i >= n:
                raise ValueError(f"Variable {name!r}: missing columns")
            cols = int(lines[i].split(":", 1)[1].strip())
            i += 1

            needed = rows * cols
            data: list[float] = []

            while i < n and len(data) < needed:
                raw = stripCommentTail(lines[i])
                if raw and not raw.startswith("#"):
                    parts = raw.split()
                    for p in parts:
                        data.append(float(p))
                i += 1

            if len(data) < needed:
                raise ValueError(
                    f"Variable {name!r}: expected {needed} values, got {len(data)}"
                )

            arr = np.array(data[:needed], dtype=np.float64).reshape((rows, cols))
            out[name] = arr

        elif vtype == "scalar":
            value = None
            while i < n and value is None:
                raw = stripCommentTail(lines[i])
                if raw and not raw.startswith("#"):
                    value = float(raw.split()[0])
                i += 1

            if value is None:
                raise ValueError(f"Variable {name!r}: scalar value not found")

            out[name] = value

        else:
            raise ValueError(f"Unsupported type {vtype!r} for variable {name!r}")

    return out


def as1d(value: Any) -> np.ndarray:
    """
    Convert a loaded value (scalar or matrix) to a 1D numpy array suitable for plotting.
    - scalar -> shape (1,)
    - matrix -> ravel()
    """
    if isinstance(value, (float, int, np.floating, np.integer)):
        return np.array([float(value)], dtype=np.float64)

    arr = np.asarray(value, dtype=np.float64)
    return arr.ravel()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot variables from a multi-variable Octave ASCII save file."
    )
    parser.add_argument("--file", required=True, help="Path to Octave ASCII file")
    parser.add_argument("--x", required=True, help="Variable name for x (e.g., v2)")
    parser.add_argument("--y", required=True, help="Variable name for y (e.g., v3)")
    parser.add_argument(
        "--list",
        action="store_true",
        help="List variables in the file and exit",
    )
    parser.add_argument("--marker-size", type=float, default=2.0)
    parser.add_argument("--title", default=None)

    args = parser.parse_args()

    d = loadOctaveAsciiMulti(args.file)

    if args.list:
        keys = sorted(d.keys())
        print("Variables found:")
        for k in keys:
            v = d[k]
            if isinstance(v, np.ndarray):
                print(f"  {k}: matrix {v.shape}")
            else:
                print(f"  {k}: scalar {v}")
        return

    if args.x not in d:
        raise KeyError(f"x variable {args.x!r} not found. Use --list to see keys.")
    if args.y not in d:
        raise KeyError(f"y variable {args.y!r} not found. Use --list to see keys.")

    x = as1d(d[args.x])
    y = as1d(d[args.y])

    if x.size != y.size and x.size != 1 and y.size != 1:
        raise ValueError(f"x and y sizes differ: x={x.size}, y={y.size}")

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, ".", markersize=args.marker_size)
    plt.xlabel(args.x)
    plt.ylabel(args.y)
    plt.title(args.title if args.title is not None else f"{args.y} vs {args.x}")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
