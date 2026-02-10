from __future__ import annotations

from pathlib import Path

import scipy.io as sio


def parsave(fname: str, v) -> None:
    """Python equivalent of MATLAB parsave(fname, v)."""
    path = Path(fname)
    path.parent.mkdir(parents=True, exist_ok=True)
    sio.savemat(path, {"v": v})


def parsave8(fname: str, v1, v2, v3, v4, v5, v6, v7, v8) -> None:
    """Python equivalent of MATLAB parsave8(...)."""
    path = Path(fname)
    path.parent.mkdir(parents=True, exist_ok=True)
    sio.savemat(
        path,
        {
            "v1": v1,
            "v2": v2,
            "v3": v3,
            "v4": v4,
            "v5": v5,
            "v6": v6,
            "v7": v7,
            "v8": v8,
        },
    )
