#!/usr/bin/env python3
"""Generate steady T/S OBC files from the CTW/SVB initial hydrography.

The production setup uses these files to prescribe tracer values at the open
boundaries while leaving velocity/free-surface files unspecified.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
NX, NY, NR = 512, 612, 100


def write_boundary(field: np.ndarray, outdir: Path, prefix: str) -> None:
    profile = field[:, 0, 0]
    max_dev = float(np.max(np.abs(field - profile[:, None, None])))
    if max_dev > 0:
        raise ValueError(f"{prefix} field is not horizontally homogeneous: max deviation {max_dev}")
    west = np.repeat(profile[:, None], NY, axis=1)
    south = np.repeat(profile[:, None], NX, axis=1)
    north = np.repeat(profile[:, None], NX, axis=1)
    west.astype(">f8").tofile(outdir / f"OBW{prefix}.bin")
    south.astype(">f8").tofile(outdir / f"OBS{prefix}.bin")
    north.astype(">f8").tofile(outdir / f"OBN{prefix}.bin")


def main() -> None:
    common = ROOT / "input" / "common_input_512x612"
    outdir = common / "obcs_ts_from_initial_hydrography"
    outdir.mkdir(exist_ok=True)
    theta = np.fromfile(common / "feb_Tref_100zlev_512x612.bin", dtype=">f8").reshape(NR, NY, NX)
    salt = np.fromfile(common / "feb_Sref_100zlev_512x612.bin", dtype=">f8").reshape(NR, NY, NX)
    write_boundary(theta, outdir, "t")
    write_boundary(salt, outdir, "s")
    print(f"wrote steady T/S boundary files to {outdir}")


if __name__ == "__main__":
    main()
