#!/usr/bin/env python3
"""Pack MITgcm-written west64/south72 grid files into one NetCDF.

This reads grid output produced by MITgcm itself from each run directory. It is
intended for bathymetry/partial-cell provenance checks, so it does not read the
input bathymetry binaries directly.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from MITgcmutils import mds
from netCDF4 import Dataset


DEFAULT_RUN_ROOT = Path("/Volumes/LaCie/CTW_SVB_reruns_telescoped_west64_south72")
DEFAULT_OUT = (
    Path(__file__).resolve().parents[1]
    / "ctw_reruns_telescoped_west64_south72"
    / "CTW_final_bathymetries_50_70.nc"
)
DEFAULT_RUNS = [
    "S50plus",
    "S50minus",
    "S70plus",
    "S70minus",
]
STATIC_2D = [
    "XC",
    "YC",
    "XG",
    "YG",
    "RAC",
    "RAS",
    "RAW",
    "RAZ",
    "DXC",
    "DXF",
    "DXG",
    "DXV",
    "DYC",
    "DYF",
    "DYG",
    "DYU",
]
STATIC_1D = {
    "RC": "Z",
    "RF": "Zp1",
    "DRC": "Zp1",
    "DRF": "Z",
    "PHrefC": "Z",
    "PHrefF": "Zp1",
    "RhoRef": "Z",
}
CASE_2D = ["Depth"]
CASE_3D = ["hFacC", "hFacS", "hFacW"]


def rd(run_root: Path, run: str, name: str) -> np.ndarray:
    return np.asarray(mds.rdmds(str(run_root / run / name)))


def squeeze_1d(a: np.ndarray) -> np.ndarray:
    return np.asarray(a).squeeze()


def validate_static(run_root: Path, runs: list[str], reference: dict[str, np.ndarray]) -> None:
    for run in runs[1:]:
        for name in STATIC_2D:
            arr = rd(run_root, run, name)
            if not np.array_equal(reference[name], arr):
                raise RuntimeError(f"{name} differs between {runs[0]} and {run}")
        for name in STATIC_1D:
            arr = squeeze_1d(rd(run_root, run, name))
            if not np.array_equal(reference[name], arr):
                raise RuntimeError(f"{name} differs between {runs[0]} and {run}")


def create_var(ds: Dataset, name: str, dtype: str, dims: tuple[str, ...], **attrs):
    kwargs = {}
    if len(dims) >= 2:
        kwargs.update(zlib=True, complevel=4, shuffle=True)
    var = ds.createVariable(name, dtype, dims, **kwargs)
    for key, value in attrs.items():
        setattr(var, key, value)
    return var


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--runs", nargs="+", default=DEFAULT_RUNS)
    args = parser.parse_args()

    run_root = args.run_root
    runs = list(args.runs)
    out = args.output
    out.parent.mkdir(parents=True, exist_ok=True)

    static: dict[str, np.ndarray] = {}
    for name in STATIC_2D:
        static[name] = rd(run_root, runs[0], name)
    for name in STATIC_1D:
        static[name] = squeeze_1d(rd(run_root, runs[0], name))
    validate_static(run_root, runs, static)

    ny, nx = static["XC"].shape
    nz = static["RC"].size
    nzp1 = static["RF"].size

    with Dataset(out, "w") as ds:
        ds.createDimension("run", len(runs))
        ds.createDimension("Y", ny)
        ds.createDimension("X", nx)
        ds.createDimension("Z", nz)
        ds.createDimension("Zp1", nzp1)

        run_var = ds.createVariable("run", str, ("run",))
        run_var[:] = np.asarray(runs, dtype=object)
        create_var(ds, "Y", "i4", ("Y",), long_name="Y index")[:] = np.arange(ny, dtype=np.int32)
        create_var(ds, "X", "i4", ("X",), long_name="X index")[:] = np.arange(nx, dtype=np.int32)
        create_var(ds, "Z", "i4", ("Z",), long_name="vertical cell-center index")[:] = np.arange(nz, dtype=np.int32)
        create_var(ds, "Zp1", "i4", ("Zp1",), long_name="vertical cell-interface index")[:] = np.arange(nzp1, dtype=np.int32)

        for name in STATIC_2D:
            var = create_var(ds, name, "f8", ("Y", "X"), source=f"{runs[0]}/{name}.data")
            var[:, :] = static[name].astype("float64")
        for name, dim in STATIC_1D.items():
            var = create_var(ds, name, "f8", (dim,), source=f"{runs[0]}/{name}.data")
            var[:] = static[name].astype("float64")

        for name in CASE_2D:
            var = create_var(
                ds,
                name,
                "f4",
                ("run", "Y", "X"),
                source=f"MITgcm-written {name}.data from each run",
                units="m",
                positive="down",
            )
            for irun, run in enumerate(runs):
                var[irun, :, :] = rd(run_root, run, name).astype("float32")

        for name in CASE_3D:
            var = create_var(
                ds,
                name,
                "f4",
                ("run", "Z", "Y", "X"),
                source=f"MITgcm-written {name}.data from each run",
            )
            for irun, run in enumerate(runs):
                var[irun, :, :, :] = rd(run_root, run, name).astype("float32")

        ds.title = "MITgcm grid and bathymetry for CTW/SVB west64/south72 runs"
        ds.source = str(run_root)
        ds.history = "Created from MITgcm-written MDS grid output files, not from input bathymetry binaries."
        ds.runs = ", ".join(runs)
        ds.grid_size = f"{nx} x {ny} x {nz}"

    print(f"wrote {out}")
    print(f"runs: {', '.join(runs)}")
    print(f"grid_size: {nx} x {ny} x {nz}")


if __name__ == "__main__":
    main()
