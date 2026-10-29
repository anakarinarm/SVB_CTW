#!/usr/bin/env python3
"""Pack MITgcm-written grid files for the four production runs into NetCDF."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr
from MITgcmutils import mds


RUN_ROOT = Path("/Volumes/LaCie/CTW_SVB_reruns_telescoped_west64_south72")
OUT = Path(__file__).resolve().parents[1] / "ctw_reruns_telescoped_west64_south72" / "production_grid_bathymetries.nc"

RUNS = ["Rplus", "Rminus", "Splus", "Sminus"]
STATIC_2D = ["XC", "YC", "XG", "YG", "RAC", "RAS", "RAW", "RAZ", "DXC", "DXF", "DXG", "DXV", "DYC", "DYF", "DYG", "DYU"]
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


def rd(run: str, name: str) -> np.ndarray:
    return np.asarray(mds.rdmds(str(RUN_ROOT / run / name)))


def squeeze_1d(a: np.ndarray) -> np.ndarray:
    return np.asarray(a).squeeze()


def validate_static(reference: dict[str, np.ndarray]) -> None:
    for run in RUNS[1:]:
        for name in STATIC_2D:
            arr = rd(run, name)
            if not np.array_equal(reference[name], arr):
                raise RuntimeError(f"{name} differs between {RUNS[0]} and {run}")
        for name in STATIC_1D:
            arr = squeeze_1d(rd(run, name))
            if not np.array_equal(reference[name], arr):
                raise RuntimeError(f"{name} differs between {RUNS[0]} and {run}")


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)

    static: dict[str, np.ndarray] = {}
    for name in STATIC_2D:
        static[name] = rd(RUNS[0], name)
    for name in STATIC_1D:
        static[name] = squeeze_1d(rd(RUNS[0], name))
    validate_static(static)

    ny, nx = static["XC"].shape
    nz = static["RC"].size
    nzp1 = static["RF"].size

    data_vars: dict[str, tuple[tuple[str, ...], np.ndarray, dict[str, str]]] = {}
    for name in STATIC_2D:
        data_vars[name] = (("Y", "X"), static[name].astype("float64"), {"source": f"{RUNS[0]}/{name}.data"})
    for name, dim in STATIC_1D.items():
        data_vars[name] = ((dim,), static[name].astype("float64"), {"source": f"{RUNS[0]}/{name}.data"})

    for name in CASE_2D:
        arr = np.stack([rd(run, name) for run in RUNS]).astype("float32")
        data_vars[name] = (("run", "Y", "X"), arr, {"source": f"MITgcm {name}.data from each run"})
    for name in CASE_3D:
        arr = np.stack([rd(run, name) for run in RUNS]).astype("float32")
        data_vars[name] = (("run", "Z", "Y", "X"), arr, {"source": f"MITgcm {name}.data from each run"})

    ds = xr.Dataset(
        data_vars=data_vars,
        coords={
            "run": RUNS,
            "Y": np.arange(ny, dtype=np.int32),
            "X": np.arange(nx, dtype=np.int32),
            "Z": np.arange(nz, dtype=np.int32),
            "Zp1": np.arange(nzp1, dtype=np.int32),
        },
        attrs={
            "title": "Production MITgcm grid and bathymetry for CTW/SVB west64/south72 runs",
            "source": str(RUN_ROOT),
            "history": "Created from MITgcm-written MDS grid output files, not from input bathymetry binaries.",
            "runs": ", ".join(RUNS),
            "grid_size": f"{nx} x {ny} x {nz}",
        },
    )

    ds["Depth"].attrs.update(
        long_name="MITgcm model bathymetry",
        units="m",
        positive="down",
        note="Read from each run's MITgcm-written Depth.data grid output.",
    )
    for name in CASE_3D:
        ds[name].attrs.update(long_name=name, note=f"Read from each run's MITgcm-written {name}.data grid output.")
    ds["XC"].attrs.update(long_name="cell-center longitude", units="degrees_east")
    ds["YC"].attrs.update(long_name="cell-center latitude", units="degrees_north")
    ds["XG"].attrs.update(long_name="cell-corner longitude", units="degrees_east")
    ds["YG"].attrs.update(long_name="cell-corner latitude", units="degrees_north")
    ds["RC"].attrs.update(long_name="vertical cell-center coordinate", units="m")
    ds["RF"].attrs.update(long_name="vertical cell-interface coordinate", units="m")

    encoding = {
        name: {"zlib": True, "complevel": 4, "shuffle": True}
        for name in ds.data_vars
        if ds[name].ndim >= 2
    }
    ds.to_netcdf(OUT, encoding=encoding)
    print(f"wrote {OUT}")
    print(ds)


if __name__ == "__main__":
    main()
