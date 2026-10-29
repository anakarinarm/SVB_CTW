#!/usr/bin/env python3
"""Build telescoped flat-abyss buffer inputs for CTW/SVB reruns.

The grid size remains 512 x 612. Old open-boundary strips are replaced by
wide, boring, flat abyssal buffers on the west, south, and north sides.
The interior grid cells outside those strips retain their original coordinates.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "input"
OUT = INPUT / "telescoped_buffer"
NX, NY = 512, 612
FLAT_DEPTH_M = 4000.0
MIN_WET_DEPTH_M = 2.0
WEST_BUFFER = 48
SOUTH_BUFFER = 40
NORTH_BUFFER = 40
TRANSITION = 30
WEST_MAX_FACTOR = 4.0
Y_MAX_FACTOR = 3.0
FORCING_RECORDS = 768

CASES = {
    "realistic_plusSVB": "bathy_realistic_plusSVB_TELEBUF.bin",
    "realistic_minusSVB": "bathy_realistic_minusSVB_TELEBUF.bin",
    "smooth30_plusSVB": "bathy_smooth30_plusSVB_TELEBUF.bin",
    "smooth30_minusSVB": "bathy_smooth30_minusSVB_TELEBUF.bin",
}
WIND_FILES = {
    "upwelling_merid_halfgauss2hr_16day_512x612.bin": "upwelling_merid_halfgauss2hr_16day_512x612_TELEBUF.bin",
    "upwelling_zonal_halfgauss2hr_16day_512x612.bin": "upwelling_zonal_halfgauss2hr_16day_512x612_TELEBUF.bin",
}


def smoothstep(x: np.ndarray) -> np.ndarray:
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def edge_to_interior_spacing(base: float, count: int, max_factor: float) -> np.ndarray:
    if count <= 1:
        return np.array([base], dtype=float)
    ratio = max_factor ** (1.0 / (count - 1))
    # Coarsest cell at the open boundary, native spacing at the science edge.
    return base * max_factor / ratio ** np.arange(count, dtype=float)


def interior_to_edge_spacing(base: float, count: int, max_factor: float) -> np.ndarray:
    return edge_to_interior_spacing(base, count, max_factor)[::-1]


def build_spacing(old_x: np.ndarray, old_y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    dx0 = float(np.median(np.diff(old_x)))
    dy0 = float(np.median(np.diff(old_y)))

    delx = np.full(NX, dx0, dtype=float)
    delx[:WEST_BUFFER] = edge_to_interior_spacing(dx0, WEST_BUFFER, WEST_MAX_FACTOR)
    x_origin = float(old_x[WEST_BUFFER] - np.sum(delx[:WEST_BUFFER]) - 0.5 * delx[WEST_BUFFER])

    dely = np.full(NY, dy0, dtype=float)
    dely[:SOUTH_BUFFER] = edge_to_interior_spacing(dy0, SOUTH_BUFFER, Y_MAX_FACTOR)
    dely[NY - NORTH_BUFFER :] = interior_to_edge_spacing(dy0, NORTH_BUFFER, Y_MAX_FACTOR)
    y_origin = float(old_y[SOUTH_BUFFER] - np.sum(dely[:SOUTH_BUFFER]) - 0.5 * dely[SOUTH_BUFFER])
    return delx, dely, x_origin, y_origin


def buffer_weight() -> np.ndarray:
    i = np.arange(NX)[None, :]
    j = np.arange(NY)[:, None]
    wx = smoothstep((i - WEST_BUFFER) / TRANSITION)
    ws = smoothstep((j - SOUTH_BUFFER) / TRANSITION)
    jn0 = NY - NORTH_BUFFER - TRANSITION
    wn = 1.0 - smoothstep((j - jn0) / TRANSITION)
    return np.minimum(np.minimum(wx, ws), wn)


def build_bathy(depth_old: np.ndarray, weight: np.ndarray) -> np.ndarray:
    old_wet = depth_old > 0.0
    depth = np.where(old_wet, FLAT_DEPTH_M * (1.0 - weight) + depth_old * weight, 0.0)
    depth[(depth > 0.0) & (depth < MIN_WET_DEPTH_M)] = MIN_WET_DEPTH_M
    return -depth


def write_big_endian(path: Path, array: np.ndarray) -> None:
    np.ascontiguousarray(array, dtype=">f8").tofile(path)


def write_tapered_wind(src: Path, dst: Path, taper: np.ndarray) -> None:
    record_values = NY * NX
    values = src.stat().st_size // 8
    if values % record_values:
        raise ValueError(f"{src} size is not an integer number of XY records")
    records = values // record_values
    if records != FORCING_RECORDS:
        raise ValueError(f"{src} has {records} records, expected {FORCING_RECORDS}")
    data = np.memmap(src, dtype=">f8", mode="r", shape=(records, NY, NX))
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    with tmp.open("wb") as out:
        for rec in range(records):
            np.ascontiguousarray(np.asarray(data[rec]) * taper, dtype=">f8").tofile(out)
    tmp.replace(dst)


def data_edge_origin(origin: float, is_lon: bool) -> float:
    if is_lon and origin < 0:
        return origin + 360.0
    return origin


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with Dataset(INPUT / "CTW_final_bathymetries.nc") as ds:
        old_x = np.asarray(ds.variables["XC"][:], dtype=float)
        old_y = np.asarray(ds.variables["YC"][:], dtype=float)
        delx, dely, x_origin, y_origin = build_spacing(old_x, old_y)
        weight = buffer_weight()
        write_big_endian(OUT / "delX_TELEBUF.bin", delx)
        write_big_endian(OUT / "delY_TELEBUF.bin", dely)
        write_big_endian(OUT / "buffer_taper_TELEBUF.bin", weight)
        for name, out_name in CASES.items():
            old_depth = np.asarray(ds.variables[f"Depth_{name}"][:], dtype=float)
            bathy = build_bathy(old_depth, weight)
            write_big_endian(OUT / out_name, bathy)

    common = ROOT / "ctw_reruns" / "common_input"
    for src_name, dst_name in WIND_FILES.items():
        write_tapered_wind(common / src_name, OUT / dst_name, weight)

    rows = [
        ("nx", NX),
        ("ny", NY),
        ("west_buffer_cells", WEST_BUFFER),
        ("south_buffer_cells", SOUTH_BUFFER),
        ("north_buffer_cells", NORTH_BUFFER),
        ("transition_cells", TRANSITION),
        ("flat_depth_m", FLAT_DEPTH_M),
        ("minimum_wet_depth_m", MIN_WET_DEPTH_M),
        ("xgOrigin_data_value", data_edge_origin(x_origin, True)),
        ("ygOrigin", y_origin),
        ("delx_min_deg", float(delx.min())),
        ("delx_max_deg", float(delx.max())),
        ("dely_min_deg", float(dely.min())),
        ("dely_max_deg", float(dely.max())),
        ("west_extension_km_approx", float((np.sum(delx[:WEST_BUFFER]) - WEST_BUFFER * np.median(np.diff(old_x))) * 111.0)),
        ("south_extension_km_approx", float((np.sum(dely[:SOUTH_BUFFER]) - SOUTH_BUFFER * np.median(np.diff(old_y))) * 111.0)),
        ("north_extension_km_approx", float((np.sum(dely[NY - NORTH_BUFFER :]) - NORTH_BUFFER * np.median(np.diff(old_y))) * 111.0)),
    ]
    with (OUT / "telescoped_buffer_summary.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["parameter", "value"])
        writer.writerows(rows)
    print(f"wrote telescoped buffer inputs to {OUT}")


if __name__ == "__main__":
    main()
