#!/usr/bin/env python3
"""Build telescoped CTW/SVB inputs with an added flat southern buffer.

The domain is 512 x 648 x 100: the original 512 x 612 domain is preserved and
36 new rows are prepended at the southern boundary. The added rows are a flat
abyssal plane over the same open-water columns as the old southern OBC.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "input"
BASE_COMMON = ROOT / "ctw_reruns" / "common_input"
OUT = INPUT / "telescoped_south_extended"

NX, OLD_NY, NY, NR = 512, 612, 648, 100
EXTRA_SOUTH = NY - OLD_NY
FLAT_DEPTH_M = 4000.0
MIN_WET_DEPTH_M = 2.0
WEST_BUFFER = 48
NORTH_BUFFER = 40
TRANSITION = 30
WEST_MAX_FACTOR = 4.0
Y_MAX_FACTOR = 3.0
SOUTH_WET_COLS = 479
FORCING_RECORDS = 768

CASES = {
    "realistic_plusSVB": "bathy_realistic_plusSVB_SOUTHEXT.bin",
    "realistic_minusSVB": "bathy_realistic_minusSVB_SOUTHEXT.bin",
    "smooth30_plusSVB": "bathy_smooth30_plusSVB_SOUTHEXT.bin",
    "smooth30_minusSVB": "bathy_smooth30_minusSVB_SOUTHEXT.bin",
}
WIND_FILES = {
    "upwelling_merid_halfgauss2hr_16day_512x612.bin": "upwelling_merid_halfgauss2hr_16day_512x648_SOUTHEXT.bin",
    "upwelling_zonal_halfgauss2hr_16day_512x612.bin": "upwelling_zonal_halfgauss2hr_16day_512x648_SOUTHEXT.bin",
}


def smoothstep(x: np.ndarray) -> np.ndarray:
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def edge_to_interior_spacing(base: float, count: int, max_factor: float) -> np.ndarray:
    if count <= 1:
        return np.array([base], dtype=float)
    ratio = max_factor ** (1.0 / (count - 1))
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
    dely[NY - NORTH_BUFFER :] = interior_to_edge_spacing(dy0, NORTH_BUFFER, Y_MAX_FACTOR)
    y_origin = float(old_y[0] - np.sum(dely[:EXTRA_SOUTH]) - 0.5 * dely[EXTRA_SOUTH])
    return delx, dely, x_origin, y_origin


def bathy_weight_old_domain() -> np.ndarray:
    i = np.arange(NX)[None, :]
    j = np.arange(OLD_NY)[:, None]
    wx = smoothstep((i - WEST_BUFFER) / TRANSITION)
    jn0 = OLD_NY - NORTH_BUFFER - TRANSITION
    wn = 1.0 - smoothstep((j - jn0) / TRANSITION)
    return np.minimum(wx, wn)


def wind_taper() -> np.ndarray:
    taper = np.ones((NY, NX), dtype=float)
    i = np.arange(NX)[None, :]
    j_old = np.arange(OLD_NY)[:, None]
    wx = smoothstep((i - WEST_BUFFER) / TRANSITION)
    jn0 = OLD_NY - NORTH_BUFFER - TRANSITION
    wn = 1.0 - smoothstep((j_old - jn0) / TRANSITION)
    taper[EXTRA_SOUTH:, :] = np.minimum(wx, wn)
    taper[:EXTRA_SOUTH, :] = smoothstep((np.arange(EXTRA_SOUTH)[:, None] + 0.5) / EXTRA_SOUTH)
    return taper


def build_bathy(depth_old: np.ndarray, weight_old: np.ndarray) -> np.ndarray:
    old_wet = depth_old > 0.0
    depth_old_buffered = np.where(old_wet, FLAT_DEPTH_M * (1.0 - weight_old) + depth_old * weight_old, 0.0)
    depth_old_buffered[(depth_old_buffered > 0.0) & (depth_old_buffered < MIN_WET_DEPTH_M)] = MIN_WET_DEPTH_M

    depth = np.zeros((NY, NX), dtype=float)
    depth[:EXTRA_SOUTH, :SOUTH_WET_COLS] = FLAT_DEPTH_M
    depth[EXTRA_SOUTH:, :] = depth_old_buffered
    return -depth


def write_big_endian(path: Path, array: np.ndarray) -> None:
    np.ascontiguousarray(array, dtype=">f8").tofile(path)


def write_tapered_extended_wind(src: Path, dst: Path, taper: np.ndarray) -> None:
    old_record_values = OLD_NY * NX
    values = src.stat().st_size // 8
    if values % old_record_values:
        raise ValueError(f"{src} size is not an integer number of old XY records")
    records = values // old_record_values
    if records != FORCING_RECORDS:
        raise ValueError(f"{src} has {records} records, expected {FORCING_RECORDS}")
    data = np.memmap(src, dtype=">f8", mode="r", shape=(records, OLD_NY, NX))
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    with tmp.open("wb") as out:
        extended = np.zeros((NY, NX), dtype=float)
        for rec in range(records):
            extended[:EXTRA_SOUTH, :] = 0.0
            extended[EXTRA_SOUTH:, :] = np.asarray(data[rec])
            np.ascontiguousarray(extended * taper, dtype=">f8").tofile(out)
    tmp.replace(dst)


def write_hydrography_and_obcs(theta_src: Path, salt_src: Path, out: Path) -> None:
    theta = np.fromfile(theta_src, dtype=">f8").reshape(NR, OLD_NY, NX)
    salt = np.fromfile(salt_src, dtype=">f8").reshape(NR, OLD_NY, NX)
    tprof = theta[:, 0, 0]
    sprof = salt[:, 0, 0]
    if np.max(np.abs(theta - tprof[:, None, None])) > 0:
        raise ValueError("Theta initial condition is not horizontally homogeneous")
    if np.max(np.abs(salt - sprof[:, None, None])) > 0:
        raise ValueError("Salt initial condition is not horizontally homogeneous")

    theta_ext = np.broadcast_to(tprof[:, None, None], (NR, NY, NX))
    salt_ext = np.broadcast_to(sprof[:, None, None], (NR, NY, NX))
    write_big_endian(out / "feb_Tref_100zlev_512x648_SOUTHEXT.bin", theta_ext)
    write_big_endian(out / "feb_Sref_100zlev_512x648_SOUTHEXT.bin", salt_ext)

    obcs = out / "obcs_ts_from_initial_hydrography_16day"
    obcs.mkdir(exist_ok=True)
    fields = {
        "OBWt.bin": np.repeat(tprof[:, None], NY, axis=1),
        "OBWs.bin": np.repeat(sprof[:, None], NY, axis=1),
        "OBSt.bin": np.repeat(tprof[:, None], NX, axis=1),
        "OBSs.bin": np.repeat(sprof[:, None], NX, axis=1),
        "OBNt.bin": np.repeat(tprof[:, None], NX, axis=1),
        "OBNs.bin": np.repeat(sprof[:, None], NX, axis=1),
    }
    for name, arr in fields.items():
        tmp = (obcs / name).with_suffix(".bin.tmp")
        with tmp.open("wb") as f:
            record = np.ascontiguousarray(arr, dtype=">f8")
            for _ in range(FORCING_RECORDS):
                record.tofile(f)
        tmp.replace(obcs / name)


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
        weight = bathy_weight_old_domain()
        taper = wind_taper()
        write_big_endian(OUT / "delX_SOUTHEXT.bin", delx)
        write_big_endian(OUT / "delY_SOUTHEXT.bin", dely)
        write_big_endian(OUT / "wind_taper_SOUTHEXT.bin", taper)
        for name, out_name in CASES.items():
            old_depth = np.asarray(ds.variables[f"Depth_{name}"][:], dtype=float)
            write_big_endian(OUT / out_name, build_bathy(old_depth, weight))

    for src_name, dst_name in WIND_FILES.items():
        write_tapered_extended_wind(BASE_COMMON / src_name, OUT / dst_name, taper)

    write_hydrography_and_obcs(
        BASE_COMMON / "feb_Tref_100zlev_512x612.bin",
        BASE_COMMON / "feb_Sref_100zlev_512x612.bin",
        OUT,
    )

    dy0 = float(np.median(dely[EXTRA_SOUTH : NY - NORTH_BUFFER]))
    rows = [
        ("nx", NX),
        ("ny", NY),
        ("old_ny", OLD_NY),
        ("extra_south_rows", EXTRA_SOUTH),
        ("extra_south_degrees", float(np.sum(dely[:EXTRA_SOUTH]))),
        ("extra_south_km_approx", float(np.sum(dely[:EXTRA_SOUTH]) * 111.0)),
        ("west_buffer_cells", WEST_BUFFER),
        ("north_buffer_cells", NORTH_BUFFER),
        ("transition_cells", TRANSITION),
        ("flat_depth_m", FLAT_DEPTH_M),
        ("south_wet_columns", SOUTH_WET_COLS),
        ("xgOrigin_data_value", data_edge_origin(x_origin, True)),
        ("ygOrigin", y_origin),
        ("native_dy_deg", dy0),
        ("dely_min_deg", float(dely.min())),
        ("dely_max_deg", float(dely.max())),
    ]
    with (OUT / "telescoped_south_extended_summary.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["parameter", "value"])
        writer.writerows(rows)
    print(f"wrote south-extended telescoped inputs to {OUT}")


if __name__ == "__main__":
    main()
