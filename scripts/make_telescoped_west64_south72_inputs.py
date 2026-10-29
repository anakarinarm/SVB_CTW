#!/usr/bin/env python3
"""Build CTW/SVB inputs with added flat western and southern buffers.

The domain is 576 x 684 x 100. The original 512 x 612 domain is preserved and
inserted 64 columns east and 72 rows north of the new open boundaries.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "input"
BASE_COMMON = ROOT / "ctw_reruns" / "common_input"
OUT = INPUT / "telescoped_west64_south72"

OLD_NX, OLD_NY, NX, NY, NR = 512, 612, 576, 684, 100
EXTRA_WEST = NX - OLD_NX
EXTRA_SOUTH = NY - OLD_NY
FLAT_DEPTH_M = 4400.0
MIN_WET_DEPTH_M = 2.0
NORTH_BUFFER = 40
TRANSITION = 30
TELESCOPE_FACTOR = 2.5
SOUTH_WET_OLD_COLS = 479
NORTH_WET_OLD_COLS = 73
FORCING_RECORDS = 768

CASES = {
    "realistic_plusSVB": "bathy_realistic_plusSVB_W64S72.bin",
    "realistic_minusSVB": "bathy_realistic_minusSVB_W64S72.bin",
    "smooth30_plusSVB": "bathy_smooth30_plusSVB_W64S72.bin",
    "smooth30_minusSVB": "bathy_smooth30_minusSVB_W64S72.bin",
}
BINARY_BATHY_CASES = {
    "bathy_smooth50_plusSVB_FINAL.bin": "bathy_smooth50_plusSVB_W64S72.bin",
    "bathy_smooth50_minusSVB_FINAL.bin": "bathy_smooth50_minusSVB_W64S72.bin",
    "bathy_smooth70_plusSVB_FINAL.bin": "bathy_smooth70_plusSVB_W64S72.bin",
    "bathy_smooth70_minusSVB_FINAL.bin": "bathy_smooth70_minusSVB_W64S72.bin",
}
WIND_FILES = {
    "upwelling_merid_halfgauss2hr_16day_512x612.bin": "upwelling_merid_halfgauss2hr_16day_576x684_W64S72.bin",
    "upwelling_zonal_halfgauss2hr_16day_512x612.bin": "upwelling_zonal_halfgauss2hr_16day_576x684_W64S72.bin",
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
    delx[:EXTRA_WEST] = edge_to_interior_spacing(dx0, EXTRA_WEST, TELESCOPE_FACTOR)
    x_origin = float(old_x[0] - np.sum(delx[:EXTRA_WEST]) - 0.5 * delx[EXTRA_WEST])

    dely = np.full(NY, dy0, dtype=float)
    dely[:EXTRA_SOUTH] = edge_to_interior_spacing(dy0, EXTRA_SOUTH, TELESCOPE_FACTOR)
    dely[NY - NORTH_BUFFER :] = interior_to_edge_spacing(dy0, NORTH_BUFFER, TELESCOPE_FACTOR)
    y_origin = float(old_y[0] - np.sum(dely[:EXTRA_SOUTH]) - 0.5 * dely[EXTRA_SOUTH])
    return delx, dely, x_origin, y_origin


def old_north_weight() -> np.ndarray:
    return np.ones((OLD_NY, 1), dtype=float)


def smooth_wet_rectangle(depth: np.ndarray, y0: int, y1: int, x0: int, x1: int, iterations: int = 250) -> None:
    """Diffuse wet depths inside a rectangle while holding its edge fixed."""
    wet = depth[y0:y1, x0:x1] > 0.0
    update = wet.copy()
    update[[0, -1], :] = False
    update[:, [0, -1]] = False
    sub = depth[y0:y1, x0:x1].copy()
    update[1:-1, 1:-1] &= wet[:-2, 1:-1] & wet[2:, 1:-1] & wet[1:-1, :-2] & wet[1:-1, 2:]
    for _ in range(iterations):
        new = sub.copy()
        new[1:-1, 1:-1] = 0.25 * (
            sub[:-2, 1:-1] + sub[2:, 1:-1] + sub[1:-1, :-2] + sub[1:-1, 2:]
        )
        sub[update] = new[update]
    depth[y0:y1, x0:x1] = np.where(wet, sub, depth[y0:y1, x0:x1])


def smooth_wet_line(values: np.ndarray, wet: np.ndarray, iterations: int = 500) -> np.ndarray:
    """Diffuse a one-dimensional wet profile while holding segment ends fixed."""
    out = values.copy()
    update = wet.copy()
    update[[0, -1]] = False
    update[1:-1] &= wet[:-2] & wet[2:]
    for _ in range(iterations):
        new = out.copy()
        new[1:-1] = 0.5 * (out[:-2] + out[2:])
        out[update] = new[update]
    return out


def build_bathy(depth_old: np.ndarray, north_weight: np.ndarray) -> np.ndarray:
    old_wet = depth_old > 0.0
    depth_old_buffered = np.where(old_wet, FLAT_DEPTH_M * (1.0 - north_weight) + depth_old * north_weight, 0.0)
    depth_old_buffered[(depth_old_buffered > 0.0) & (depth_old_buffered < MIN_WET_DEPTH_M)] = MIN_WET_DEPTH_M

    depth = np.zeros((NY, NX), dtype=float)
    depth[EXTRA_SOUTH:, EXTRA_WEST:] = depth_old_buffered

    # Smooth the depth jump across the added western buffer without modifying
    # the original domain.
    west_ramp = smoothstep((np.arange(EXTRA_WEST, dtype=float) + 0.5) / EXTRA_WEST)
    old_west_depth = np.where(depth_old_buffered[:, 0] > 0.0, depth_old_buffered[:, 0], 0.0)
    for ix, weight in enumerate(west_ramp):
        blended = FLAT_DEPTH_M * (1.0 - weight) + old_west_depth * weight
        depth[EXTRA_SOUTH:, ix] = np.where(old_west_depth > 0.0, blended, FLAT_DEPTH_M)

    # Smooth the depth jump across the added southern buffer. The target row is
    # the completed first active-domain row, including the western-buffer blend,
    # so the southwest and southeast corners remain continuous.
    south_ramp = smoothstep((np.arange(EXTRA_SOUTH, dtype=float) + 0.5) / EXTRA_SOUTH)
    south_wet = np.arange(NX) < EXTRA_WEST + SOUTH_WET_OLD_COLS
    south_target = depth[EXTRA_SOUTH, :]
    for jy, weight in enumerate(south_ramp):
        blended = FLAT_DEPTH_M * (1.0 - weight) + south_target * weight
        depth[jy, south_wet] = np.where(south_target[south_wet] > 0.0, blended[south_wet], FLAT_DEPTH_M)

    # Keep the northern buffer/taper band locally smooth. The finalized source
    # bathymetry already encodes the scientific +SVB/-SVB differences, but the
    # northern open-boundary buffer should not add an artificial shelf/slope
    # step. Continue the local depth from just south of the buffer over the wet
    # north-boundary band while preserving the land mask.
    north_start = NY - NORTH_BUFFER - TRANSITION
    north_target = depth[north_start - 1, :]
    north_existing_wet = depth[north_start:, :] > 0.0
    north_target_wet = north_target > 0.0
    north_profile = north_target.copy()
    north_wet_span = slice(0, EXTRA_WEST + 1)
    north_profile[north_wet_span] = smooth_wet_line(
        north_target[north_wet_span],
        north_target_wet[north_wet_span],
    )
    depth[north_start:, :] = np.where(
        north_existing_wet & north_target_wet[None, :],
        north_target[None, :],
        depth[north_start:, :],
    )
    north_ramp = smoothstep((np.arange(NY - north_start, dtype=float) + 0.5) / (NY - north_start))
    for offset, weight in enumerate(north_ramp):
        target = (1.0 - weight) * north_target + weight * north_profile
        replace = north_existing_wet[offset] & north_target_wet
        depth[north_start + offset, replace] = target[replace]

    # Remove residual creases where the one-dimensional west, south, and north
    # sponge constructions meet. These patches are outside the pristine bay and
    # hold their outer/active edges fixed while smoothing only wet interiors.
    smooth_wet_rectangle(depth, 0, EXTRA_SOUTH + 1, 0, EXTRA_WEST + 1)
    return -depth


def wind_taper() -> np.ndarray:
    west_ramp = smoothstep((np.arange(NX)[None, :] + 0.5) / EXTRA_WEST)
    south_ramp = smoothstep((np.arange(NY)[:, None] + 0.5) / EXTRA_SOUTH)
    taper = np.minimum(west_ramp, south_ramp)
    taper[EXTRA_SOUTH:, EXTRA_WEST:] = 1.0

    j_old = np.arange(OLD_NY)[:, None]
    jn0 = OLD_NY - NORTH_BUFFER - TRANSITION
    north = 1.0 - smoothstep((j_old - jn0) / TRANSITION)
    north_full = np.ones((NY, NX), dtype=float)
    north_full[EXTRA_SOUTH:, :] = north
    taper = np.minimum(taper, north_full)
    return taper


def write_big_endian(path: Path, array: np.ndarray) -> None:
    np.ascontiguousarray(array, dtype=">f8").tofile(path)


def extend_wind_record(old: np.ndarray) -> np.ndarray:
    extended = np.zeros((NY, NX), dtype=float)
    extended[EXTRA_SOUTH:, EXTRA_WEST:] = old
    extended[EXTRA_SOUTH:, :EXTRA_WEST] = old[:, :1]
    extended[:EXTRA_SOUTH, EXTRA_WEST:] = old[:1, :]
    extended[:EXTRA_SOUTH, :EXTRA_WEST] = old[0, 0]
    return extended


def write_tapered_extended_wind(src: Path, dst: Path, taper: np.ndarray) -> None:
    old_record_values = OLD_NY * OLD_NX
    values = src.stat().st_size // 8
    if values % old_record_values:
        raise ValueError(f"{src} size is not an integer number of old XY records")
    records = values // old_record_values
    if records != FORCING_RECORDS:
        raise ValueError(f"{src} has {records} records, expected {FORCING_RECORDS}")
    data = np.memmap(src, dtype=">f8", mode="r", shape=(records, OLD_NY, OLD_NX))
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    with tmp.open("wb") as out:
        for rec in range(records):
            np.ascontiguousarray(extend_wind_record(np.asarray(data[rec])) * taper, dtype=">f8").tofile(out)
    tmp.replace(dst)


def write_hydrography_and_obcs(theta_src: Path, salt_src: Path, out: Path) -> None:
    theta = np.fromfile(theta_src, dtype=">f8").reshape(NR, OLD_NY, OLD_NX)
    salt = np.fromfile(salt_src, dtype=">f8").reshape(NR, OLD_NY, OLD_NX)
    tprof = theta[:, 0, 0]
    sprof = salt[:, 0, 0]
    if np.max(np.abs(theta - tprof[:, None, None])) > 0:
        raise ValueError("Theta initial condition is not horizontally homogeneous")
    if np.max(np.abs(salt - sprof[:, None, None])) > 0:
        raise ValueError("Salt initial condition is not horizontally homogeneous")

    theta_ext = np.broadcast_to(tprof[:, None, None], (NR, NY, NX))
    salt_ext = np.broadcast_to(sprof[:, None, None], (NR, NY, NX))
    write_big_endian(out / "feb_Tref_100zlev_576x684_W64S72.bin", theta_ext)
    write_big_endian(out / "feb_Sref_100zlev_576x684_W64S72.bin", salt_ext)

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
        north_weight = old_north_weight()
        taper = wind_taper()
        write_big_endian(OUT / "delX_W64S72.bin", delx)
        write_big_endian(OUT / "delY_W64S72.bin", dely)
        write_big_endian(OUT / "wind_taper_W64S72.bin", taper)
        for name, out_name in CASES.items():
            old_depth = np.asarray(ds.variables[f"Depth_{name}"][:], dtype=float)
            write_big_endian(OUT / out_name, build_bathy(old_depth, north_weight))
        for src_name, out_name in BINARY_BATHY_CASES.items():
            old_depth = -np.fromfile(INPUT / src_name, dtype=">f8").reshape(OLD_NY, OLD_NX)
            write_big_endian(OUT / out_name, build_bathy(old_depth, north_weight))

    for src_name, dst_name in WIND_FILES.items():
        write_tapered_extended_wind(BASE_COMMON / src_name, OUT / dst_name, taper)

    write_hydrography_and_obcs(
        BASE_COMMON / "feb_Tref_100zlev_512x612.bin",
        BASE_COMMON / "feb_Sref_100zlev_512x612.bin",
        OUT,
    )

    rows = [
        ("nx", NX),
        ("ny", NY),
        ("old_nx", OLD_NX),
        ("old_ny", OLD_NY),
        ("extra_west_cols", EXTRA_WEST),
        ("extra_south_rows", EXTRA_SOUTH),
        ("extra_west_degrees", float(np.sum(delx[:EXTRA_WEST]))),
        ("extra_south_degrees", float(np.sum(dely[:EXTRA_SOUTH]))),
        ("extra_west_km_at_30N_approx", float(np.sum(delx[:EXTRA_WEST]) * 111.0 * np.cos(np.deg2rad(30.0)))),
        ("extra_south_km_approx", float(np.sum(dely[:EXTRA_SOUTH]) * 111.0)),
        ("north_buffer_cells", NORTH_BUFFER),
        ("transition_cells", TRANSITION),
        ("flat_depth_m", FLAT_DEPTH_M),
        ("south_wet_columns", EXTRA_WEST + SOUTH_WET_OLD_COLS),
        ("north_wet_columns", EXTRA_WEST + NORTH_WET_OLD_COLS),
        ("west_wet_rows", NY),
        ("xgOrigin_data_value", data_edge_origin(x_origin, True)),
        ("ygOrigin", y_origin),
        ("delx_min_deg", float(delx.min())),
        ("delx_max_deg", float(delx.max())),
        ("dely_min_deg", float(dely.min())),
        ("dely_max_deg", float(dely.max())),
    ]
    with (OUT / "telescoped_west64_south72_summary.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["parameter", "value"])
        writer.writerows(rows)
    print(f"wrote west64/south72 telescoped inputs to {OUT}")


if __name__ == "__main__":
    main()
