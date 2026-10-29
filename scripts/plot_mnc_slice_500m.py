#!/usr/bin/env python3
"""Assemble CTW MNC tiles and plot U, V, and RHOAnoma near 500 m."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


DEL_R = np.array(
    [
        2, 2, 2, 2, 2, 2.2, 2.4, 2.6, 2.8,
        3.0, 3.2, 3.5, 3.8, 4.2, 4.5, 4.8, 5, 5, 5,
        5, 5, 5, 5.5, 5.5, 5.5, 5.5, 5.5, 6, 6.5,
        7, 7.5, 8, 8.5, 9.2, 10, 10, 10, 10, 10, 10,
        10, 11, 12, 13, 14, 15, 16, 17, 18, 20,
        22, 23, 25, 25, 25, 25, 25, 25, 25, 25,
        27, 28, 30, 30, 32, 34, 36, 38, 40, 43,
        46, 50, 55, 59, 65, 70, 75, 80, 85, 92,
        100, 100, 100, 100, 100, 100, 105, 115, 125, 132,
        145, 160, 170, 180, 190, 200, 210, 230, 245, 269.8,
    ],
    dtype=float,
)
Z_C = np.cumsum(DEL_R) - 0.5 * DEL_R


def choose_index(diag_levels: np.ndarray, target_depth: float) -> tuple[int, int, float]:
    model_levels = diag_levels.astype(int)
    depths = Z_C[model_levels - 1]
    local = int(np.argmin(np.abs(depths - target_depth)))
    return local, int(model_levels[local]), float(depths[local])


def fill(target: np.ndarray, value: np.ndarray, x0: int, y0: int) -> None:
    ny, nx = value.shape
    target[y0 : y0 + ny, x0 : x0 + nx] = value


def centered_from_staggered_x(u: np.ndarray) -> np.ndarray:
    return 0.5 * (u[:, :-1] + u[:, 1:])


def centered_from_staggered_y(v: np.ndarray) -> np.ndarray:
    return 0.5 * (v[:-1, :] + v[1:, :])


def finite_limits(a: np.ndarray, symmetric: bool = False) -> tuple[float, float]:
    vals = a[np.isfinite(a)]
    if vals.size == 0:
        return -1.0, 1.0
    if symmetric:
        lim = float(np.nanpercentile(np.abs(vals), 99))
        return -lim, lim
    return tuple(float(x) for x in np.nanpercentile(vals, [1, 99]))


def tile_offset(ds: nc.Dataset, axis: str) -> int:
    for name in (f"tBase{axis}", f"exch2_tBase{axis}"):
        if name in ds.ncattrs():
            return int(getattr(ds, name))
    global_name = f"exch2_t{axis}Globalo"
    if global_name in ds.ncattrs():
        return int(getattr(ds, global_name)) - 1
    raise AttributeError(f"No tile base attribute found for {axis}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--prefix", default="dyn3d_avg")
    parser.add_argument("--time-index", type=int, default=-1)
    parser.add_argument("--depth", type=float, default=500.0)
    parser.add_argument("--bathy", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    files = sorted(args.run_dir.glob(f"mnc_*_*/*{args.prefix}*.nc"))
    files = [p for p in files if not p.name.startswith("._")]
    if not files:
        raise FileNotFoundError(f"No dyn3d_avg MNC files found under {args.run_dir}")

    rho = np.full((612, 512), np.nan)
    u_stag = np.full((612, 513), np.nan)
    v_stag = np.full((613, 512), np.nan)
    lon = np.full(512, np.nan)
    lat = np.full(612, np.nan)
    chosen = None
    time_value = None
    iter_value = None

    for path in files:
        with nc.Dataset(path) as ds:
            if "RHOAnoma" not in ds.variables:
                continue
            local_k, model_k, depth = choose_index(ds.variables["diag_levels"][:], args.depth)
            time_len = len(ds.dimensions["T"])
            time_idx = args.time_index if args.time_index >= 0 else time_len + args.time_index
            if not 0 <= time_idx < time_len:
                raise IndexError(f"time index {args.time_index} is out of range for {path}")

            x0 = tile_offset(ds, "x")
            y0 = tile_offset(ds, "y")
            if chosen is None:
                chosen = (model_k, depth)
                time_value = float(ds.variables["T"][time_idx])
                iter_value = int(ds.variables["iter"][time_idx])
            fill(rho, ds.variables["RHOAnoma"][time_idx, local_k, :, :], x0, y0)
            fill(u_stag, ds.variables["UVEL"][time_idx, local_k, :, :], x0, y0)
            fill(v_stag, ds.variables["VVEL"][time_idx, local_k, :, :], x0, y0)
            lon[x0 : x0 + len(ds.variables["X"][:])] = ds.variables["X"][:]
            lat[y0 : y0 + len(ds.variables["Y"][:])] = ds.variables["Y"][:]

    u = centered_from_staggered_x(u_stag)
    v = centered_from_staggered_y(v_stag)
    lon = 238.0 + (np.arange(512) + 0.5) * 0.0156
    lat = 27.0 + (np.arange(612) + 0.5) * 0.0136
    if args.bathy:
        bathy = np.fromfile(args.bathy, ">f8").reshape(612, 512)
        wet = bathy <= -depth
        u = np.where(wet, u, np.nan)
        v = np.where(wet, v, np.nan)
        rho = np.where(wet, rho, np.nan)
    xx, yy = np.meshgrid(lon, lat)
    model_k, depth = chosen

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), constrained_layout=True)
    fields = [
        ("UVEL centered (m/s)", u, "RdBu_r", True),
        ("VVEL centered (m/s)", v, "RdBu_r", True),
        ("rhoNil + RHOAnoma (kg/m^3)", 999.8 + rho, "viridis", False),
    ]
    for ax, (title, field, cmap, symmetric) in zip(axes, fields, strict=True):
        vmin, vmax = finite_limits(field, symmetric=symmetric)
        mesh = ax.pcolormesh(xx, yy, field, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        fig.colorbar(mesh, ax=ax, shrink=0.85)
        ax.set_title(title)
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
        ax.set_aspect("equal", adjustable="box")

    fig.suptitle(
        f"{args.run_dir.name}: iter={iter_value}, t={time_value / 3600:.2f} h, "
        f"level={model_k}, z={depth:.1f} m"
    )
    output = args.output or args.run_dir / f"slice_500m_iter{iter_value}.png"
    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
