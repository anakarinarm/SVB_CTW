#!/usr/bin/env python3
"""Make a quick-look Rplus animation from MITgcm MNC tile output."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import cmocean as cmo
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import TwoSlopeNorm
from netCDF4 import Dataset


NX, NY = 512, 612
TILE_RE = re.compile(r"\.t(\d{3})\.nc$")
DELR = np.array(
    [
        2, 2, 2, 2, 2, 2.2, 2.4, 2.6, 2.8, 3.0,
        3.2, 3.5, 3.8, 4.2, 4.5, 4.8, 5, 5, 5, 5,
        5, 5, 5.5, 5.5, 5.5, 5.5, 5.5, 6, 6.5, 7,
        7.5, 8, 8.5, 9.2, 10, 10, 10, 10, 10, 10,
        10, 11, 12, 13, 14, 15, 16, 17, 18, 20,
        22, 23, 25, 25, 25, 25, 25, 25, 25, 25,
        27, 28, 30, 30, 32, 34, 36, 38, 40, 43,
        46, 50, 55, 59, 65, 70, 75, 80, 85, 92,
        100, 100, 100, 100, 100, 100, 105, 115, 125, 132,
        145, 160, 170, 180, 190, 200, 210, 230, 245, 269.8,
    ],
    dtype=float,
)
ZC = np.cumsum(DELR) - 0.5 * DELR


def tile_id(path: Path) -> int:
    match = TILE_RE.search(path.name)
    if not match:
        raise ValueError(f"cannot parse tile id from {path}")
    return int(match.group(1))


def chunk_key(path: Path, stem: str) -> str:
    name = path.name
    prefix = f"{stem}."
    if not name.startswith(prefix):
        return "single"
    rest = name[len(prefix):]
    parts = rest.split(".")
    return parts[0] if len(parts) > 2 else "single"


def open_stream(root: Path, stem: str) -> tuple[dict[int, dict[str, Dataset]], dict[str, np.ndarray]]:
    files = sorted(p for p in root.glob(f"mnc_Rplus_*/{stem}*.nc") if not p.name.startswith("._"))
    by_tile: dict[int, dict[str, Path]] = {}
    for path in files:
        by_tile.setdefault(tile_id(path), {})[chunk_key(path, stem)] = path
    first_by_key: dict[str, Path] = {}
    for paths in by_tile.values():
        for key, path in paths.items():
            first_by_key.setdefault(key, path)
    times_by_key = {}
    for key, path in first_by_key.items():
        with Dataset(str(path)) as ds:
            times_by_key[key] = np.asarray(ds.variables["T"][:], dtype=float) if "T" in ds.variables else np.array([])
    stream = {}
    for tid, paths in by_tile.items():
        stream[tid] = {key: Dataset(str(path)) for key, path in paths.items()}
    return stream, times_by_key


def close_stream(stream: dict[int, dict[str, Dataset]]) -> None:
    for datasets in stream.values():
        for ds in datasets.values():
            ds.close()


def read_xy(ds: Dataset) -> tuple[slice, slice]:
    x = np.asarray(ds.variables["X"][:], dtype=int)
    y = np.asarray(ds.variables["Y"][:], dtype=int)
    return slice(int(y[0]) - 1, int(y[-1])), slice(int(x[0]) - 1, int(x[-1]))


def dataset_for_time(datasets: dict[str, Dataset], times_by_key: dict[str, np.ndarray], time_s: float) -> tuple[Dataset, int]:
    for key, times in times_by_key.items():
        hit = np.nonzero(np.isclose(times, time_s))[0]
        if hit.size:
            return datasets[key], int(hit[0])
    raise ValueError(f"time {time_s} not found in tile stream")


def mosaic_grid(grid: dict[int, dict[str, Dataset]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lon = np.full((NY, NX), np.nan)
    lat = np.full((NY, NX), np.nan)
    hfac = np.full((NY, NX), np.nan)
    for datasets in grid.values():
        ds = next(iter(datasets.values()))
        ys, xs = read_xy(ds)
        lon[ys, xs] = ds.variables["xC"][:]
        lat[ys, xs] = ds.variables["yC"][:]
        hfac[ys, xs] = ds.variables["hFacC"][LEVEL_INDEX, :, :]
    return lon, lat, hfac > 0


def mosaic_etan(stream: dict[int, dict[str, Dataset]], times_by_key: dict[str, np.ndarray], time_s: float) -> np.ma.MaskedArray:
    out = np.full((NY, NX), np.nan)
    for datasets in stream.values():
        ds, it = dataset_for_time(datasets, times_by_key, time_s)
        ys, xs = read_xy(ds)
        out[ys, xs] = ds.variables["ETAN"][it, 0, :, :]
    return np.ma.masked_invalid(out)


def mosaic_dyn(stream: dict[int, dict[str, Dataset]], times_by_key: dict[str, np.ndarray], time_s: float) -> tuple[np.ma.MaskedArray, np.ma.MaskedArray, np.ma.MaskedArray]:
    u = np.full((NY, NX), np.nan)
    v = np.full((NY, NX), np.nan)
    w = np.full((NY, NX), np.nan)
    for datasets in stream.values():
        ds, it = dataset_for_time(datasets, times_by_key, time_s)
        ys, xs = read_xy(ds)
        u_edge = ds.variables["UVEL"][it, LEVEL_IN_FILE, :, :]
        v_edge = ds.variables["VVEL"][it, LEVEL_IN_FILE, :, :]
        u[ys, xs] = 0.5 * (u_edge[:, :-1] + u_edge[:, 1:])
        v[ys, xs] = 0.5 * (v_edge[:-1, :] + v_edge[1:, :])
        w[ys, xs] = ds.variables["WVEL"][it, LEVEL_IN_FILE, :, :]
    return np.ma.masked_invalid(u), np.ma.masked_invalid(v), np.ma.masked_invalid(w)


def robust_limits(
    surf: dict[int, dict[str, Dataset]],
    surf_times: dict[str, np.ndarray],
    dyn: dict[int, dict[str, Dataset]],
    dyn_times: dict[str, np.ndarray],
    times: np.ndarray,
    mask: np.ndarray,
) -> tuple[float, float, float]:
    etan_vals = []
    w_vals = []
    spd_vals = []
    sample_times = times[:: max(1, len(times) // 24)]
    for time_s in sample_times:
        etan = np.ma.masked_where(~mask, mosaic_etan(surf, surf_times, float(time_s)))
        u, v, w = mosaic_dyn(dyn, dyn_times, float(time_s))
        u = np.ma.masked_where(~mask, u)
        v = np.ma.masked_where(~mask, v)
        w = np.ma.masked_where(~mask, w)
        spd = np.sqrt(u * u + v * v)
        etan_vals.append(np.asarray(etan.compressed()))
        w_vals.append(np.asarray(w.compressed()))
        spd_vals.append(np.asarray(spd.compressed()))
    etan_all = np.concatenate(etan_vals)
    w_all = np.concatenate(w_vals)
    spd_all = np.concatenate(spd_vals)
    etan_lim = float(np.nanpercentile(np.abs(etan_all), 99.5))
    w_lim = float(np.nanpercentile(np.abs(w_all), 99.5))
    spd_lim = float(np.nanpercentile(spd_all, 99.5))
    return max(etan_lim, 1e-6), max(w_lim, 1e-8), max(spd_lim, 1e-6)


def make_animation(args: argparse.Namespace) -> None:
    root = args.root
    out = args.output
    out.parent.mkdir(parents=True, exist_ok=True)

    surf, surf_times_by_key = open_stream(root, "surf_avg")
    dyn, dyn_times_by_key = open_stream(root, "dyn3d_avg")
    grid, _grid_times_by_key = open_stream(root, "grid")
    try:
        global LEVEL_INDEX, LEVEL_IN_FILE
        first_dyn = next(iter(next(iter(dyn.values())).values()))
        saved_levels = np.asarray(first_dyn.variables["diag_levels"][:], dtype=int)
        available_depths = ZC[saved_levels - 1]
        LEVEL_IN_FILE = int(np.argmin(np.abs(available_depths - args.depth_m)))
        LEVEL_INDEX = int(saved_levels[LEVEL_IN_FILE] - 1)
        depth_m = float(ZC[LEVEL_INDEX])

        surf_times = np.concatenate(list(surf_times_by_key.values()))
        dyn_times = np.concatenate(list(dyn_times_by_key.values()))
        common = np.intersect1d(surf_times, dyn_times)
        times = common[2::3]
        times = times[times <= args.days * 86400]

        print(f"using {len(times)} hourly frames from {times[0]} to {times[-1]} s", flush=True)
        print(f"requested depth {args.depth_m:g} m; using saved model level {LEVEL_INDEX + 1} at {depth_m:.1f} m", flush=True)

        lon, lat, mask = mosaic_grid(grid)
        extent = [float(np.nanmin(lon)), float(np.nanmax(lon)), float(np.nanmin(lat)), float(np.nanmax(lat))]

        etan_lim, w_lim, spd_lim = robust_limits(surf, surf_times_by_key, dyn, dyn_times_by_key, times, mask)
        print(f"color limits: ETAN +/-{etan_lim:.4g} m, WVEL +/-{w_lim:.4g} m/s, speed 0-{spd_lim:.4g} m/s", flush=True)

        fig, axes = plt.subplots(1, 3, figsize=(16, 6), constrained_layout=True)
        fig.patch.set_facecolor("white")
        for ax in axes:
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")

        etan0 = np.ma.masked_where(~mask, mosaic_etan(surf, surf_times_by_key, float(times[0])))
        u0, v0, w0 = mosaic_dyn(dyn, dyn_times_by_key, float(times[0]))
        u0 = np.ma.masked_where(~mask, u0)
        v0 = np.ma.masked_where(~mask, v0)
        w0 = np.ma.masked_where(~mask, w0)
        speed0 = np.sqrt(u0 * u0 + v0 * v0)

        im0 = axes[0].imshow(etan0, origin="lower", extent=extent, cmap=cmo.cm.curl, norm=TwoSlopeNorm(vcenter=0, vmin=-etan_lim, vmax=etan_lim))
        im1 = axes[1].imshow(w0, origin="lower", extent=extent, cmap=cmo.cm.balance, norm=TwoSlopeNorm(vcenter=0, vmin=-w_lim, vmax=w_lim))
        im2 = axes[2].imshow(speed0, origin="lower", extent=extent, cmap=cmo.cm.speed, vmin=0, vmax=spd_lim)

        qstep = args.quiver_step
        qlon = lon[::qstep, ::qstep]
        qlat = lat[::qstep, ::qstep]
        q = axes[2].quiver(qlon, qlat, u0[::qstep, ::qstep], v0[::qstep, ::qstep], color="k", scale=args.quiver_scale, width=0.0022, alpha=0.75)

        cb0 = fig.colorbar(im0, ax=axes[0], shrink=0.82)
        cb1 = fig.colorbar(im1, ax=axes[1], shrink=0.82)
        cb2 = fig.colorbar(im2, ax=axes[2], shrink=0.82)
        cb0.set_label("ETAN (m)")
        cb1.set_label("WVEL (m s$^{-1}$)")
        cb2.set_label("Speed (m s$^{-1}$)")

        axes[0].set_title("SSH / ETAN")
        axes[1].set_title(f"Vertical velocity, level {LEVEL_INDEX + 1} ({depth_m:.1f} m)")
        axes[2].set_title(f"Horizontal speed + velocity, level {LEVEL_INDEX + 1} ({depth_m:.1f} m)")
        title = fig.suptitle("")

        writer = FFMpegWriter(fps=args.fps, metadata={"title": "Rplus MITgcm QC animation"})
        with writer.saving(fig, str(out), dpi=args.dpi):
            for i, time_s in enumerate(times):
                etan = np.ma.masked_where(~mask, mosaic_etan(surf, surf_times_by_key, float(time_s)))
                u, v, w = mosaic_dyn(dyn, dyn_times_by_key, float(time_s))
                u = np.ma.masked_where(~mask, u)
                v = np.ma.masked_where(~mask, v)
                w = np.ma.masked_where(~mask, w)
                speed = np.sqrt(u * u + v * v)

                im0.set_data(etan)
                im1.set_data(w)
                im2.set_data(speed)
                q.set_UVC(u[::qstep, ::qstep], v[::qstep, ::qstep])
                title.set_text(f"Rplus QC, day {time_s / 86400:.2f} ({time_s / 3600:.1f} h)")
                writer.grab_frame()
                if (i + 1) % 24 == 0 or i == len(times) - 1:
                    print(f"rendered {i + 1}/{len(times)} frames", flush=True)

        plt.close(fig)
        print(f"wrote {out}", flush=True)
    finally:
        close_stream(surf)
        close_stream(dyn)
        close_stream(grid)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("/Volumes/LaCie/CTW_SVB_reruns/Rplus"))
    parser.add_argument("--output", type=Path, default=Path("/Volumes/LaCie/CTW_SVB_reruns/Rplus/Rplus_qc_hourly_etan_wvel_speed.mp4"))
    parser.add_argument("--depth-m", type=float, default=500.0)
    parser.add_argument("--days", type=float, default=15.0)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--dpi", type=int, default=130)
    parser.add_argument("--quiver-step", type=int, default=18)
    parser.add_argument("--quiver-scale", type=float, default=4.0)
    return parser.parse_args()


if __name__ == "__main__":
    make_animation(parse_args())
