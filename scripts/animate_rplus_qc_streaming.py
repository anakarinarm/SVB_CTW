#!/usr/bin/env python3
"""Make a streaming quick-look Rplus animation from MITgcm MNC tile output."""

from __future__ import annotations

import argparse
import re
from collections import OrderedDict
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
    if stem == "grid":
        return "grid"
    rest = path.name[len(stem) + 1 :]
    parts = rest.split(".")
    return parts[0] if len(parts) > 2 else "single"


def read_xy(ds: Dataset) -> tuple[slice, slice]:
    nx = int(getattr(ds, "sNx"))
    ny = int(getattr(ds, "sNy"))
    tiles_x = int(getattr(ds, "Nx")) // nx
    tile = int(getattr(ds, "tile_number")) - 1
    x0 = (tile % tiles_x) * nx
    y0 = (tile // tiles_x) * ny
    return slice(y0, y0 + ny), slice(x0, x0 + nx)


class DatasetCache:
    def __init__(self, max_open: int = 256) -> None:
        self.max_open = max_open
        self._cache: OrderedDict[Path, Dataset] = OrderedDict()

    def get(self, path: Path) -> Dataset:
        ds = self._cache.pop(path, None)
        if ds is None:
            ds = Dataset(str(path))
        self._cache[path] = ds
        while len(self._cache) > self.max_open:
            _, old = self._cache.popitem(last=False)
            old.close()
        return ds

    def close(self) -> None:
        for ds in self._cache.values():
            ds.close()
        self._cache.clear()


def build_index(root: Path, stem: str) -> dict[int, dict[str, dict[str, object]]]:
    files = sorted(p for p in root.glob(f"mnc_*/{stem}*.nc") if not p.name.startswith("._"))
    index: dict[int, dict[str, dict[str, object]]] = {}
    first_by_chunk: dict[str, Path] = {}
    for path in files:
        chunk = chunk_key(path, stem)
        index.setdefault(tile_id(path), {})[chunk] = {"path": path}
        first_by_chunk.setdefault(chunk, path)
    if not index:
        raise FileNotFoundError(f"no files found for {stem} under {root}")
    times_by_chunk = {}
    for chunk, path in first_by_chunk.items():
        with Dataset(str(path)) as ds:
            times_by_chunk[chunk] = np.asarray(ds.variables["T"][:], dtype=float) if "T" in ds.variables else np.array([])
    for chunks in index.values():
        for chunk, item in chunks.items():
            item["times"] = times_by_chunk[chunk]
    return index


def find_time(item: dict[str, object], time_s: float) -> tuple[Path, int] | None:
    times = item["times"]
    hit = np.nonzero(np.isclose(times, time_s))[0]
    if hit.size:
        return item["path"], int(hit[0])
    return None


def tile_for_time(chunks: dict[str, dict[str, object]], time_s: float) -> tuple[Path, int]:
    for item in chunks.values():
        found = find_time(item, time_s)
        if found is not None:
            return found
    raise ValueError(f"time {time_s} not found")


def complete_times(*indices: dict[int, dict[str, dict[str, object]]]) -> np.ndarray:
    complete: set[float] | None = None
    for index in indices:
        for chunks in index.values():
            tile_times = set()
            for item in chunks.values():
                tile_times.update(float(t) for t in item["times"])
            complete = tile_times if complete is None else complete.intersection(tile_times)
    return np.array(sorted(complete or []), dtype=float)


def mosaic_grid(grid_index: dict[int, dict[str, dict[str, object]]], level_index: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lon = np.full((NY, NX), np.nan)
    lat = np.full((NY, NX), np.nan)
    hfac = np.full((NY, NX), np.nan)
    for chunks in grid_index.values():
        path = next(iter(chunks.values()))["path"]
        with Dataset(str(path)) as ds:
            ys, xs = read_xy(ds)
            lon[ys, xs] = ds.variables["xC"][:]
            lat[ys, xs] = ds.variables["yC"][:]
            hfac[ys, xs] = ds.variables["hFacC"][level_index, :, :]
    return lon, lat, hfac > 0


def mosaic_etan(index: dict[int, dict[str, dict[str, object]]], cache: DatasetCache, time_s: float) -> np.ma.MaskedArray:
    out = np.full((NY, NX), np.nan)
    for chunks in index.values():
        path, it = tile_for_time(chunks, time_s)
        ds = cache.get(path)
        ys, xs = read_xy(ds)
        out[ys, xs] = ds.variables["ETAN"][it, 0, :, :]
    return np.ma.masked_invalid(out)


def mosaic_dyn(
    index: dict[int, dict[str, dict[str, object]]],
    cache: DatasetCache,
    time_s: float,
    level_in_file: int,
) -> tuple[np.ma.MaskedArray, np.ma.MaskedArray, np.ma.MaskedArray]:
    u = np.full((NY, NX), np.nan)
    v = np.full((NY, NX), np.nan)
    w = np.full((NY, NX), np.nan)
    for chunks in index.values():
        path, it = tile_for_time(chunks, time_s)
        ds = cache.get(path)
        ys, xs = read_xy(ds)
        u_edge = ds.variables["UVEL"][it, level_in_file, :, :]
        v_edge = ds.variables["VVEL"][it, level_in_file, :, :]
        u[ys, xs] = 0.5 * (u_edge[:, :-1] + u_edge[:, 1:])
        v[ys, xs] = 0.5 * (v_edge[:-1, :] + v_edge[1:, :])
        w[ys, xs] = ds.variables["WVEL"][it, level_in_file, :, :]
    return np.ma.masked_invalid(u), np.ma.masked_invalid(v), np.ma.masked_invalid(w)


def robust_limits(
    surf: dict[int, dict[str, dict[str, object]]],
    dyn: dict[int, dict[str, dict[str, object]]],
    times: np.ndarray,
    mask: np.ndarray,
    level_in_file: int,
    cache: DatasetCache,
) -> tuple[float, float, float]:
    etan_vals = []
    w_vals = []
    spd_vals = []
    sample_times = times[:: max(1, len(times) // 24)]
    for time_s in sample_times:
        etan = np.ma.masked_where(~mask, mosaic_etan(surf, cache, float(time_s)))
        u, v, w = mosaic_dyn(dyn, cache, float(time_s), level_in_file)
        u = np.ma.masked_where(~mask, u)
        v = np.ma.masked_where(~mask, v)
        w = np.ma.masked_where(~mask, w)
        etan_vals.append(np.asarray(etan.compressed()))
        w_vals.append(np.asarray(w.compressed()))
        spd_vals.append(np.asarray(np.sqrt(u * u + v * v).compressed()))
    etan_lim = float(np.nanpercentile(np.concatenate(etan_vals), 99.5))
    w_lim = float(np.nanpercentile(np.abs(np.concatenate(w_vals)), 99.5))
    spd_lim = float(np.nanpercentile(np.concatenate(spd_vals), 99.5))
    return max(abs(etan_lim), 1e-6), max(w_lim, 1e-8), max(spd_lim, 1e-6)


def make_animation(args: argparse.Namespace) -> None:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    surf = build_index(args.root, "surf_avg")
    dyn = build_index(args.root, "dyn3d_avg")
    grid = build_index(args.root, "grid")
    cache = DatasetCache(args.max_open)
    try:
        first_dyn_path = next(iter(next(iter(dyn.values())).values()))["path"]
        with Dataset(str(first_dyn_path)) as ds:
            saved_levels = np.asarray(ds.variables["diag_levels"][:], dtype=int)
        available_depths = ZC[saved_levels - 1]
        level_in_file = int(np.argmin(np.abs(available_depths - args.depth_m)))
        level_index = int(saved_levels[level_in_file] - 1)
        depth_m = float(ZC[level_index])

        common = complete_times(surf, dyn)
        hourly = common[2::3] if len(common) >= 3 else common
        times = hourly[hourly <= args.days * 86400]
        if args.max_frames and len(times) > args.max_frames:
            times = times[: args.max_frames]
        if len(times) == 0:
            raise RuntimeError("no common surface/3-D times available for animation")

        print(f"using {len(times)} hourly frames from {times[0]} to {times[-1]} s", flush=True)
        print(f"requested depth {args.depth_m:g} m; using saved model level {level_index + 1} at {depth_m:.1f} m", flush=True)

        lon, lat, mask = mosaic_grid(grid, level_index)
        extent = [float(np.nanmin(lon)), float(np.nanmax(lon)), float(np.nanmin(lat)), float(np.nanmax(lat))]
        etan_lim, w_lim, spd_lim = robust_limits(surf, dyn, times, mask, level_in_file, cache)
        print(f"color limits: ETAN +/-{etan_lim:.4g} m, WVEL +/-{w_lim:.4g} m/s, speed 0-{spd_lim:.4g} m/s", flush=True)

        etan0 = np.ma.masked_where(~mask, mosaic_etan(surf, cache, float(times[0])))
        u0, v0, w0 = mosaic_dyn(dyn, cache, float(times[0]), level_in_file)
        u0 = np.ma.masked_where(~mask, u0)
        v0 = np.ma.masked_where(~mask, v0)
        w0 = np.ma.masked_where(~mask, w0)
        speed0 = np.sqrt(u0 * u0 + v0 * v0)

        fig, axes = plt.subplots(1, 3, figsize=(16, 6), constrained_layout=True)
        for ax in axes:
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")

        im0 = axes[0].imshow(etan0, origin="lower", extent=extent, cmap=cmo.cm.curl, norm=TwoSlopeNorm(vcenter=0, vmin=-etan_lim, vmax=etan_lim))
        im1 = axes[1].imshow(w0, origin="lower", extent=extent, cmap=cmo.cm.balance, norm=TwoSlopeNorm(vcenter=0, vmin=-w_lim, vmax=w_lim))
        im2 = axes[2].imshow(speed0, origin="lower", extent=extent, cmap=cmo.cm.speed, vmin=0, vmax=spd_lim)

        qstep = args.quiver_step
        q = axes[2].quiver(
            lon[::qstep, ::qstep],
            lat[::qstep, ::qstep],
            u0[::qstep, ::qstep],
            v0[::qstep, ::qstep],
            color="k",
            scale=args.quiver_scale,
            width=0.0022,
            alpha=0.75,
        )
        fig.colorbar(im0, ax=axes[0], shrink=0.82).set_label("ETAN (m)")
        fig.colorbar(im1, ax=axes[1], shrink=0.82).set_label("WVEL (m s$^{-1}$)")
        fig.colorbar(im2, ax=axes[2], shrink=0.82).set_label("Speed (m s$^{-1}$)")
        axes[0].set_title("SSH / ETAN")
        axes[1].set_title(f"Vertical velocity, level {level_index + 1} ({depth_m:.1f} m)")
        axes[2].set_title(f"Horizontal speed + velocity, level {level_index + 1} ({depth_m:.1f} m)")
        title = fig.suptitle("")

        writer = FFMpegWriter(fps=args.fps, metadata={"title": "Rplus telescoped MITgcm QC animation"})
        with writer.saving(fig, str(args.output), dpi=args.dpi):
            for i, time_s in enumerate(times):
                etan = np.ma.masked_where(~mask, mosaic_etan(surf, cache, float(time_s)))
                u, v, w = mosaic_dyn(dyn, cache, float(time_s), level_in_file)
                u = np.ma.masked_where(~mask, u)
                v = np.ma.masked_where(~mask, v)
                w = np.ma.masked_where(~mask, w)
                im0.set_data(etan)
                im1.set_data(w)
                im2.set_data(np.sqrt(u * u + v * v))
                q.set_UVC(u[::qstep, ::qstep], v[::qstep, ::qstep])
                title.set_text(f"Rplus telescoped QC, day {time_s / 86400:.2f} ({time_s / 3600:.1f} h)")
                writer.grab_frame()
                if (i + 1) % 24 == 0 or i == len(times) - 1:
                    print(f"rendered {i + 1}/{len(times)} frames", flush=True)
        plt.close(fig)
        print(f"wrote {args.output}", flush=True)
    finally:
        cache.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--depth-m", type=float, default=500.0)
    parser.add_argument("--days", type=float, default=15.0)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--dpi", type=int, default=130)
    parser.add_argument("--quiver-step", type=int, default=18)
    parser.add_argument("--quiver-scale", type=float, default=4.0)
    parser.add_argument("--max-open", type=int, default=256)
    parser.add_argument("--max-frames", type=int, default=0)
    return parser.parse_args()


if __name__ == "__main__":
    make_animation(parse_args())
