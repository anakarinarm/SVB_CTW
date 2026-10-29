#!/usr/bin/env python3
"""Render a six-panel active-domain R/S50 difference animation."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import cmocean as cmo
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from MITgcmutils import mds
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import TwoSlopeNorm


DELR = np.array([2, 2, 2, 2, 2, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.5, 3.8, 4.2, 4.5, 4.8, 5, 5, 5, 5, 5, 5, 5.5, 5.5, 5.5, 5.5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9.2, 10, 10, 10, 10, 10, 10, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 22, 23, 25, 25, 25, 25, 25, 25, 25, 25, 27, 28, 30, 30, 32, 34, 36, 38, 40, 43, 46, 50, 55, 59, 65, 70, 75, 80, 85, 92, 100, 100, 100, 100, 100, 100, 105, 115, 125, 132, 145, 160, 170, 180, 190, 200, 210, 230, 245, 269.8], dtype=float)
SAVED_LEVELS = np.array(list(range(1, 66)) + [67, 69, 71, 73, 75, 77, 79, 84, 90, 100])
ZC = np.cumsum(DELR) - DELR / 2
ITER_RE = re.compile(r"\.(\d{10})\.meta$")
ACTIVE_Y = slice(72, 684)
ACTIVE_X = slice(64, 576)


def complete(root: Path, stem: str, size: int) -> set[int]:
    result = set()
    for meta in root.glob(f"{stem}.*.meta"):
        match = ITER_RE.search(meta.name)
        data = meta.with_suffix(".data")
        if match and data.exists() and data.stat().st_size == size:
            result.add(int(match.group(1)))
    return result


def surf(root: Path, iteration: int, ny: int, nx: int) -> np.ndarray:
    data = np.memmap(root / f"surf_avg.{iteration:010d}.data", dtype=">f8", mode="r", shape=(2, ny, nx))
    return np.asarray(data[0]).copy()


def dyn(root: Path, iteration: int, level: int, ny: int, nx: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.memmap(root / f"dyn3d_avg.{iteration:010d}.data", dtype=">f8", mode="r", shape=(5, len(SAVED_LEVELS), ny, nx))
    return tuple(np.asarray(data[index, level]).copy() for index in (0, 1, 2))


def mask(values: np.ndarray, wet: np.ndarray) -> np.ma.MaskedArray:
    return np.ma.masked_where(~wet | ~np.isfinite(values), values)


def limits(groups: list[tuple[Path, Path]], iters: list[int], masks: tuple[np.ndarray, np.ndarray], level: int, ny: int, nx: int, ys: slice, xs: slice) -> tuple[float, float, float]:
    ssh, wvel, speed = [], [], []
    for iteration in iters[::max(1, len(iters) // 24)]:
        for plus, minus in groups:
            ds = 1e3 * (surf(plus, iteration, ny, nx) - surf(minus, iteration, ny, nx))[ys, xs]
            up, vp, wp = dyn(plus, iteration, level, ny, nx)
            um, vm, wm = dyn(minus, iteration, level, ny, nx)
            du = mask((up - um)[ys, xs], masks[1])
            dv = mask((vp - vm)[ys, xs], masks[1])
            dw = mask((wp - wm)[ys, xs], masks[1])
            ssh.append(mask(ds, masks[0]).compressed())
            wvel.append(dw.compressed())
            speed.append((1e3 * np.ma.sqrt(du * du + dv * dv)).compressed())
    return tuple(max(float(np.nanpercentile(np.abs(np.concatenate(values)) if index < 2 else np.concatenate(values), 99.5)), 1e-10) for index, values in enumerate((ssh, wvel, speed)))


def main(args: argparse.Namespace) -> None:
    groups = [(args.rplus, args.rminus), (args.s50plus, args.s50minus)]
    lon = mds.rdmds(str(args.rplus / "XC"))
    lat = mds.rdmds(str(args.rplus / "YC"))
    hfac = [mds.rdmds(str(root / "hFacC")) for pair in groups for root in pair]
    ny, nx = lon.shape
    active = np.logical_and.reduce([field[0] > 0 for field in hfac])
    ys, xs = ACTIVE_Y, ACTIVE_X
    lon, lat = lon[ys, xs], lat[ys, xs]
    surface_mask = active[ys, xs]
    level_index = int(np.argmin(np.abs(ZC[SAVED_LEVELS - 1] - args.depth_m)))
    level_mask = np.logical_and.reduce([field[SAVED_LEVELS[level_index] - 1] > 0 for field in hfac])[ys, xs]
    surf_bytes, dyn_bytes = 2 * ny * nx * 8, 5 * len(SAVED_LEVELS) * ny * nx * 8
    common = None
    for plus, minus in groups:
        run_iters = complete(plus, "surf_avg", surf_bytes) & complete(plus, "dyn3d_avg", dyn_bytes) & complete(minus, "surf_avg", surf_bytes) & complete(minus, "dyn3d_avg", dyn_bytes)
        common = run_iters if common is None else common & run_iters
    iters = [it for it in sorted(common) if it * args.delta_t <= args.days * 86400]
    if args.max_frames:
        iters = iters[:args.max_frames]
    if not iters:
        raise RuntimeError("no complete common diagnostic records found")
    ssh_lim, w_lim, speed_lim = limits(groups, iters, (surface_mask, level_mask), level_index, ny, nx, ys, xs)
    depth = ZC[SAVED_LEVELS[level_index] - 1]
    extent = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True, sharex=True, sharey=True)
    fig.patch.set_facecolor("white")
    images, quivers = [], []
    for row, (plus, minus) in enumerate(groups):
        ds = mask(1e3 * (surf(plus, iters[0], ny, nx) - surf(minus, iters[0], ny, nx))[ys, xs], surface_mask)
        up, vp, wp = dyn(plus, iters[0], level_index, ny, nx)
        um, vm, wm = dyn(minus, iters[0], level_index, ny, nx)
        du, dv = mask((up - um)[ys, xs], level_mask), mask((vp - vm)[ys, xs], level_mask)
        dw = mask((wp - wm)[ys, xs], level_mask)
        fields = [(ds, cmo.cm.curl.copy(), TwoSlopeNorm(vcenter=0, vmin=-ssh_lim, vmax=ssh_lim)), (dw, cmo.cm.balance.copy(), TwoSlopeNorm(vcenter=0, vmin=-w_lim, vmax=w_lim)), (1e3 * np.ma.sqrt(du * du + dv * dv), cmo.cm.speed.copy(), None)]
        for col, (field, cmap, norm) in enumerate(fields):
            cmap.set_bad("tan")
            image = axes[row, col].imshow(field, origin="lower", extent=extent, cmap=cmap, norm=norm, vmin=None if norm else 0, vmax=None if norm else speed_lim)
            images.append(image)
            axes[row, col].set_title(("SSH difference", f"WVEL difference at {depth:.1f} m", f"Horizontal speed difference + quivers at {depth:.1f} m")[col])
            axes[row, col].set_aspect("auto")
            axes[row, col].set_facecolor("tan")
            fig.colorbar(image, ax=axes[row, col], shrink=0.82).set_label(("SSH (mm)", "WVEL (m s$^{-1}$)", "speed (mm s$^{-1}$)")[col])
        step = args.quiver_step
        quivers.append(axes[row, 2].quiver(lon[::step, ::step], lat[::step, ::step], 1e3 * du[::step, ::step], 1e3 * dv[::step, ::step], color="k", scale=args.quiver_scale, width=0.0035, headwidth=4.5, headlength=5.5, headaxislength=4.5, alpha=0.78))
        axes[row, 0].text(0.03, 0.95, ("Rplus - Rminus", "S50plus - S50minus")[row], transform=axes[row, 0].transAxes, ha="left", va="top", bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "0.25", "pad": 4})
    for ax in axes[:, 0]:
        ax.set_xlabel("longitude")
    for ax in axes[1, :]:
        ax.set_xlabel("latitude")
    title = fig.suptitle("")
    writer = FFMpegWriter(fps=args.fps, metadata={"title": "R/S50 active-domain difference animation"})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with writer.saving(fig, str(args.output), dpi=args.dpi):
        for index, iteration in enumerate(iters):
            for row, (plus, minus) in enumerate(groups):
                ds = mask(1e3 * (surf(plus, iteration, ny, nx) - surf(minus, iteration, ny, nx))[ys, xs], surface_mask)
                up, vp, wp = dyn(plus, iteration, level_index, ny, nx)
                um, vm, wm = dyn(minus, iteration, level_index, ny, nx)
                du, dv = mask((up - um)[ys, xs], level_mask), mask((vp - vm)[ys, xs], level_mask)
                dw = mask((wp - wm)[ys, xs], level_mask)
                images[row * 3].set_data(ds); images[row * 3 + 1].set_data(dw); images[row * 3 + 2].set_data(1e3 * np.ma.sqrt(du * du + dv * dv)); quivers[row].set_UVC(1e3 * du[::args.quiver_step, ::args.quiver_step], 1e3 * dv[::args.quiver_step, ::args.quiver_step])
            title.set_text(f"Active-domain differences, day {iteration * args.delta_t / 86400:.2f}")
            writer.grab_frame()
            if (index + 1) % 60 == 0 or index == len(iters) - 1:
                print(f"rendered {index + 1}/{len(iters)} frames", flush=True)
    plt.close(fig)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rplus", type=Path, required=True); parser.add_argument("--rminus", type=Path, required=True); parser.add_argument("--s50plus", type=Path, required=True); parser.add_argument("--s50minus", type=Path, required=True); parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--depth-m", type=float, default=500.0); parser.add_argument("--days", type=float, default=15.0); parser.add_argument("--delta-t", type=float, default=60.0); parser.add_argument("--fps", type=int, default=12); parser.add_argument("--dpi", type=int, default=120); parser.add_argument("--quiver-step", type=int, default=20); parser.add_argument("--quiver-scale", type=float, default=180.0); parser.add_argument("--max-frames", type=int, default=0)
    main(parser.parse_args())