#!/usr/bin/env python3
"""Make a quick-look Rplus animation from MITgcm MDS binary diagnostics."""

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


DELR = np.array(
    [
        2, 2, 2, 2, 2, 2.2, 2.4, 2.6, 2.8, 3.0,
        3.2, 3.5, 3.8, 4.2, 4.5, 4.8, 5, 5, 5,
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
SAVED_LEVELS = np.array(
    list(range(1, 66)) + [67, 69, 71, 73, 75, 77, 79, 84, 90, 100],
    dtype=int,
)
ZC = np.cumsum(DELR) - 0.5 * DELR
ITER_RE = re.compile(r"\.(\d{10})\.meta$")


def completed_iters(root: Path, stem: str, expected_bytes: int) -> list[int]:
    out = []
    for meta in root.glob(f"{stem}.*.meta"):
        match = ITER_RE.search(meta.name)
        if not match:
            continue
        data = meta.with_suffix(".data")
        if data.exists() and data.stat().st_size == expected_bytes:
            out.append(int(match.group(1)))
    return sorted(set(out))


def read_surf(root: Path, iteration: int, ny: int, nx: int,
              yslice: slice = slice(None), xslice: slice = slice(None)) -> np.ndarray:
    path = root / f"surf_avg.{iteration:010d}.data"
    mm = np.memmap(path, dtype=">f8", mode="r", shape=(2, ny, nx))
    return np.asarray(mm[0, yslice, xslice]).copy()


def read_dyn_level(
    root: Path,
    iteration: int,
    level_in_file: int,
    ny: int,
    nx: int,
    yslice: slice = slice(None),
    xslice: slice = slice(None),
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = root / f"dyn3d_avg.{iteration:010d}.data"
    mm = np.memmap(path, dtype=">f8", mode="r", shape=(5, len(SAVED_LEVELS), ny, nx))
    u = np.asarray(mm[0, level_in_file, yslice, xslice]).copy()
    v = np.asarray(mm[1, level_in_file, yslice, xslice]).copy()
    w = np.asarray(mm[2, level_in_file, yslice, xslice]).copy()
    return u, v, w


def masked(a: np.ndarray, mask: np.ndarray) -> np.ma.MaskedArray:
    return np.ma.masked_where(~mask | ~np.isfinite(a), a)


def speed(u: np.ma.MaskedArray, v: np.ma.MaskedArray) -> np.ma.MaskedArray:
    return np.ma.sqrt(u * u + v * v)


def robust_limits(
    root: Path,
    iters: list[int],
    surface_mask: np.ndarray,
    level_mask: np.ndarray,
    level_in_file: int,
    ny: int,
    nx: int,
    yslice: slice = slice(None),
    xslice: slice = slice(None),
) -> tuple[float, float, float]:
    etan_vals = []
    w_vals = []
    spd_vals = []
    sample_iters = iters[:: max(1, len(iters) // 24)]
    for iteration in sample_iters:
        etan = masked(read_surf(root, iteration, ny, nx, yslice, xslice), surface_mask)
        u, v, w = read_dyn_level(root, iteration, level_in_file, ny, nx, yslice, xslice)
        u = masked(u, level_mask)
        v = masked(v, level_mask)
        w = masked(w, level_mask)
        etan_vals.append(np.asarray(etan.compressed()))
        w_vals.append(np.asarray(w.compressed()))
        spd_vals.append(np.asarray(speed(u, v).compressed()))

    etan_all = np.concatenate(etan_vals)
    w_all = np.concatenate(w_vals)
    spd_all = np.concatenate(spd_vals)
    etan_lim = float(np.nanpercentile(np.abs(etan_all), 99.5))
    w_lim = float(np.nanpercentile(np.abs(w_all), 99.5))
    spd_lim = float(np.nanpercentile(spd_all, 99.5))
    return max(etan_lim, 1e-6), max(w_lim, 1e-8), max(spd_lim, 1e-6)


def make_animation(args: argparse.Namespace) -> None:
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    run_label = args.label or args.root.name
    lon = mds.rdmds(str(args.root / "XC"))
    lat = mds.rdmds(str(args.root / "YC"))
    hfac = mds.rdmds(str(args.root / "hFacC"))
    ny, nx = lon.shape
    xend, yend = args.x_end or nx, args.y_end or ny
    if not (0 <= args.x_start < xend <= nx and 0 <= args.y_start < yend <= ny):
        raise ValueError("Invalid grid crop")
    xslice, yslice = slice(args.x_start, xend), slice(args.y_start, yend)
    lon, lat = lon[yslice, xslice], lat[yslice, xslice]
    print(f"plotting X={args.x_start}:{xend}, Y={args.y_start}:{yend}, shape={lon.shape}", flush=True)

    level_index = int(np.argmin(np.abs(ZC - args.depth_m)))
    saved_hit = np.nonzero(SAVED_LEVELS == level_index + 1)[0]
    if not saved_hit.size:
        level_in_file = int(np.argmin(np.abs(ZC[SAVED_LEVELS - 1] - args.depth_m)))
        level_index = int(SAVED_LEVELS[level_in_file] - 1)
    else:
        level_in_file = int(saved_hit[0])
    depth_m = float(ZC[level_index])
    surface_mask = hfac[0, yslice, xslice] > 0
    level_mask = hfac[level_index, yslice, xslice] > 0

    surf_bytes = 2 * ny * nx * 8
    dyn_bytes = 5 * len(SAVED_LEVELS) * ny * nx * 8
    iters = sorted(set(completed_iters(args.root, "surf_avg", surf_bytes)).intersection(
        completed_iters(args.root, "dyn3d_avg", dyn_bytes)
    ))
    iters = [it for it in iters if it * args.delta_t <= args.days * 86400]
    iters = [it for it in iters if (it * args.delta_t) % 3600 == 0]
    if args.max_frames and len(iters) > args.max_frames:
        iters = iters[: args.max_frames]
    if not iters:
        raise RuntimeError("no complete common surf_avg/dyn3d_avg records found")

    print(
        f"using {len(iters)} hourly frames from iteration {iters[0]} to {iters[-1]} "
        f"({iters[-1] * args.delta_t / 86400:.2f} days)",
        flush=True,
    )
    print(f"requested depth {args.depth_m:g} m; using model level {level_index + 1} at {depth_m:.1f} m", flush=True)
    etan_lim, w_lim, spd_lim = robust_limits(args.root, iters, surface_mask, level_mask, level_in_file, ny, nx, yslice, xslice)
    print(f"color limits: ETAN +/-{etan_lim:.4g} m, WVEL +/-{w_lim:.4g} m/s, speed 0-{spd_lim:.4g} m/s", flush=True)

    etan0 = masked(read_surf(args.root, iters[0], ny, nx, yslice, xslice), surface_mask)
    u0, v0, w0 = read_dyn_level(args.root, iters[0], level_in_file, ny, nx, yslice, xslice)
    u0 = masked(u0, level_mask)
    v0 = masked(v0, level_mask)
    w0 = masked(w0, level_mask)

    fig, axes = plt.subplots(1, 3, figsize=(16, 6), constrained_layout=True)
    fig.patch.set_facecolor("white")
    for ax in axes:
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    # Geographic coordinates preserve the spacing of the telescoping buffers.
    im0 = axes[0].pcolormesh(lon, lat, etan0, shading="auto", cmap=cmo.cm.curl, norm=TwoSlopeNorm(vcenter=0, vmin=-etan_lim, vmax=etan_lim))
    im1 = axes[1].pcolormesh(lon, lat, w0, shading="auto", cmap=cmo.cm.balance, norm=TwoSlopeNorm(vcenter=0, vmin=-w_lim, vmax=w_lim))
    im2 = axes[2].pcolormesh(lon, lat, speed(u0, v0), shading="auto", cmap=cmo.cm.speed, vmin=0, vmax=spd_lim)
    for ax in axes:
        ax.set_aspect("equal")

    qstep = args.quiver_step
    q = axes[2].quiver(
        lon[::qstep, ::qstep],
        lat[::qstep, ::qstep],
        u0[::qstep, ::qstep],
        v0[::qstep, ::qstep],
        color="k",
        scale=args.quiver_scale if args.quiver_scale > 0 else 20 * spd_lim,
        width=0.0022,
        alpha=0.75,
    )
    fig.colorbar(im0, ax=axes[0], shrink=0.82).set_label("ETAN (m)")
    fig.colorbar(im1, ax=axes[1], shrink=0.82).set_label("WVEL (m s$^{-1}$)")
    fig.colorbar(im2, ax=axes[2], shrink=0.82).set_label("Speed (m s$^{-1}$)")
    axes[0].set_title("SSH / ETAN")
    axes[1].set_title(f"Vertical velocity ({depth_m:.1f} m)")
    axes[2].set_title(f"Horizontal speed + velocity ({depth_m:.1f} m)")
    title = fig.suptitle("")

    writer = FFMpegWriter(fps=args.fps, metadata={"title": f"{run_label} MITgcm QC animation"})
    with writer.saving(fig, str(args.output), dpi=args.dpi):
        for i, iteration in enumerate(iters):
            etan = masked(read_surf(args.root, iteration, ny, nx, yslice, xslice), surface_mask)
            u, v, w = read_dyn_level(args.root, iteration, level_in_file, ny, nx, yslice, xslice)
            u = masked(u, level_mask)
            v = masked(v, level_mask)
            w = masked(w, level_mask)
            im0.set_array(etan.ravel())
            im1.set_array(w.ravel())
            im2.set_array(speed(u, v).ravel())
            q.set_UVC(u[::qstep, ::qstep], v[::qstep, ::qstep])
            title.set_text(f"{run_label} QC, day {iteration * args.delta_t / 86400:.2f} ({iteration * args.delta_t / 3600:.1f} h)")
            writer.grab_frame()
            if (i + 1) % 24 == 0 or i == len(iters) - 1:
                print(f"rendered {i + 1}/{len(iters)} frames", flush=True)
    plt.close(fig)
    print(f"wrote {args.output}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--label", default="")
    parser.add_argument("--depth-m", type=float, default=500.0)
    parser.add_argument("--days", type=float, default=15.0)
    parser.add_argument("--delta-t", type=float, default=60.0)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--dpi", type=int, default=120)
    parser.add_argument("--quiver-step", type=int, default=20)
    parser.add_argument("--quiver-scale", type=float, default=4.0, help="0 scales arrows to the fixed velocity color limit")
    parser.add_argument("--max-frames", type=int, default=0)
    parser.add_argument("--x-start", type=int, default=0)
    parser.add_argument("--x-end", type=int, default=0, help="exclusive; 0 means full width")
    parser.add_argument("--y-start", type=int, default=0)
    parser.add_argument("--y-end", type=int, default=0, help="exclusive; 0 means full height")
    return parser.parse_args()


if __name__ == "__main__":
    make_animation(parse_args())
