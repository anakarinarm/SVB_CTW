#!/usr/bin/env python3
"""Make a quick-look Rplus-Rminus animation from MITgcm MDS diagnostics."""

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


def completed_iters(root: Path, stem: str, expected_bytes: int) -> set[int]:
    out = set()
    for meta in root.glob(f"{stem}.*.meta"):
        match = ITER_RE.search(meta.name)
        if not match:
            continue
        data = meta.with_suffix(".data")
        if data.exists() and data.stat().st_size == expected_bytes:
            out.add(int(match.group(1)))
    return out


def read_surf(root: Path, iteration: int, ny: int, nx: int) -> np.ndarray:
    mm = np.memmap(root / f"surf_avg.{iteration:010d}.data", dtype=">f8", mode="r", shape=(2, ny, nx))
    return np.asarray(mm[0]).copy()


def read_dyn_level(root: Path, iteration: int, level_in_file: int, ny: int, nx: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mm = np.memmap(root / f"dyn3d_avg.{iteration:010d}.data", dtype=">f8", mode="r", shape=(5, len(SAVED_LEVELS), ny, nx))
    return (
        np.asarray(mm[0, level_in_file]).copy(),
        np.asarray(mm[1, level_in_file]).copy(),
        np.asarray(mm[2, level_in_file]).copy(),
    )


def masked(a: np.ndarray, mask: np.ndarray) -> np.ma.MaskedArray:
    return np.ma.masked_where(~mask | ~np.isfinite(a), a)


def robust_limits(
    plus: Path,
    minus: Path,
    iters: list[int],
    surface_mask: np.ndarray,
    level_mask: np.ndarray,
    level_in_file: int,
    ny: int,
    nx: int,
    yslice: slice,
    xslice: slice,
) -> tuple[float, float, float]:
    etan_vals = []
    w_vals = []
    spd_vals = []
    for iteration in iters[:: max(1, len(iters) // 24)]:
        detan = (read_surf(plus, iteration, ny, nx) - read_surf(minus, iteration, ny, nx))[yslice, xslice]
        up, vp, wp = read_dyn_level(plus, iteration, level_in_file, ny, nx)
        um, vm, wm = read_dyn_level(minus, iteration, level_in_file, ny, nx)
        du = masked((up - um)[yslice, xslice], level_mask)
        dv = masked((vp - vm)[yslice, xslice], level_mask)
        dw = masked((wp - wm)[yslice, xslice], level_mask)
        etan_vals.append(np.asarray(masked(detan, surface_mask).compressed()))
        w_vals.append(np.asarray(dw.compressed()))
        spd_vals.append(np.asarray(np.ma.sqrt(du * du + dv * dv).compressed()))
    etan_lim = float(np.nanpercentile(np.abs(np.concatenate(etan_vals)), 99.5))
    w_lim = float(np.nanpercentile(np.abs(np.concatenate(w_vals)), 99.5))
    spd_lim = float(np.nanpercentile(np.concatenate(spd_vals), 99.5))
    return max(etan_lim, 1e-8), max(w_lim, 1e-10), max(spd_lim, 1e-8)


def make_animation(args: argparse.Namespace) -> None:
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plus_label = args.plus_label or args.plus.name
    minus_label = args.minus_label or args.minus.name
    diff_label = f"{plus_label} - {minus_label}"
    lon = mds.rdmds(str(args.plus / "XC"))
    lat = mds.rdmds(str(args.plus / "YC"))
    hfac_plus = mds.rdmds(str(args.plus / "hFacC"))
    hfac_minus = mds.rdmds(str(args.minus / "hFacC"))
    ny, nx = lon.shape
    xslice = slice(args.x_start, args.x_end if args.x_end > 0 else nx)
    yslice = slice(args.y_start, args.y_end if args.y_end > 0 else ny)
    lon_crop = lon[yslice, xslice]
    lat_crop = lat[yslice, xslice]

    level_in_file = int(np.argmin(np.abs(ZC[SAVED_LEVELS - 1] - args.depth_m)))
    level_index = int(SAVED_LEVELS[level_in_file] - 1)
    depth_m = float(ZC[level_index])
    surface_mask = (hfac_plus[0] > 0)[yslice, xslice]
    level_mask = ((hfac_plus[level_index] > 0) & (hfac_minus[level_index] > 0))[yslice, xslice]

    surf_bytes = 2 * ny * nx * 8
    dyn_bytes = 5 * len(SAVED_LEVELS) * ny * nx * 8
    iters = sorted(
        completed_iters(args.plus, "surf_avg", surf_bytes)
        & completed_iters(args.plus, "dyn3d_avg", dyn_bytes)
        & completed_iters(args.minus, "surf_avg", surf_bytes)
        & completed_iters(args.minus, "dyn3d_avg", dyn_bytes)
    )
    iters = [it for it in iters if it * args.delta_t <= args.days * 86400]
    iters = [it for it in iters if (it * args.delta_t) % 3600 == 0]
    if args.max_frames and len(iters) > args.max_frames:
        iters = iters[: args.max_frames]
    if not iters:
        raise RuntimeError("no complete common plus/minus diagnostics found")

    print(f"using {len(iters)} hourly frames from iteration {iters[0]} to {iters[-1]} ({iters[-1] * args.delta_t / 86400:.2f} days)", flush=True)
    print(f"requested depth {args.depth_m:g} m; using model level {level_index + 1} at {depth_m:.1f} m", flush=True)
    etan_lim, w_lim, spd_lim = robust_limits(args.plus, args.minus, iters, surface_mask, level_mask, level_in_file, ny, nx, yslice, xslice)
    print(f"color limits: dETAN +/-{etan_lim:.4g} m, dWVEL +/-{w_lim:.4g} m/s, d-speed 0-{spd_lim:.4g} m/s", flush=True)

    extent = [float(np.nanmin(lon_crop)), float(np.nanmax(lon_crop)), float(np.nanmin(lat_crop)), float(np.nanmax(lat_crop))]
    fig, axes = plt.subplots(1, 3, figsize=(16, 6), constrained_layout=True)
    fig.patch.set_facecolor("white")
    for ax in axes:
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    def frame(iteration: int) -> tuple[np.ma.MaskedArray, np.ma.MaskedArray, np.ma.MaskedArray, np.ma.MaskedArray, np.ma.MaskedArray]:
        detan = masked((read_surf(args.plus, iteration, ny, nx) - read_surf(args.minus, iteration, ny, nx))[yslice, xslice], surface_mask)
        up, vp, wp = read_dyn_level(args.plus, iteration, level_in_file, ny, nx)
        um, vm, wm = read_dyn_level(args.minus, iteration, level_in_file, ny, nx)
        du = masked((up - um)[yslice, xslice], level_mask)
        dv = masked((vp - vm)[yslice, xslice], level_mask)
        dw = masked((wp - wm)[yslice, xslice], level_mask)
        dspd = np.ma.sqrt(du * du + dv * dv)
        return detan, dw, dspd, du, dv

    detan0, dw0, dspd0, du0, dv0 = frame(iters[0])
    im0 = axes[0].imshow(detan0, origin="lower", extent=extent, cmap=cmo.cm.curl, norm=TwoSlopeNorm(vcenter=0, vmin=-etan_lim, vmax=etan_lim))
    im1 = axes[1].imshow(dw0, origin="lower", extent=extent, cmap=cmo.cm.balance, norm=TwoSlopeNorm(vcenter=0, vmin=-w_lim, vmax=w_lim))
    im2 = axes[2].imshow(dspd0, origin="lower", extent=extent, cmap=cmo.cm.speed, vmin=0, vmax=spd_lim)
    qstep = args.quiver_step
    quiver_scale = args.quiver_scale if args.quiver_scale > 0 else 20 * spd_lim
    q = axes[2].quiver(lon_crop[::qstep, ::qstep], lat_crop[::qstep, ::qstep], du0[::qstep, ::qstep], dv0[::qstep, ::qstep], color="k", scale=quiver_scale, width=0.0022, alpha=0.75)
    fig.colorbar(im0, ax=axes[0], shrink=0.82).set_label("SSH difference (m)")
    fig.colorbar(im1, ax=axes[1], shrink=0.82).set_label("WVEL difference (m s$^{-1}$)")
    fig.colorbar(im2, ax=axes[2], shrink=0.82).set_label("Horizontal velocity difference magnitude (m s$^{-1}$)")
    axes[0].set_title("SSH difference")
    axes[1].set_title(f"Vertical velocity difference ({depth_m:.1f} m)")
    axes[2].set_title(f"Horizontal velocity difference ({depth_m:.1f} m)")
    title = fig.suptitle("")

    writer = FFMpegWriter(fps=args.fps, metadata={"title": f"{diff_label} MITgcm QC animation"})
    with writer.saving(fig, str(args.output), dpi=args.dpi):
        for i, iteration in enumerate(iters):
            detan, dw, dspd, du, dv = frame(iteration)
            im0.set_data(detan)
            im1.set_data(dw)
            im2.set_data(dspd)
            q.set_UVC(du[::qstep, ::qstep], dv[::qstep, ::qstep])
            title.set_text(f"{diff_label} QC, day {iteration * args.delta_t / 86400:.2f} ({iteration * args.delta_t / 3600:.1f} h)")
            writer.grab_frame()
            if (i + 1) % 24 == 0 or i == len(iters) - 1:
                print(f"rendered {i + 1}/{len(iters)} frames", flush=True)
    plt.close(fig)
    print(f"wrote {args.output}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plus", type=Path, required=True)
    parser.add_argument("--minus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plus-label", default="")
    parser.add_argument("--minus-label", default="")
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
