#!/usr/bin/env python3
"""Write figures and a LaTeX summary for the final CTW/SVB rerun setup."""

from __future__ import annotations

from pathlib import Path

import cmocean as cmo
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from MITgcmutils import mds


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "ctw_reruns_telescoped_west64_south72"
INPUT = ROOT / "input" / "telescoped_west64_south72"
LACIE = Path("/Volumes/LaCie/CTW_SVB_reruns_telescoped_west64_south72")
OUT = RUN_ROOT / "model_configuration_summary"
FIG = OUT / "figures"

NX = 576
NY = 684
NR = 100
NFORCING = 768
DELTA_T = 60.0
DIAG_FREQ = 1200.0
END_TIME = 1296000.0
FORCING_PERIOD = 1800.0
FORCING_CYCLE = 1382400.0
RHO_NIL = 999.8
T_ALPHA = 2.0e-4
S_BETA = 7.4e-4
T_REF = 0.0
S_REF = 35.0

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
ZC = np.cumsum(DELR) - 0.5 * DELR
SAVED_LEVELS = list(range(1, 66)) + [67, 69, 71, 73, 75, 77, 79, 84, 90, 100]


def read_bin(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    return np.memmap(path, dtype=">f8", mode="r", shape=shape)


def initial_profiles() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    theta = read_bin(INPUT / "feb_Tref_100zlev_576x684_W64S72.bin", (NR, NY, NX))
    salt = read_bin(INPUT / "feb_Sref_100zlev_576x684_W64S72.bin", (NR, NY, NX))
    tprof = np.asarray(theta[:, NY // 2, NX // 2]).copy()
    sprof = np.asarray(salt[:, NY // 2, NX // 2]).copy()
    rho = RHO_NIL * (1.0 - T_ALPHA * (tprof - T_REF) + S_BETA * (sprof - S_REF))
    return ZC, tprof, sprof, rho


def wind_fields() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    lon = mds.rdmds(str(LACIE / "Rplus" / "XC"))
    lat = mds.rdmds(str(LACIE / "Rplus" / "YC"))
    tau_x = read_bin(INPUT / "upwelling_zonal_halfgauss2hr_16day_576x684_W64S72.bin", (NFORCING, NY, NX))
    tau_y = read_bin(INPUT / "upwelling_merid_halfgauss2hr_16day_576x684_W64S72.bin", (NFORCING, NY, NX))

    max_mag = np.zeros((NY, NX), dtype=np.float64)
    domain_max_by_time = np.zeros(NFORCING, dtype=np.float64)
    imax = 0
    best = -np.inf
    for it in range(NFORCING):
        mag = np.hypot(tau_x[it], tau_y[it])
        max_mag = np.maximum(max_mag, mag)
        domain_max_by_time[it] = float(np.nanmax(mag))
        if domain_max_by_time[it] > best:
            best = domain_max_by_time[it]
            imax = it
    return lon, lat, max_mag, np.asarray(tau_x[imax]), np.asarray(tau_y[imax]), imax


def write_profile_figure() -> None:
    z, theta, salt, rho = initial_profiles()
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 5.5), sharey=True, constrained_layout=True)
    axes[0].plot(theta, -z, color="#b23a48", lw=2)
    axes[1].plot(salt, -z, color="#1f78b4", lw=2)
    axes[2].plot(rho, -z, color="#333333", lw=2)
    axes[0].set_xlabel(r"Potential temperature ($^\circ$C)")
    axes[1].set_xlabel("Salinity")
    axes[2].set_xlabel(r"Linear EOS density (kg m$^{-3}$)")
    axes[0].set_ylabel("Depth (m)")
    for ax in axes:
        ax.grid(True, color="0.86", lw=0.8)
        ax.set_ylim(-4500, 0)
    fig.savefig(FIG / "initial_profiles.pdf")
    fig.savefig(FIG / "initial_profiles.png", dpi=180)
    plt.close(fig)


def bathymetry_fields() -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    lon = mds.rdmds(str(LACIE / "Rplus" / "XC"))
    lat = mds.rdmds(str(LACIE / "Rplus" / "YC"))
    files = {
        "Rplus": "bathy_realistic_plusSVB_W64S72.bin",
        "Rminus": "bathy_realistic_minusSVB_W64S72.bin",
        "Splus": "bathy_smooth30_plusSVB_W64S72.bin",
        "Sminus": "bathy_smooth30_minusSVB_W64S72.bin",
    }
    fields = {}
    for label, fname in files.items():
        raw = read_bin(INPUT / fname, (NY, NX))
        fields[label] = np.where(raw < 0.0, -np.asarray(raw), np.nan)
    return lon, lat, fields


def write_bathymetry_figure() -> None:
    lon, lat, fields = bathymetry_fields()
    cmap_bathy = cmo.cm.deep.copy()
    cmap_bathy.set_bad("0.82")
    cmap_diff = cmo.cm.balance.copy()
    cmap_diff.set_bad("0.82")

    fig, axes = plt.subplots(3, 2, figsize=(11, 13), sharex=True, sharey=True, constrained_layout=True)
    panels = [
        ("Rplus: realistic +SVB", fields["Rplus"], cmap_bathy, 0.0, 4700.0, "Depth (m)"),
        ("Rminus: realistic -SVB", fields["Rminus"], cmap_bathy, 0.0, 4700.0, "Depth (m)"),
        ("Splus: smooth30 +SVB", fields["Splus"], cmap_bathy, 0.0, 4700.0, "Depth (m)"),
        ("Sminus: smooth30 -SVB", fields["Sminus"], cmap_bathy, 0.0, 4700.0, "Depth (m)"),
        ("Rplus - Splus", fields["Rplus"] - fields["Splus"], cmap_diff, -500.0, 500.0, "Depth difference (m)"),
        ("Rplus - Rminus", fields["Rplus"] - fields["Rminus"], cmap_diff, -250.0, 250.0, "Depth difference (m)"),
    ]
    for ax, (title, data, cmap, vmin, vmax, label) in zip(axes.ravel(), panels):
        ax.set_facecolor("0.82")
        im = ax.pcolormesh(lon, lat, np.ma.masked_invalid(data), cmap=cmap, shading="auto", vmin=vmin, vmax=vmax)
        ax.contour(lon, lat, np.isfinite(fields["Rplus"]).astype(float), levels=[0.5], colors="k", linewidths=0.45)
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        fig.colorbar(im, ax=ax, shrink=0.82).set_label(label)
    fig.savefig(FIG / "bathymetry_panels.pdf")
    fig.savefig(FIG / "bathymetry_panels.png", dpi=180)
    plt.close(fig)


def write_wind_figure() -> tuple[float, float]:
    lon, lat, max_mag, tx, ty, imax = wind_fields()
    qstep = 24
    fig, ax = plt.subplots(figsize=(8.2, 7.4), constrained_layout=True)
    im = ax.pcolormesh(lon, lat, max_mag, shading="auto", cmap=cmo.cm.speed)
    ax.quiver(
        lon[::qstep, ::qstep],
        lat[::qstep, ::qstep],
        tx[::qstep, ::qstep],
        ty[::qstep, ::qstep],
        color="k",
        scale=4.0,
        width=0.0021,
        alpha=0.75,
    )
    fig.colorbar(im, ax=ax, shrink=0.86).set_label(r"max $|\tau|$ (N m$^{-2}$)")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Maximum wind-stress magnitude; arrows at domain maximum")
    fig.savefig(FIG / "windstress_max_map.pdf")
    fig.savefig(FIG / "windstress_max_map.png", dpi=180)
    plt.close(fig)
    return float(np.nanmax(max_mag)), imax * FORCING_PERIOD / 3600.0


def write_tex(wind_max: float, wind_hour: float) -> None:
    tex = rf"""\documentclass[11pt]{{article}}
\usepackage[margin=1in]{{geometry}}
\usepackage{{graphicx}}
\usepackage{{booktabs}}
\usepackage{{siunitx}}
\usepackage{{amsmath}}
\usepackage{{hyperref}}

\title{{MITgcm configuration summary for the CTW--SVB reruns}}
\author{{}}
\date{{}}

\begin{{document}}
\maketitle

\section*{{Model domain and grid}}
The final rerun suite uses MITgcm on a spherical-polar C-grid with
$N_x \times N_y \times N_r = 576 \times 684 \times 100$ grid cells.  The grid
origin supplied to MITgcm is approximately
\[
(x_g,y_g)=(236.36437^\circ\mathrm{{E}},25.39553^\circ\mathrm{{N}}).
\]
This is the original $512 \times 612$ active grid augmented by 64 telescoped
columns to the west and 72 telescoped rows to the south.  The added wet buffer
is a flat \SI{{4400}}{{m}} abyssal plane with straight artificial coastlines.
Horizontal spacing ranges from $0.0155945^\circ$--$0.0389862^\circ$ in
longitude and $0.0135994^\circ$--$0.0339985^\circ$ in latitude.  The vertical
grid has 100 levels and reaches \SI{{4700}}{{m}} depth.

\section*{{Bathymetry}}
The experiment suite uses four bathymetry files that differ only in the
resolved topography and presence or absence of Sebasti\'an Vizca\'ino Bay (SVB):
\texttt{{Rplus}} uses realistic bathymetry with SVB, \texttt{{Rminus}} uses
realistic bathymetry with SVB removed, \texttt{{Splus}} uses smoothed bathymetry
with SVB, and \texttt{{Sminus}} uses smoothed bathymetry with SVB removed.  The
MITgcm bathymetry files are big-endian binary arrays with negative ocean depth
and zero over land.

The realistic +SVB bathymetry is the cleaned baseline bathymetry.  The -SVB
members modify only the shelf/bay source geometry: the continental slope is
preserved seaward of the shelf-break guide, and the 250-m shelf-break spike near
$28.6^\circ$N is regularized to avoid an artificial land spike.  The smoothed
members use the approved 30-km half-power native-grid local-tangent smoothing,
not the historical repeated smoother.  Realistic and smoothed members retain the
same wet/dry mask within each +SVB or -SVB pair, and the immediate coastal
depths are held fixed so that the comparison isolates short-scale bathymetric
roughness rather than coastline or partial-cell changes.

For the west64/south72 reruns, these four finalized bathymetries were inserted
into the expanded domain.  The added western and southern buffer regions are
flat \SI{{4400}}{{m}} abyssal ocean at the open boundaries, with straight artificial
coastlines.  The western added buffer is smoothly blended toward the local
western-edge depths, and the southern added buffer is smoothly blended toward
the first active-domain row, so the southwest, southeast, and northwest sponge
corners are continuous.  In the northern buffer/taper band, existing wet cells
continue the local depth from the row immediately south of the band while
preserving the land mask and the configured north-boundary wet span.  The
southwest added-buffer corner is additionally smoothed with a local 2-D
wet-cell relaxation, and the northwest added-buffer corner uses a smoothed
along-boundary profile while retaining a gradual northward transition.  The
open-boundary bathymetry treatment in the finalized source bathymetries uses a
15-cell flat normal-depth core with a 15-cell smooth transition on the west,
south, and north boundaries.

\begin{{figure}}[h]
\centering
\includegraphics[width=0.95\linewidth]{{figures/bathymetry_panels.pdf}}
\caption{{Production bathymetries used by the four west64/south72 reruns and
two diagnostic differences.  Depths are shown positive downward; land is gray.}}
\end{{figure}}

\section*{{Time stepping and output}}
The simulations start from rest at $t=0$ and run for \SI{{1296000}}{{s}}
(15 days), with $\Delta t=\SI{{60}}{{s}}$.  Checkpoints are written every
\SI{{43200}}{{s}} and model dumps every \SI{{172800}}{{s}}.  Diagnostics are
written in MITgcm MDS binary format rather than MNC NetCDF.

\begin{{table}}[h]
\centering
\begin{{tabular}}{{lll}}
\toprule
File prefix & Frequency & Variables \\
\midrule
\texttt{{dyn3d\_avg}} & \SI{{1200}}{{s}} averages &
\texttt{{UVEL}}, \texttt{{VVEL}}, \texttt{{WVEL}}, \texttt{{RHOAnoma}}, \texttt{{PHIHYD}} \\
\texttt{{surf\_avg}} & \SI{{1200}}{{s}} averages &
\texttt{{ETAN}}, \texttt{{PHI\_SURF}} \\
\texttt{{surf\_snap}} & \SI{{1200}}{{s}} snapshots &
\texttt{{ETAN}}, \texttt{{PHI\_SURF}} \\
\bottomrule
\end{{tabular}}
\end{{table}}

The 3-D diagnostics are saved on 75 vertical levels:
\texttt{{1--65, 67, 69, 71, 73, 75, 77, 79, 84, 90, 100}}.

\section*{{Equation of state, initial conditions, and open boundaries}}
The runs use a linear equation of state with $\rho_0=\SI{{999.8}}{{kg\,m^{{-3}}}}$,
$\alpha=2.0\times10^{{-4}}~^\circ\mathrm{{C}}^{{-1}}$, and
$\beta=7.4\times10^{{-4}}$.  MITgcm reference values are set to
$T_\mathrm{{ref}}=0^\circ\mathrm{{C}}$ and $S_\mathrm{{ref}}=35$, so that the
linear equation of state is satisfied by constant reference values.  The
initial temperature and salinity fields are vertically varying and horizontally
uniform; the same profiles are prescribed at the open boundaries using OBCS
temperature and salinity files.  No open-boundary sponge relaxation is applied
to temperature or salinity.

\begin{{figure}}[h]
\centering
\includegraphics[width=0.95\linewidth]{{figures/initial_profiles.pdf}}
\caption{{Initial vertically varying, horizontally uniform temperature,
salinity, and linear-equation-of-state density profiles.}}
\end{{figure}}

\section*{{Open boundaries and sponge}}
Open boundaries are active on the west, south, and north boundaries.  The OBCS
masks are
\texttt{{OB\_Iwest = 684*1}},
\texttt{{OB\_Jsouth = 543*1,33*0}}, and
\texttt{{OB\_Jnorth = 137*-1,439*0}}.  Orlanski radiation is enabled on all
three open boundaries, with \texttt{{cvelTimeScale=2000.0}} and
\texttt{{CMAX=1.0}}.  OBCS volume balance is enabled with unit balance factors
on west, south, and north boundaries.

A linear OBCS sponge is active for momentum only.  The sponge thickness is
48 grid cells, with inner relaxation time scale \SI{{43200}}{{s}} and boundary
relaxation time scale \SI{{1800}}{{s}} for both $u$ and $v$.  The sponge acts
on west, south, and north boundaries and relaxes both east--west and north--south
velocity components.  Temperature and salinity sponge switches are disabled.

\section*{{Forcing and experiment set}}
Wind stress is periodic with a \SI{{1800}}{{s}} forcing period and a
\SI{{1382400}}{{s}} forcing cycle.  The forcing files prescribe zonal and
meridional wind stress over the $576 \times 684$ domain.  The maximum wind
stress magnitude over the full forcing cycle is
\SI{{{wind_max:.4g}}}{{N\,m^{{-2}}}}, occurring at hour {wind_hour:.1f} of the
forcing record.  The four experiments are \texttt{{Rplus}}, \texttt{{Rminus}},
\texttt{{Splus}}, and \texttt{{Sminus}}, as described in the bathymetry section
above.

\begin{{figure}}[h]
\centering
\includegraphics[width=0.82\linewidth]{{figures/windstress_max_map.pdf}}
\caption{{Maximum wind-stress magnitude at each grid point over the forcing
cycle.  Arrows show the wind-stress vector at the time of the domain-wide
maximum.}}
\end{{figure}}

\end{{document}}
"""
    (OUT / "model_configuration_summary.tex").write_text(tex)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    write_profile_figure()
    write_bathymetry_figure()
    wind_max, wind_hour = write_wind_figure()
    write_tex(wind_max, wind_hour)
    print(f"wrote {OUT / 'model_configuration_summary.tex'}")
    print(f"wrote {FIG / 'initial_profiles.pdf'}")
    print(f"wrote {FIG / 'bathymetry_panels.pdf'}")
    print(f"wrote {FIG / 'windstress_max_map.pdf'}")


if __name__ == "__main__":
    main()
