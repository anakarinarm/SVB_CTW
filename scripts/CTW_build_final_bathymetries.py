#!/usr/bin/env python3
"""Build the final CTW bathymetry suite from GEBCO + the MITgcm horizontal grid.

This script reproduces the finalized workflow developed for the CTW/SVB reruns:

A. COMMON numerical/topological cleanup
   1. Direct GEBCO -> actual MITgcm XC/YC interpolation.
   2. Intentional Gulf of California / Ojo de Liebre land edits from the old setup.
   3. Clip depths at 4700 m and keep the eastern model wall dry.
   4. Conservative deep-ocean pass-line cleanup (common to all experiments).
   5. Remove disconnected wet puddles, tiny/near-mainland dry islands, and
      one-cell wet dead ends.

B. REALISTIC +SVB / -SVB source geometry
   1. +SVB is the cleaned common baseline.
   2. -SVB keeps the continental slope exactly unchanged seaward of the shelf-break
      guide and replaces only the shelf with the profile-A shelf shape.
   3. The awkward 28.6 N shelf-break spike is deliberately removed by a linear
      bridge between its two flanks (28.53--28.6932 N).

C. SMOOTH experiment
   A seam-free native-grid Gaussian filter is applied only along the local tangent
   to broad bathymetric contours. The source region, immediate coast, deep ocean,
   and region south of profile A are protected.

D. OPEN BOUNDARIES
   1. Southern coastline is made exactly straight/perpendicular to the OBC over
      rows 0--30, using the natural coast at rows 30/31.
   2. West, south, and north OBC bathymetry: 15-cell flat core + 15-cell smooth
      transition. East is closed.

E. PARTIAL-CELL PREPARATION
   1. Common minimum positive depth Hmin = 2 m.
   2. Exact production delR (100 levels, total 4700 m).
   3. Diagnostics for hFacMin=0.10, hFacMinDr=2 m,
      useMin4hFacEdges=.TRUE., pCellMix_select=0.

Outputs
-------
- Final NetCDF containing final bathymetries and useful intermediate fields.
- Four MITgcm big-endian float64 bathymetry binaries (negative ocean depth, land=0).
- Exact 2-D partial-cell/effective-depth diagnostic NetCDF.
- CSV audit tables and a small data namelist snippet.

Example
-------
python CTW_build_final_bathymetries.py \
  --gebco gebco_2020_n36.0_s26.0_w-122.0_e-113.5.nc \
  --grid SVB_10_grid_reduced.nc \
  --outdir CTW_final_bathymetry_package
"""

from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import xarray as xr
import h5py
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import (
    distance_transform_edt,
    gaussian_filter,
    label,
    map_coordinates,
    median_filter,
)

# -------------------------- fixed experiment choices -------------------------
MAX_DEPTH_M = 4700.0
HMIN_FINAL_M = 2.0

# Common deep GEBCO acquisition-line cleanup.
DEEP_CLEAN_START_M = 2500.0
DEEP_CLEAN_FULL_M = 3000.0
DEEP_LAND_PROTECT_KM = 10.0
DEEP_LAND_FULL_KM = 20.0
DEEP_MEDIAN_SIZE = 5
DEEP_CORRECTION_CAP_M = 250.0

# Common topology.
MAX_WET_DEAD_END_ITER = 20

# no-SVB source construction.
LAT_A = 29.81
SHELF_BREAK_M = 250.0
DIAG_SMOOTH_KM = 5.0
SOUTH_JOIN_MIN = 27.55
SOUTH_JOIN_MAX = 27.90
SOUTH_DEPTH_RAMP_KM = 10.0
SPIKE_BRIDGE_SOUTH = 28.530
SPIKE_BRIDGE_NORTH = 28.6932

# Smooth bathymetry.
ORIENTATION_SIGMA_KM = 10.0
LAMBDA50_POWER_KM = 30.0
FILTER_SAMPLE_STEP_KM = 1.0
COAST_EXACT_KM = 3.0
COAST_FULL_FILTER_KM = 9.0
SOURCE_RAMP_KM = 50.0
FILTER_FULL_DEPTH_M = 2000.0
FILTER_ZERO_DEPTH_M = 3000.0
# Keep the previously accepted broad source-protection band exact.
SOURCE_PROTECT_OFFSHORE_KM = -30.0

# OBC.
N_OBC_CORE = 15
N_OBC_BLEND = 15
N_OBC_TOTAL = N_OBC_CORE + N_OBC_BLEND

# Partial cells.
HFAC_MIN = 0.10
HFAC_MIN_DR_M = 2.0
PCELL_MIX_SELECT = 0
DEL_R = np.array([
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
], dtype=float)


# -------------------------------- utilities ---------------------------------
def smoothstep(x: np.ndarray | float) -> np.ndarray:
    q = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    return q*q*(3.0 - 2.0*q)


def model_dx_dy_km(lon: np.ndarray, lat: np.ndarray, lat_ref: float = 30.0) -> tuple[float, float]:
    dx = abs(float(np.median(np.diff(lon)))) * 111.32 * np.cos(np.deg2rad(lat_ref))
    dy = abs(float(np.median(np.diff(lat)))) * 111.32
    return dx, dy


def wet_neighbor_count(mask: np.ndarray) -> np.ndarray:
    n = np.zeros(mask.shape, dtype=np.int16)
    n[1:, :] += mask[:-1, :]
    n[:-1, :] += mask[1:, :]
    n[:, 1:] += mask[:, :-1]
    n[:, :-1] += mask[:, 1:]
    return n


def remove_disconnected_wet_puddles(H: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, pd.DataFrame]:
    """Keep the largest 8-connected wet ocean component; dry the rest."""
    wet = H > 0
    labs, nlab = label(wet, structure=np.ones((3, 3), dtype=int))
    if nlab <= 1:
        return H.copy(), pd.DataFrame()
    sizes = np.bincount(labs.ravel())
    main = 1 + int(np.argmax(sizes[1:]))
    out = H.copy()
    rows: list[dict] = []
    for lab_id in range(1, nlab+1):
        if lab_id == main:
            continue
        jj, ii = np.where(labs == lab_id)
        if len(jj) == 0:
            continue
        rows.append({
            "component": lab_id,
            "n_cells": int(len(jj)),
            "mean_lat": float(np.mean(lat[jj])),
            "mean_lon": float(np.mean(lon[ii])),
            "min_depth_m": float(np.min(H[jj, ii])),
            "max_depth_m": float(np.max(H[jj, ii])),
        })
        out[jj, ii] = 0.0
    return out, pd.DataFrame(rows)


def dry_components_4(mask_wet: np.ndarray):
    dry = ~mask_wet
    structure = np.array([[0,1,0],[1,1,1],[0,1,0]], dtype=int)
    labs, nlab = label(dry, structure=structure)
    comps = []
    ny, nx = dry.shape
    for lab_id in range(1, nlab+1):
        jj, ii = np.where(labs == lab_id)
        if len(jj) == 0:
            continue
        touches_edge = np.any((jj == 0) | (jj == ny-1) | (ii == 0) | (ii == nx-1))
        comps.append((lab_id, jj, ii, touches_edge))
    return labs, comps


def local_plane_fill(H: np.ndarray, cells: list[tuple[int, int]], lon: np.ndarray, lat: np.ndarray, radius: int = 2):
    out = H.copy()
    jj0 = np.array([c[0] for c in cells], dtype=int)
    ii0 = np.array([c[1] for c in cells], dtype=int)
    jmin, jmax = max(0, jj0.min()-radius), min(H.shape[0]-1, jj0.max()+radius)
    imin, imax = max(0, ii0.min()-radius), min(H.shape[1]-1, ii0.max()+radius)

    latc = float(np.mean(lat[jj0]))
    lonc = float(np.mean(lon[ii0]))
    xscale = 111.32*np.cos(np.deg2rad(latc))
    yscale = 111.32
    cellset = set(cells)

    pts, vals, weights = [], [], []
    for j in range(jmin, jmax+1):
        for i in range(imin, imax+1):
            if (j, i) in cellset or H[j, i] <= 0:
                continue
            x = (lon[i]-lonc)*xscale
            y = (lat[j]-latc)*yscale
            d = np.hypot(x, y)
            pts.append([1.0, x, y])
            vals.append(H[j, i])
            weights.append(1.0/(d*d + 0.25))

    pts = np.asarray(pts, dtype=float)
    vals = np.asarray(vals, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if len(vals) == 0:
        raise RuntimeError("No wet neighbors available to fill a tiny dry island.")

    if len(vals) < 4:
        fills = [float(np.median(vals))]*len(cells)
    else:
        sw = np.sqrt(weights)
        beta, *_ = np.linalg.lstsq(pts*sw[:, None], vals*sw, rcond=None)
        vmin, vmax = float(vals.min()), float(vals.max())
        fills = []
        for j, i in cells:
            x = (lon[i]-lonc)*xscale
            y = (lat[j]-latc)*yscale
            pred = float(beta @ np.array([1.0, x, y]))
            pred = float(np.clip(pred, vmin, vmax))
            if not np.isfinite(pred) or pred <= 0:
                pred = float(np.median(vals))
            fills.append(pred)

    for (j, i), fill in zip(cells, fills):
        out[j, i] = fill
    return out, fills


def remove_small_dry_islands(H: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, pd.DataFrame]:
    """Final rule: open all 1-2 cell islands; open 3-5 cell islands if within
    Chebyshev distance <=2 cells of another land component (mainland or island).
    """
    out = H.copy()
    wet = out > 0
    labs, comps = dry_components_4(wet)
    to_open = []
    for lab_id, jj, ii, touches_edge in comps:
        if touches_edge:
            continue
        size = len(jj)
        cells = list(zip(jj.tolist(), ii.tolist()))
        remove = False
        reason = ""
        if size <= 2:
            remove = True
            reason = "enclosed dry island <= 2 cells"
        elif size <= 5:
            min_cheb = 99
            nearest_lab = None
            for j, i in cells:
                for dj in range(-2, 3):
                    for di in range(-2, 3):
                        if dj == 0 and di == 0:
                            continue
                        jt, it = j+dj, i+di
                        if not (0 <= jt < H.shape[0] and 0 <= it < H.shape[1]):
                            continue
                        other = int(labs[jt, it])
                        if other > 0 and other != lab_id:
                            dist = max(abs(dj), abs(di))
                            if dist < min_cheb:
                                min_cheb = dist
                                nearest_lab = other
            if min_cheb <= 2:
                remove = True
                reason = (
                    f"{size}-cell island within Chebyshev distance {min_cheb} "
                    f"of dry component {nearest_lab}"
                )
        if remove:
            to_open.append((lab_id, cells, reason))

    rows = []
    for lab_id, cells, reason in to_open:
        out, fills = local_plane_fill(out, cells, lon, lat, radius=2)
        for (j, i), fill in zip(cells, fills):
            rows.append({
                "component": int(lab_id),
                "component_size": int(len(cells)),
                "j": int(j), "i": int(i),
                "lat": float(lat[j]), "lon": float(lon[i]),
                "reason": reason,
                "filled_depth_m": float(fill),
            })
    return out, pd.DataFrame(rows)


def remove_wet_dead_ends(H: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, pd.DataFrame]:
    out = H.copy()
    rows = []
    for iteration in range(1, MAX_WET_DEAD_END_ITER+1):
        wet = out > 0
        nwet = wet_neighbor_count(wet)
        bad = wet & (nwet <= 1)
        bad[[0, -1], :] = False
        bad[:, [0, -1]] = False
        if not np.any(bad):
            break
        for j, i in np.argwhere(bad):
            rows.append({
                "iteration": iteration,
                "j": int(j), "i": int(i),
                "lat": float(lat[j]), "lon": float(lon[i]),
                "depth_before_m": float(out[j, i]),
                "wet_cardinal_neighbors": int(nwet[j, i]),
            })
        out[bad] = 0.0
    return out, pd.DataFrame(rows)


def topology_cleanup(H: np.ndarray, lon: np.ndarray, lat: np.ndarray):
    logs = {}
    out, logs["puddles_1"] = remove_disconnected_wet_puddles(H, lon, lat)
    out, logs["dry_islands_1"] = remove_small_dry_islands(out, lon, lat)
    out, logs["wet_dead_ends"] = remove_wet_dead_ends(out, lon, lat)
    out, logs["dry_islands_2"] = remove_small_dry_islands(out, lon, lat)
    out, logs["puddles_2"] = remove_disconnected_wet_puddles(out, lon, lat)
    return out, logs


def topology_audit(mask_wet: np.ndarray) -> dict[str, int]:
    nwet = wet_neighbor_count(mask_wet)
    dangling = mask_wet & (nwet <= 1)
    dangling[[0, -1], :] = False
    dangling[:, [0, -1]] = False

    labs, comps = dry_components_4(mask_wet)
    bad12 = 0
    bad35 = 0
    for lab_id, jj, ii, touches_edge in comps:
        if touches_edge:
            continue
        size = len(jj)
        if size <= 2:
            bad12 += 1
        elif size <= 5:
            min_cheb = 99
            for j, i in zip(jj, ii):
                for dj in range(-2, 3):
                    for di in range(-2, 3):
                        if dj == 0 and di == 0:
                            continue
                        jt, it = j+dj, i+di
                        if (
                            0 <= jt < mask_wet.shape[0] and 0 <= it < mask_wet.shape[1]
                            and labs[jt, it] > 0 and labs[jt, it] != lab_id
                        ):
                            min_cheb = min(min_cheb, max(abs(dj), abs(di)))
            if min_cheb <= 2:
                bad35 += 1
    return {
        "wet_columns": int(mask_wet.sum()),
        "dangling_wet_columns": int(dangling.sum()),
        "enclosed_dry_islands_size_1_2": int(bad12),
        "near_land_dry_islands_size_3_5": int(bad35),
    }


def contour_crossings(row: np.ndarray, lon: np.ndarray, target: float, require_wet: bool = True) -> list[float]:
    out: list[float] = []
    for i in range(len(row)-1):
        a, b = row[i], row[i+1]
        if not (np.isfinite(a) and np.isfinite(b)):
            continue
        if require_wet and not (a > 0 and b > 0):
            continue
        if (a-target)*(b-target) <= 0 and b != a:
            f = (target-a)/(b-a)
            out.append(float(lon[i] + f*(lon[i+1]-lon[i])))
    return out


def persistent_mainland_coast(Hrow: np.ndarray, lon: np.ndarray, start_lon: float, dry_run: int = 6) -> float:
    i0 = int(np.searchsorted(lon, start_lon))
    i0 = int(np.clip(i0, 0, len(lon)-2))
    dry = Hrow <= 0
    for i in range(i0, len(lon)-dry_run+1):
        if dry[i:i+dry_run].all():
            k = i-1
            while k >= i0 and Hrow[k] <= 0:
                k -= 1
            return float(lon[k]) if k >= i0 else float(lon[i])
    return np.nan


# ------------------------------- GEBCO -> H0 --------------------------------
def build_common_H0(gebco_path: Path, grid_path: Path):
    gebco = xr.open_dataset(gebco_path)
    if "elevation" not in gebco:
        raise ValueError("GEBCO file must contain an 'elevation' variable.")
    glon = gebco["lon"].values.astype(float)
    glat = gebco["lat"].values.astype(float)
    elev = gebco["elevation"].values.astype(float).copy()

    # Ocean/land and model-depth cap.
    elev[elev > 0] = 0.0
    elev[elev < -MAX_DEPTH_M] = -MAX_DEPTH_M

    # Exact intentional raw-GEBCO edits inherited from the original workflow.
    lon1, lon2, lat1, lat2 = 0, 1919, 240, 2232
    # Flatten the northeast Gulf of California corner to land.
    elev[lat2-640:lat2+3, lon2-550:lon2+3] = 0.0
    elev[lat2-1110:lat2-640, lon2-300:lon2+3] = 0.0
    elev[lat2-1375:lat2-640, lon2-200:lon2+3] = 0.0
    elev[lat2-760:lat2+3, lon2-430:lon2+3] = 0.0
    # Remove Ojo de Liebre lagoon.
    elev[lat1+145:lat1+180, lon2-80:lon2+3] = 0.0
    elev[lat1+180:lat1+210, lon2-80:lon2+3] = 0.0
    elev[lat1+180:lat1+225, lon2-75:lon2+3] = 0.0

    # The reduced MITgcm grid file is HDF5/NetCDF4. Use h5py directly so the
    # script does not require the optional netCDF4/h5netcdf Python packages.
    with h5py.File(grid_path, "r") as grid:
        lon = np.asarray(grid["XC"][:], dtype=float)
        lat = np.asarray(grid["YC"][:], dtype=float)
    if lon.ndim != 1 or lat.ndim != 1:
        raise ValueError("Grid file must contain 1-D XC and YC coordinates.")
    lon180 = np.where(lon > 180.0, lon-360.0, lon)

    interp = RegularGridInterpolator((glat, glon), elev, bounds_error=False, fill_value=np.nan)
    LO, LA = np.meshgrid(lon180, lat)
    Ec = interp(np.column_stack([LA.ravel(), LO.ravel()])).reshape(LA.shape)
    if np.isnan(Ec).any():
        raise RuntimeError("The GEBCO subset does not fully cover the MITgcm grid.")
    H = np.maximum(-Ec, 0.0)
    H[:, -1] = 0.0  # closed eastern wall

    # Initial disconnected-puddle cleanup.
    H, puddle_log = remove_disconnected_wet_puddles(H, lon180, lat)

    # Conservative deep-ocean pass-line cleanup, common to every experiment.
    dx_km, dy_km = model_dx_dy_km(lon180, lat, 30.0)
    wet = H > 0
    dist_land = distance_transform_edt(wet, sampling=(dy_km, dx_km))
    Hmed = median_filter(H, size=DEEP_MEDIAN_SIZE, mode="nearest")
    Wdepth = smoothstep((H-DEEP_CLEAN_START_M)/(DEEP_CLEAN_FULL_M-DEEP_CLEAN_START_M))
    Wland = smoothstep((dist_land-DEEP_LAND_PROTECT_KM)/(DEEP_LAND_FULL_KM-DEEP_LAND_PROTECT_KM))
    correction = np.clip(Hmed-H, -DEEP_CORRECTION_CAP_M, DEEP_CORRECTION_CAP_M)
    Hart = H + Wdepth*Wland*correction
    Hart[~wet] = 0.0

    # Final common topology cleanup.
    H0, topo_logs = topology_cleanup(Hart, lon180, lat)
    return lon180, lat, H, Hart, H0, puddle_log, topo_logs


# ------------------------- revised shelf-only -SVB ---------------------------
def build_noSVB_shelf_only(Hplus: np.ndarray, lon: np.ndarray, lat: np.ndarray):
    dx_km, dy_km = model_dx_dy_km(lon, lat, 29.0)
    wet = Hplus > 0
    sigma_cells = DIAG_SMOOTH_KM / np.mean([dx_km, dy_km])
    num = gaussian_filter(Hplus*wet, sigma=sigma_cells, mode="nearest")
    den = gaussian_filter(wet.astype(float), sigma=sigma_cells, mode="nearest")
    Hdiag = np.full_like(Hplus, np.nan)
    np.divide(num, den, out=Hdiag, where=den > 0.5)

    jA = int(np.argmin(np.abs(lat-LAT_A)))
    jSouthSearch = int(np.argmin(np.abs(lat-27.50)))

    diagA = [x for x in contour_crossings(Hdiag[jA], lon, SHELF_BREAK_M, require_wet=False)
             if -117.2 < x < -114.8]
    if not diagA:
        raise RuntimeError("Could not find broad 250-m crossing at profile A.")
    prev = min(diagA, key=lambda x: abs(x+115.88))

    diag_guide = {}
    for j in range(jA, jSouthSearch-1, -1):
        cs = [x for x in contour_crossings(Hdiag[j], lon, SHELF_BREAK_M, require_wet=False)
              if -117.2 < x < -114.8]
        if not cs:
            raise RuntimeError(f"No broad 250-m crossing at latitude {lat[j]:.4f}.")
        chosen = min(cs, key=lambda x: abs(x-prev))
        diag_guide[j] = chosen
        prev = chosen

    # Actual 250-m crossing nearest the broad branch.
    actual_sb = {}
    for j in range(jA, jSouthSearch-1, -1):
        cs = [x for x in contour_crossings(Hplus[j], lon, SHELF_BREAK_M, require_wet=True)
              if -117.2 < x < -114.5]
        actual_sb[j] = diag_guide[j] if not cs else min(cs, key=lambda x: abs(x-diag_guide[j]))

    sbA = actual_sb[jA]
    kmdegA = 111.32*np.cos(np.deg2rad(lat[jA]))
    qA = (lon-sbA)*kmdegA
    HA = Hplus[jA].copy()
    coastA = persistent_mainland_coast(HA, lon, sbA)
    qcoastA = (coastA-sbA)*kmdegA

    # Southern join where shifted A coast naturally meets Punta Eugenia coast.
    join_rows = []
    for j, sb in actual_sb.items():
        if not (SOUTH_JOIN_MIN <= lat[j] <= SOUTH_JOIN_MAX):
            continue
        kmdeg = 111.32*np.cos(np.deg2rad(lat[j]))
        predicted = sb + qcoastA/kmdeg
        natural = persistent_mainland_coast(Hplus[j], lon, sb)
        if np.isfinite(natural):
            join_rows.append({
                "j": j, "lat": float(lat[j]), "shelfbreak_lon": float(sb),
                "predicted_shifted_A_coast_lon": float(predicted),
                "natural_coast_lon": float(natural),
                "coast_separation_km": float((predicted-natural)*kmdeg),
            })
    join_df = pd.DataFrame(join_rows)
    if join_df.empty:
        raise RuntimeError("Could not identify the southern no-SVB join.")
    best_join = join_df.iloc[int(np.argmin(np.abs(join_df["coast_separation_km"].values)))]
    jJoin = int(best_join["j"])

    # Smooth the 28.6 N spike in the guide.
    guide_df = pd.DataFrame({
        "j": list(range(jJoin, jA+1)),
        "lat": [float(lat[j]) for j in range(jJoin, jA+1)],
        "actual_250m_shelfbreak_lon": [float(actual_sb[j]) for j in range(jJoin, jA+1)],
    })
    j0 = int(guide_df.loc[(guide_df["lat"]-SPIKE_BRIDGE_SOUTH).abs().idxmin(), "j"])
    j1 = int(guide_df.loc[(guide_df["lat"]-SPIKE_BRIDGE_NORTH).abs().idxmin(), "j"])
    lat0, lat1 = float(lat[j0]), float(lat[j1])
    x0 = float(guide_df.loc[guide_df["j"] == j0, "actual_250m_shelfbreak_lon"].iloc[0])
    x1 = float(guide_df.loc[guide_df["j"] == j1, "actual_250m_shelfbreak_lon"].iloc[0])
    guide_df["smoothed_shelfbreak_lon"] = guide_df["actual_250m_shelfbreak_lon"].copy()
    bridge = (guide_df["j"] >= j0) & (guide_df["j"] <= j1)
    guide_df.loc[bridge, "smoothed_shelfbreak_lon"] = (
        x0 + (x1-x0)*(guide_df.loc[bridge, "lat"]-lat0)/(lat1-lat0)
    )
    guide_df["guide_adjustment_km"] = (
        (guide_df["smoothed_shelfbreak_lon"]-guide_df["actual_250m_shelfbreak_lon"])
        * 111.32*np.cos(np.deg2rad(guide_df["lat"]))
    )

    # Shelf-only replacement.
    Hminus = Hplus.copy()
    noSVB_source_mask = np.zeros_like(Hplus, dtype=bool)
    # Broader source protection used for the accepted smooth +SVB filter.
    filter_source_protect = np.zeros_like(Hplus, dtype=bool)
    ramp_rows = max(1, int(np.ceil(SOUTH_DEPTH_RAMP_KM/dy_km)))
    construction_rows = []

    for _, rec in guide_df.iterrows():
        j = int(rec["j"])
        sb = float(rec["smoothed_shelfbreak_lon"])
        sb_literal = float(rec["actual_250m_shelfbreak_lon"])
        kmdeg = 111.32*np.cos(np.deg2rad(lat[j]))
        q = (lon-sb)*kmdeg
        q_literal = (lon-sb_literal)*kmdeg

        Href = np.interp(q, qA, HA, left=float(HA[0]), right=0.0)
        shelf = q >= 0.0
        noSVB_source_mask[j, shelf] = True
        # Preserve broad source/slope geometry from filtering exactly as in the
        # previously accepted smooth experiment: 30 km offshore of the literal
        # local shelf-break branch through profile A.
        filter_source_protect[j, q_literal >= SOURCE_PROTECT_OFFSHORE_KM] = True

        srcwet = Href > 0
        natwet = Hplus[j] > 0
        walong = 1.0
        if j < jJoin+ramp_rows:
            walong = float(smoothstep((j-jJoin+1)/ramp_rows))

        close = shelf & (~srcwet) & natwet
        Hminus[j, close] = 0.0
        both = shelf & srcwet & natwet
        Hminus[j, both] = (1.0-walong)*Hplus[j, both] + walong*Href[both]

        construction_rows.append({
            "j": j, "lat": float(lat[j]),
            "literal_250m_lon": sb_literal,
            "smoothed_guide_lon": sb,
            "guide_adjustment_km": float(rec["guide_adjustment_km"]),
            "predicted_shifted_A_coast_lon": float(sb+qcoastA/kmdeg),
            "alongshore_depth_weight": walong,
        })

    # Only a minimal one-cell wet-dead-end cleanup is needed after this surgery.
    Hminus, dead_log = remove_wet_dead_ends(Hminus, lon, lat)
    noSVB_source_mask |= (np.abs(Hminus-Hplus) > 1e-10) | ((Hminus > 0) != (Hplus > 0))
    filter_source_protect |= noSVB_source_mask

    # Invariant: continental slope is unchanged seaward of the smoothed guide.
    max_seaward_change = 0.0
    for rec in construction_rows:
        j = int(rec["j"])
        sb = float(rec["smoothed_guide_lon"])
        q = (lon-sb)*111.32*np.cos(np.deg2rad(lat[j]))
        if np.any(q < 0):
            max_seaward_change = max(
                max_seaward_change,
                float(np.max(np.abs(Hminus[j, q < 0]-Hplus[j, q < 0]))),
            )
    if max_seaward_change != 0.0:
        raise RuntimeError(f"no-SVB construction changed the slope by {max_seaward_change} m.")

    return (
        Hminus,
        noSVB_source_mask,
        filter_source_protect,
        pd.DataFrame(construction_rows),
        join_df,
        dead_log,
        max_seaward_change,
    )


# ------------------------- native-grid tangent filter ------------------------
def build_smooth_plus(Hplus: np.ndarray, lon: np.ndarray, lat: np.ndarray, filter_source_protect: np.ndarray):
    dx_km, dy_km = model_dx_dy_km(lon, lat, 30.0)
    wet = Hplus > 0

    # Orientation only: normalized 10-km Gaussian.
    sy = ORIENTATION_SIGMA_KM/dy_km
    sx = ORIENTATION_SIGMA_KM/dx_km
    num = gaussian_filter(Hplus*wet, sigma=(sy, sx), mode="nearest")
    den = gaussian_filter(wet.astype(float), sigma=(sy, sx), mode="nearest")
    Horient = np.full_like(Hplus, np.nan)
    np.divide(num, den, out=Horient, where=den > 0.25)

    gy, gx = np.gradient(Horient, dy_km, dx_km)
    gmag = np.hypot(gx, gy)
    ok = np.isfinite(gmag) & (gmag > 1e-6)
    tx = np.zeros_like(Hplus)
    ty = np.zeros_like(Hplus)
    tx[ok] = -gy[ok]/gmag[ok]
    ty[ok] =  gx[ok]/gmag[ok]

    coast_dist = distance_transform_edt(wet, sampling=(dy_km, dx_km))
    Wcoast = smoothstep((coast_dist-COAST_EXACT_KM)/(COAST_FULL_FILTER_KM-COAST_EXACT_KM))
    jA = int(np.argmin(np.abs(lat-LAT_A)))
    north_distance = np.maximum(0.0, (lat[:, None]-lat[jA])*111.32)
    Wsource = smoothstep(north_distance/SOURCE_RAMP_KM)
    Wsource[lat[:, None] < lat[jA]] = 0.0

    Wdepth = np.ones_like(Hplus)
    Wdepth[Hplus >= FILTER_ZERO_DEPTH_M] = 0.0
    tr = (Hplus > FILTER_FULL_DEPTH_M) & (Hplus < FILTER_ZERO_DEPTH_M)
    Wdepth[tr] = smoothstep(
        (FILTER_ZERO_DEPTH_M-Hplus[tr])/(FILTER_ZERO_DEPTH_M-FILTER_FULL_DEPTH_M)
    )

    W = Wcoast*Wsource*Wdepth
    W[~wet] = 0.0
    W[filter_source_protect] = 0.0
    W[~ok] = 0.0

    jj, ii = np.indices(Hplus.shape, dtype=float)
    sigma_km = LAMBDA50_POWER_KM*np.sqrt(np.log(2.0))/(2*np.pi)
    radii = np.arange(
        FILTER_SAMPLE_STEP_KM,
        float(np.ceil(3*sigma_km/FILTER_SAMPLE_STEP_KM))*FILTER_SAMPLE_STEP_KM + 0.5*FILTER_SAMPLE_STEP_KM,
        FILTER_SAMPLE_STEP_KM,
    )

    acc = Hplus.copy().astype(float)
    wsum = np.ones_like(Hplus, dtype=float)
    for sign in (-1.0, 1.0):
        alive = wet.copy()
        for r in radii:
            wt = np.exp(-0.5*(r/sigma_km)**2)
            yq = jj + sign*r*ty/dy_km
            xq = ii + sign*r*tx/dx_km
            wf = map_coordinates(wet.astype(float), [yq, xq], order=1, mode="constant", cval=0.0)
            vn = map_coordinates(Hplus, [yq, xq], order=1, mode="constant", cval=0.0)
            alive &= wf > 0.80
            q = alive
            acc[q] += wt*(vn[q]/wf[q])
            wsum[q] += wt
    Hlp = acc/wsum

    Hsmooth = Hplus + W*(Hlp-Hplus)
    Hsmooth[~wet] = 0.0
    return Hsmooth, W, tx, ty


# -------------------------- OBC + straight south -----------------------------
def apply_obc_with_straight_south(fields_pre: dict[str, np.ndarray], lon: np.ndarray, lat: np.ndarray):
    common = fields_pre["realistic_plusSVB"]
    wet = common > 0
    ny, nx = common.shape

    # Straight coast through the full 0..30 OBC zone using rows 30/31.
    wet30 = np.where(wet[N_OBC_TOTAL])[0]
    wet31 = np.where(wet[N_OBC_TOTAL+1])[0]
    if len(wet30) == 0 or len(wet31) == 0:
        raise RuntimeError("Could not identify southern mainland coast at rows 30/31.")
    i_coast = int(wet30.max())
    if int(wet31.max()) != i_coast:
        raise RuntimeError("Natural coastline differs at rows 30 and 31; cannot join straight OBC cleanly.")
    if np.any(~wet[:N_OBC_TOTAL+1, :i_coast+1]):
        raise RuntimeError("Unexpected dry hole west of proposed straight southern coastline.")

    straight_change = np.zeros_like(wet, dtype=bool)
    masked = {}
    for name, H in fields_pre.items():
        Hm = H.copy()
        q = np.zeros_like(wet, dtype=bool)
        q[:N_OBC_TOTAL+1, i_coast+1:] = True
        close = q & (Hm > 0)
        Hm[close] = 0.0
        if not np.any(straight_change):
            straight_change = close.copy()
        else:
            if not np.array_equal(close, straight_change):
                raise RuntimeError("Southern mask straightening is not identical across cases.")
        masked[name] = Hm

    common_ref = masked["realistic_plusSVB"].copy()
    common_mask = common_ref > 0
    active_w = np.all(common_mask[:, :N_OBC_TOTAL+1], axis=1)
    active_s = np.all(common_mask[:N_OBC_TOTAL+1, :], axis=0)
    active_n = np.all(common_mask[-(N_OBC_TOTAL+1):, :], axis=0)

    def one_boundary(Hcase, Hanchor_ref, boundary, active):
        H = Hcase.copy()
        if boundary == "west":
            anchor = Hanchor_ref[:, N_OBC_CORE].copy()
            rows = np.where(active)[0]
            H[np.ix_(rows, np.arange(N_OBC_CORE))] = anchor[rows, None]
            original = Hcase.copy()
            for i in range(N_OBC_CORE, N_OBC_TOTAL+1):
                w = float(smoothstep((i-N_OBC_CORE)/N_OBC_BLEND))
                H[rows, i] = (1-w)*anchor[rows] + w*original[rows, i]
        elif boundary == "south":
            anchor = Hanchor_ref[N_OBC_CORE, :].copy()
            cols = np.where(active)[0]
            H[np.ix_(np.arange(N_OBC_CORE), cols)] = anchor[None, cols]
            original = Hcase.copy()
            for j in range(N_OBC_CORE, N_OBC_TOTAL+1):
                w = float(smoothstep((j-N_OBC_CORE)/N_OBC_BLEND))
                H[j, cols] = (1-w)*anchor[cols] + w*original[j, cols]
        elif boundary == "north":
            anchor_j = ny-1-N_OBC_CORE
            anchor = Hanchor_ref[anchor_j, :].copy()
            cols = np.where(active)[0]
            for d in range(N_OBC_CORE):
                j = ny-1-d
                H[j, cols] = anchor[cols]
            original = Hcase.copy()
            for d in range(N_OBC_CORE, N_OBC_TOTAL+1):
                j = ny-1-d
                w = float(smoothstep((d-N_OBC_CORE)/N_OBC_BLEND))
                H[j, cols] = (1-w)*anchor[cols] + w*original[j, cols]
        else:
            raise ValueError(boundary)
        return H

    fields = {k: v.copy() for k, v in masked.items()}
    ref = common_ref.copy()
    for boundary, active in [("west", active_w), ("south", active_s), ("north", active_n)]:
        ref_before = ref.copy()
        fields = {k: one_boundary(H, ref_before, boundary, active) for k, H in fields.items()}
        ref = one_boundary(ref_before, ref_before, boundary, active)

    for name in fields:
        if not np.array_equal(fields[name] > 0, masked[name] > 0):
            raise RuntimeError(f"OBC depth treatment changed the mask in {name}.")

    return fields, straight_change, i_coast, active_w, active_s, active_n


# -------------------------- exact partial-cell audit -------------------------
def hfac_exact_2d(H: np.ndarray):
    rF = np.r_[0.0, np.cumsum(DEL_R)]
    wet_input = H > 0
    Heff = np.zeros_like(H)
    wetcol = np.zeros_like(wet_input)
    bottom_hfac = np.zeros_like(H)
    bottom_thickness = np.zeros_like(H)
    bottom_k = np.full(H.shape, -1, dtype=np.int16)
    min_open_W = np.inf
    min_open_S = np.inf

    for k, dz in enumerate(DEL_R):
        hmin = max(HFAC_MIN, min(HFAC_MIN_DR_M/dz, 1.0))
        raw = np.clip((H-rF[k])/dz, 0.0, 1.0)
        hC = raw.copy()
        hC[hC < hmin/2.0] = np.nan
        hC[~wet_input] = np.nan
        hC = np.maximum(hC, hmin)
        hC[np.isnan(hC)] = 0.0

        Heff += hC*dz
        q = hC > 0
        wetcol |= q
        bottom_hfac[q] = hC[q]
        bottom_thickness[q] = hC[q]*dz
        bottom_k[q] = k

        # useMin4hFacEdges=.TRUE.
        hW = np.minimum(hC[:, 1:], hC[:, :-1])
        hS = np.minimum(hC[1:, :], hC[:-1, :])
        if np.any(hW > 0):
            min_open_W = min(min_open_W, float(np.min(hW[hW > 0]*dz)))
        if np.any(hS > 0):
            min_open_S = min(min_open_S, float(np.min(hS[hS > 0]*dz)))

    return {
        "effective_depth": Heff,
        "wet_column": wetcol,
        "bottom_hfac": bottom_hfac,
        "bottom_thickness": bottom_thickness,
        "bottom_k": bottom_k,
        "min_open_W_face_thickness_m": min_open_W,
        "min_open_S_face_thickness_m": min_open_S,
    }


# --------------------------------- output -----------------------------------
def write_outputs(
    outdir: Path,
    smooth_label: str,
    lon: np.ndarray,
    lat: np.ndarray,
    H_gebco_direct: np.ndarray,
    H_artifact_clean: np.ndarray,
    H_common: np.ndarray,
    preobc: dict[str, np.ndarray],
    final: dict[str, np.ndarray],
    noSVB_source_mask: np.ndarray,
    filter_source_protect: np.ndarray,
    filter_weight: np.ndarray,
    tx: np.ndarray,
    ty: np.ndarray,
    straight_change: np.ndarray,
    construction: pd.DataFrame,
    join_df: pd.DataFrame,
    exact_diag: dict[str, dict],
    audits: pd.DataFrame,
    partial_summary: pd.DataFrame,
):
    outdir.mkdir(parents=True, exist_ok=True)
    rF = np.r_[0.0, np.cumsum(DEL_R)]
    rC = rF[:-1] + 0.5*DEL_R

    ds = xr.Dataset(
        data_vars={
            "Depth_GEBCO_direct": (("YC", "XC"), H_gebco_direct.astype(np.float32)),
            "Depth_common_artifact_clean": (("YC", "XC"), H_artifact_clean.astype(np.float32)),
            "Depth_common_topology_clean": (("YC", "XC"), H_common.astype(np.float32)),
            "Depth_realistic_plusSVB_preOBC": (("YC", "XC"), preobc["realistic_plusSVB"].astype(np.float32)),
            f"Depth_{smooth_label}_plusSVB_preOBC": (("YC", "XC"), preobc[f"{smooth_label}_plusSVB"].astype(np.float32)),
            "Depth_realistic_minusSVB_preOBC": (("YC", "XC"), preobc["realistic_minusSVB"].astype(np.float32)),
            f"Depth_{smooth_label}_minusSVB_preOBC": (("YC", "XC"), preobc[f"{smooth_label}_minusSVB"].astype(np.float32)),
            "Depth_realistic_plusSVB": (("YC", "XC"), final["realistic_plusSVB"].astype(np.float32)),
            f"Depth_{smooth_label}_plusSVB": (("YC", "XC"), final[f"{smooth_label}_plusSVB"].astype(np.float32)),
            "Depth_realistic_minusSVB": (("YC", "XC"), final["realistic_minusSVB"].astype(np.float32)),
            f"Depth_{smooth_label}_minusSVB": (("YC", "XC"), final[f"{smooth_label}_minusSVB"].astype(np.float32)),
            "noSVB_source_mask": (("YC", "XC"), noSVB_source_mask.astype(np.int8)),
            "filter_source_protected_mask": (("YC", "XC"), filter_source_protect.astype(np.int8)),
            "filter_weight": (("YC", "XC"), filter_weight.astype(np.float32)),
            "orientation_tx": (("YC", "XC"), tx.astype(np.float32)),
            "orientation_ty": (("YC", "XC"), ty.astype(np.float32)),
            "southern_straight_coast_mask_change": (("YC", "XC"), straight_change.astype(np.int8)),
            "delR": (("Z",), DEL_R.astype(np.float64)),
            "rF_depth": (("Zp1",), rF.astype(np.float64)),
            "rC_depth": (("Z",), rC.astype(np.float64)),
        },
        coords={"XC": lon, "YC": lat, "Z": np.arange(1, 101), "Zp1": np.arange(0, 101)},
        attrs={
            "description": "Final CTW bathymetry suite reproduced from GEBCO.",
            "smooth_filter_lambda50_power_km": LAMBDA50_POWER_KM,
            "noSVB_shelfbreak_reference_m": SHELF_BREAK_M,
            "noSVB_spike_bridge_lat_south": SPIKE_BRIDGE_SOUTH,
            "noSVB_spike_bridge_lat_north": SPIKE_BRIDGE_NORTH,
            "minimum_positive_depth_m": HMIN_FINAL_M,
            "hFacMin": HFAC_MIN,
            "hFacMinDr_m": HFAC_MIN_DR_M,
            "useMin4hFacEdges": "TRUE",
            "pCellMix_select": PCELL_MIX_SELECT,
            "MITgcm_binary_convention": "big-endian float64; ocean depths negative; land zero",
        },
    )
    final_nc = outdir/"CTW_final_bathymetries.nc"
    ds.to_netcdf(final_nc)

    # Final MITgcm bathymetry binaries.
    binary_paths = []
    for name, H in final.items():
        p = outdir/f"bathy_{name}_FINAL.bin"
        (-H).astype(">f8").tofile(p)
        binary_paths.append(p)

    # Exact 2-D partial-cell diagnostics.
    diag_vars = {
        "delR": (("Z",), DEL_R.astype(np.float64)),
        "rF_depth": (("Zp1",), rF.astype(np.float64)),
        "rC_depth": (("Z",), rC.astype(np.float64)),
    }
    for name, H in final.items():
        d = exact_diag[name]
        diag_vars[f"Depth_input_{name}"] = (("YC", "XC"), H.astype(np.float32))
        diag_vars[f"Depth_effective_{name}"] = (("YC", "XC"), d["effective_depth"].astype(np.float32))
        diag_vars[f"wetColumn_{name}"] = (("YC", "XC"), d["wet_column"].astype(np.int8))
        diag_vars[f"bottom_hFacC_{name}"] = (("YC", "XC"), d["bottom_hfac"].astype(np.float32))
        diag_vars[f"bottom_thickness_m_{name}"] = (("YC", "XC"), d["bottom_thickness"].astype(np.float32))
        diag_vars[f"bottom_k_0based_{name}"] = (("YC", "XC"), d["bottom_k"].astype(np.int16))
    dsd = xr.Dataset(
        data_vars=diag_vars,
        coords={"XC": lon, "YC": lat, "Z": np.arange(1, 101), "Zp1": np.arange(0, 101)},
        attrs={
            "hFacMin": HFAC_MIN,
            "hFacMinDr_m": HFAC_MIN_DR_M,
            "minimum_positive_depth_m": HMIN_FINAL_M,
            "useMin4hFacEdges": "TRUE",
            "status": "Exact 2-D partial-cell audit using production delR",
        },
    )
    diag_nc = outdir/"CTW_final_partial_cell_diagnostics.nc"
    dsd.to_netcdf(diag_nc)

    construction.to_csv(outdir/"noSVB_final_smoothed_guide.csv", index=False)
    join_df.to_csv(outdir/"noSVB_southern_join_candidates.csv", index=False)
    audits.to_csv(outdir/"final_topology_and_consistency_audit.csv", index=False)
    partial_summary.to_csv(outdir/"final_partial_cell_summary.csv", index=False)

    # Final namelist fragment.
    vals = [f"{x:g}" for x in DEL_R]
    lines = [" " + ", ".join(vals[i:i+10]) + "," for i in range(0, 100, 10)]
    snippet = (
        "# FINAL CTW vertical-grid / partial-cell settings\n\n"
        "&PARM01\n"
        f" hFacMin   = {HFAC_MIN:.2f},\n"
        f" hFacMinDr = {HFAC_MIN_DR_M:.1f},\n"
        "/\n\n"
        "&PARM04\n"
        " delR =\n" + "\n".join(lines) + "\n"
        " useMin4hFacEdges = .TRUE.,\n"
        f" pCellMix_select   = {PCELL_MIX_SELECT},\n"
        "/\n"
    )
    (outdir/"FINAL_data_partial_cells_and_delR.txt").write_text(snippet)
    return final_nc, diag_nc, binary_paths


# ---------------------------------- main ------------------------------------
def main() -> None:
    global LAMBDA50_POWER_KM
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gebco", type=Path, required=True, help="GEBCO_2020 subset NetCDF")
    parser.add_argument("--grid", type=Path, required=True, help="MITgcm grid NetCDF containing 1-D XC, YC")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory")
    parser.add_argument(
        "--smooth-length-km",
        type=float,
        default=LAMBDA50_POWER_KM,
        help="Along-isobath smoothing length in km for the smooth bathymetry.",
    )
    parser.add_argument(
        "--smooth-label",
        default=None,
        help="Output label for the smooth bathymetry, e.g. smooth30 or smooth70.",
    )
    args = parser.parse_args()
    LAMBDA50_POWER_KM = args.smooth_length_km
    smooth_label = args.smooth_label or f"smooth{args.smooth_length_km:g}".replace(".", "p")

    if len(DEL_R) != 100 or abs(DEL_R.sum()-4700.0) > 1e-8:
        raise RuntimeError("Production delR must contain 100 levels summing to 4700 m.")

    print("[1/7] GEBCO -> common cleaned H0")
    lon, lat, Hdirect, Hart, H0, puddle_log, topo_logs = build_common_H0(args.gebco, args.grid)

    print("[2/7] Revised shelf-only -SVB geometry with 28.6 N spike removed")
    Hminus, noSVB_source_mask, filter_source_protect, construction, join_df, dead_log, max_seaward = build_noSVB_shelf_only(H0, lon, lat)

    print("[3/7] Seam-free native-grid tangent filter")
    Hsmooth_plus, W, tx, ty = build_smooth_plus(H0, lon, lat, filter_source_protect)
    # Exact revised -SVB source inside the source patch; common smooth propagation outside.
    Hsmooth_minus = Hsmooth_plus.copy()
    Hsmooth_minus[noSVB_source_mask] = Hminus[noSVB_source_mask]
    Hsmooth_minus[~(Hminus > 0)] = 0.0
    if not np.array_equal(Hsmooth_minus > 0, Hminus > 0):
        raise RuntimeError("Smooth -SVB mask differs from realistic -SVB mask before OBC.")
    if np.max(np.abs(Hsmooth_plus-H0)[noSVB_source_mask & (H0 > 0)]) != 0.0:
        raise RuntimeError("Smooth +SVB is not exactly protected over the revised no-SVB source patch.")

    preobc = {
        "realistic_plusSVB": H0,
        f"{smooth_label}_plusSVB": Hsmooth_plus,
        "realistic_minusSVB": Hminus,
        f"{smooth_label}_minusSVB": Hsmooth_minus,
    }

    print("[4/7] Straight southern mask + common OBC depth treatment")
    postobc, straight_change, i_coast, active_w, active_s, active_n = apply_obc_with_straight_south(preobc, lon, lat)

    print("[5/7] Common minimum positive depth Hmin=2 m")
    final = {}
    for name, H0case in postobc.items():
        H = H0case.copy()
        q = (H > 0) & (H < HMIN_FINAL_M)
        H[q] = HMIN_FINAL_M
        final[name] = H

    print("[6/7] Exact hFac / topology / consistency audit")
    exact_diag = {}
    audit_rows = []
    partial_rows = []
    for name, H in final.items():
        d = hfac_exact_2d(H)
        exact_diag[name] = d
        wet = H > 0
        a = topology_audit(d["wet_column"])
        a.update({
            "case": name,
            "mask_differences_input_vs_hFac": int(np.sum(wet != d["wet_column"])),
        })
        audit_rows.append(a)
        err = d["effective_depth"][wet]-H[wet]
        partial_rows.append({
            "case": name,
            "input_wet_columns": int(wet.sum()),
            "newly_dry_columns": int(np.sum(wet & ~d["wet_column"])),
            "min_bottom_thickness_m": float(np.min(d["bottom_thickness"][wet])),
            "p01_bottom_thickness_m": float(np.percentile(d["bottom_thickness"][wet], 1)),
            "median_bottom_thickness_m": float(np.median(d["bottom_thickness"][wet])),
            "min_bottom_hFacC": float(np.min(d["bottom_hfac"][wet])),
            "max_abs_effective_depth_error_m": float(np.max(np.abs(err))),
            "p95_abs_effective_depth_error_m": float(np.percentile(np.abs(err), 95)),
            "median_abs_effective_depth_error_m": float(np.median(np.abs(err))),
            "min_open_W_face_thickness_m": d["min_open_W_face_thickness_m"],
            "min_open_S_face_thickness_m": d["min_open_S_face_thickness_m"],
        })
    audits = pd.DataFrame(audit_rows)
    partial_summary = pd.DataFrame(partial_rows)

    # Critical pairwise consistency checks.
    for pair, rname, sname in [
        ("plusSVB", "realistic_plusSVB", f"{smooth_label}_plusSVB"),
        ("minusSVB", "realistic_minusSVB", f"{smooth_label}_minusSVB"),
    ]:
        if not np.array_equal(final[rname] > 0, final[sname] > 0):
            raise RuntimeError(f"Final realistic/smooth masks differ for {pair}.")
        if np.max(np.abs(exact_diag[rname]["effective_depth"][:15]-exact_diag[sname]["effective_depth"][:15])) != 0.0:
            raise RuntimeError(f"Final effective southern OBC depths differ for {pair}.")

    if (audits[["mask_differences_input_vs_hFac", "dangling_wet_columns", "enclosed_dry_islands_size_1_2", "near_land_dry_islands_size_3_5"]].to_numpy() != 0).any():
        raise RuntimeError("Final partial-cell topology audit failed.")

    print("[7/7] Writing NetCDFs, binaries, and audit tables")
    final_nc, diag_nc, binaries = write_outputs(
        args.outdir, smooth_label, lon, lat, Hdirect, Hart, H0, preobc, final,
        noSVB_source_mask, filter_source_protect, W, tx, ty,
        straight_change, construction, join_df, exact_diag, audits, partial_summary,
    )

    # Additional compact build summary.
    summary = pd.DataFrame([{
        "grid_nx": len(lon),
        "grid_ny": len(lat),
        "straight_south_last_wet_i": i_coast,
        "straight_south_coast_lon": float(lon[i_coast]),
        "west_active_rows": int(active_w.sum()),
        "south_active_cols": int(active_s.sum()),
        "north_active_cols": int(active_n.sum()),
        "max_noSVB_change_seaward_smoothed_guide_m": max_seaward,
        "max_noSVB_guide_adjustment_km": float(np.max(np.abs(construction["guide_adjustment_km"]))),
        "hFacMin": HFAC_MIN,
        "hFacMinDr_m": HFAC_MIN_DR_M,
        "minimum_positive_depth_m": HMIN_FINAL_M,
    }])
    summary.to_csv(args.outdir/"final_build_summary.csv", index=False)

    # Save cleanup logs.
    logframes = []
    for key, d in topo_logs.items():
        if len(d):
            q = d.copy(); q["stage"] = "common_H0"; q["pass"] = key; logframes.append(q)
    if len(dead_log):
        q = dead_log.copy(); q["stage"] = "revised_noSVB"; q["pass"] = "wet_dead_ends"; logframes.append(q)
    if logframes:
        pd.concat(logframes, ignore_index=True, sort=False).to_csv(args.outdir/"topology_cleanup_log.csv", index=False)
    else:
        pd.DataFrame().to_csv(args.outdir/"topology_cleanup_log.csv", index=False)
    puddle_log.to_csv(args.outdir/"initial_disconnected_puddles.csv", index=False)

    print("\nFinal outputs")
    print(final_nc)
    print(diag_nc)
    for p in binaries:
        print(p)
    print(args.outdir/"final_build_summary.csv")
    print(args.outdir/"final_topology_and_consistency_audit.csv")
    print(args.outdir/"final_partial_cell_summary.csv")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        main()
