# MITgcm Configuration Summary for the CTW-SVB Reruns

## Model Domain and Grid

The final rerun suite uses MITgcm on a spherical-polar C-grid with
`Nx x Ny x Nr = 576 x 684 x 100` grid cells. The grid origin supplied to
MITgcm is approximately:

```text
(xgOrigin, ygOrigin) = (236.36437 deg E, 25.39553 deg N)
```

This is the original `512 x 612` active grid augmented by 64 telescoped columns
to the west and 72 telescoped rows to the south. The added wet buffer is a flat
4400 m abyssal plane with straight artificial coastlines.

Horizontal spacing ranges from `0.0155945 deg` to `0.0389862 deg` in longitude
and from `0.0135994 deg` to `0.0339985 deg` in latitude. The vertical grid has
100 levels and reaches 4700 m depth.

## Bathymetry

The experiment suite uses four bathymetry files that differ only in the resolved
topography and presence or absence of Sebastian Vizcaino Bay (SVB):

| Run | Bathymetry | SVB |
| --- | --- | --- |
| `Rplus` | realistic | present |
| `Rminus` | realistic | removed |
| `Splus` | smoothed | present |
| `Sminus` | smoothed | removed |

The MITgcm bathymetry files are big-endian binary arrays with negative ocean
depth and zero over land.

The realistic +SVB bathymetry is the cleaned baseline bathymetry. The -SVB
members modify only the shelf/bay source geometry: the continental slope is
preserved seaward of the shelf-break guide, and the 250 m shelf-break spike near
28.6 deg N is regularized to avoid an artificial land spike.

The smoothed members use the approved 30 km half-power native-grid local-tangent
smoothing, not the historical repeated smoother. Realistic and smoothed members
retain the same wet/dry mask within each +SVB or -SVB pair, and the immediate
coastal depths are held fixed so that the comparison isolates short-scale
bathymetric roughness rather than coastline or partial-cell changes.

For the west64/south72 reruns, these four finalized bathymetries were inserted
into the expanded domain. The added western and southern buffer regions are flat
4400 m abyssal ocean at the open boundaries, with straight artificial
coastlines. The western added buffer is smoothly blended toward the local
western-edge depths, and the southern added buffer is smoothly blended toward
the first active-domain row, so the southwest, southeast, and northwest sponge
corners are continuous. In the northern buffer/taper band, existing wet cells
continue the local depth from the row immediately south of the band while
preserving the land mask and the configured north-boundary wet span. The
southwest added-buffer corner is additionally smoothed with a local 2-D wet-cell
relaxation, and the northwest added-buffer corner uses a smoothed
along-boundary profile while retaining a gradual northward transition. The
open-boundary bathymetry treatment in the finalized source bathymetries uses a
15-cell flat normal-depth core with a 15-cell smooth transition on the west,
south, and north boundaries.

![Production bathymetries and diagnostic differences](figures/bathymetry_panels.png)

## Time Stepping and Output

The simulations start from rest at `t = 0` and run for 1,296,000 s, or 15 days,
with `deltaT = 60 s`. Checkpoints are written every 43,200 s and model dumps
every 172,800 s. Diagnostics are written in MITgcm MDS binary format rather than
MNC NetCDF.

| File prefix | Frequency | Variables |
| --- | ---: | --- |
| `dyn3d_avg` | 1200 s averages | `UVEL`, `VVEL`, `WVEL`, `RHOAnoma`, `PHIHYD` |
| `surf_avg` | 1200 s averages | `ETAN`, `PHI_SURF` |
| `surf_snap` | 1200 s snapshots | `ETAN`, `PHI_SURF` |

The 3-D diagnostics are saved on 75 vertical levels:

```text
1-65, 67, 69, 71, 73, 75, 77, 79, 84, 90, 100
```

## Equation of State, Initial Conditions, and Open Boundaries

The runs use a linear equation of state with:

```text
rhoNil = 999.8 kg m^-3
tAlpha = 2.0e-4 degC^-1
sBeta  = 7.4e-4
```

MITgcm reference values are set to `tRef = 100*0.` and `sRef = 100*35.`, so
that the linear equation of state is satisfied by constant reference values. The
initial temperature and salinity fields are vertically varying and horizontally
uniform. The same profiles are prescribed at the open boundaries using OBCS
temperature and salinity files. No open-boundary sponge relaxation is applied to
temperature or salinity.

![Initial temperature, salinity, and density profiles](figures/initial_profiles.png)

## Open Boundaries and Sponge

Open boundaries are active on the west, south, and north boundaries. The OBCS
masks are:

```text
OB_Iwest  = 684*1
OB_Jsouth = 543*1,33*0
OB_Jnorth = 137*-1,439*0
```

Orlanski radiation is enabled on all three open boundaries, with
`cvelTimeScale = 2000.0` and `CMAX = 1.0`. OBCS volume balance is enabled with
unit balance factors on west, south, and north boundaries.

A linear OBCS sponge is active for momentum only. The sponge thickness is 48
grid cells, with inner relaxation time scale 43,200 s and boundary relaxation
time scale 1,800 s for both `u` and `v`. The sponge acts on the west, south, and
north boundaries and relaxes both east-west and north-south velocity components.
Temperature and salinity sponge switches are disabled.

## Forcing and Experiment Set

Wind stress is periodic with a 1,800 s forcing period and a 1,382,400 s forcing
cycle. The forcing files prescribe zonal and meridional wind stress over the
`576 x 684` domain. The maximum wind-stress magnitude over the full forcing
cycle is `0.08934 N m^-2`, occurring at hour 3.0 of the forcing record.

The four experiments are `Rplus`, `Rminus`, `Splus`, and `Sminus`, as described
in the bathymetry section above.

![Maximum wind-stress magnitude and vectors](figures/windstress_max_map.png)
