#!/usr/bin/env python3
"""Create the four CTW/SVB rerun input directories from one template."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "real_SVB"
RERUNS = ROOT / "ctw_reruns"

RUNS = {
    "Rplus": "bathy_realistic_plusSVB_FINAL.bin",
    "Rminus": "bathy_realistic_minusSVB_FINAL.bin",
    "Splus": "bathy_smooth30_plusSVB_FINAL.bin",
    "Sminus": "bathy_smooth30_minusSVB_FINAL.bin",
}

LEVELS = list(range(1, 66)) + [67, 69, 71, 73, 75, 77, 79, 84, 90, 100]
LEVEL_LINES = """  levels(1:20,1) = 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.,
                   11., 12., 13., 14., 15., 16., 17., 18., 19., 20.,
  levels(21:40,1) = 21., 22., 23., 24., 25., 26., 27., 28., 29., 30.,
                    31., 32., 33., 34., 35., 36., 37., 38., 39., 40.,
  levels(41:60,1) = 41., 42., 43., 44., 45., 46., 47., 48., 49., 50.,
                    51., 52., 53., 54., 55., 56., 57., 58., 59., 60.,
  levels(61:75,1) = 61., 62., 63., 64., 65., 67., 69., 71., 73., 75.,
                    77., 79., 84., 90., 100.,"""
NX, NY, NR = 512, 612, 100

DATA_TEMPLATE = """# Model parameters
# Continuous equation parameters
 &PARM01
 viscAh=100.,
 viscAr=1.E-5,
 no_slip_sides=.FALSE.,
 no_slip_bottom=.FALSE.,
 diffKhT=100.,
 diffKrT=1.E-5,
 diffKhS=100.,
 diffKrS=1.E-5,
 ivdc_kappa=1.,
 implicitDiffusion=.TRUE.,
 eosType='LINEAR',
 tAlpha=2.E-4,
 sBeta=7.4E-4,
 rhoNil=999.8,
 gravity=9.81,
 rigidLid=.FALSE.,
 implicitFreeSurface=.TRUE.,
 exactConserv=.TRUE.,
 saltStepping=.TRUE.,
 readBinaryPrec = 64,
 writeBinaryPrec = 64,
 hFacMin   = 0.10,
 hFacMinDr = 2.0,
 tRef=100*0.,
 sRef=100*35.,
 useSingleCpuIO=.TRUE.,
 &

# Elliptic solver parameters
 &PARM02
 cg2dTargetResidual=1.E-7,
 cg2dMaxIters=1000,
 &

# Time stepping parameters
 &PARM03
 startTime=0.,
 endTime=1296000.,
# pickupSuff='ckptB',
 nIter0=0,
 deltaT=60.,
 chkptFreq=43200.,
 dumpFreq=172800.,
 monitorFreq=1200.,
 monitorSelect=2,
 periodicExternalForcing=.TRUE.,
 externForcingPeriod=1800.,
 externForcingCycle=1382400.,
 &

# Gridding parameters
 &PARM04
 usingSphericalPolarGrid=.TRUE.,
 xgOrigin=238.0,
 ygOrigin=27.0,
 delR= 2, 2, 2, 2, 2, 2.2, 2.4, 2.6, 2.8,
       3.0, 3.2, 3.5, 3.8, 4.2, 4.5, 4.8, 5, 5, 5,
       5, 5, 5, 5.5, 5.5, 5.5, 5.5, 5.5, 6, 6.5,
       7, 7.5, 8, 8.5, 9.2, 10, 10, 10, 10, 10, 10,
       10, 11, 12, 13, 14, 15, 16, 17, 18, 20,
       22, 23, 25, 25, 25, 25, 25, 25, 25, 25,
       27, 28, 30, 30, 32, 34, 36, 38, 40, 43,
       46, 50, 55, 59, 65, 70, 75, 80, 85, 92,
       100, 100, 100, 100, 100, 100, 105, 115, 125, 132,
       145, 160, 170, 180, 190, 200, 210, 230, 245, 269.8,
 delY = 612*0.0136,
 delX = 512*0.0156,
 useMin4hFacEdges = .TRUE.,
 pCellMix_select  = 0,
 &

 &PARM05
 bathyFile='../bathymetry/{bathy}',
 meridWindFile='../common_input/upwelling_merid_halfgauss2hr_16day_512x612.bin',
 zonalWindFile='../common_input/upwelling_zonal_halfgauss2hr_16day_512x612.bin',
 hydrogThetaFile='../common_input/feb_Tref_100zlev_512x612.bin',
 hydrogSaltFile='../common_input/feb_Sref_100zlev_512x612.bin',
 &
"""

DATA_OBCS = """# Open-boundaries
 &OBCS_PARM01
 OB_Iwest = 612*1,
 OB_Jsouth = 479*1,33*0,
 OB_Jnorth = 73*-1,439*0,
 useOrlanskiWest=.TRUE.,
 useOrlanskiSouth=.TRUE.,
 useOrlanskiNorth=.TRUE.,
 useOBCSsponge=.FALSE.,
 useLinearSponge=.FALSE.,
 OBCSsponge_Salt=.FALSE.,
 OBCSsponge_Theta=.FALSE.,
 useOBCSbalance=.TRUE.,
 useOBCSprescribe=.TRUE.,
 OBWtFile='../common_input/obcs_ts_from_initial_hydrography_16day/OBWt.bin',
 OBStFile='../common_input/obcs_ts_from_initial_hydrography_16day/OBSt.bin',
 OBNtFile='../common_input/obcs_ts_from_initial_hydrography_16day/OBNt.bin',
 OBWsFile='../common_input/obcs_ts_from_initial_hydrography_16day/OBWs.bin',
 OBSsFile='../common_input/obcs_ts_from_initial_hydrography_16day/OBSs.bin',
 OBNsFile='../common_input/obcs_ts_from_initial_hydrography_16day/OBNs.bin',
 OBCS_balanceFacW=1.,
 OBCS_balanceFacS=1.,
 OBCS_balanceFacN=1.,
 OBCS_monitorFreq=1200.,
 &

 &OBCS_PARM02
 cvelTimeScale=2000.0,
 CMAX=1.0,
 &

 &OBCS_PARM03
 spongeThickness=12,
 Urelaxobcsinner=600.E0,
 Vrelaxobcsinner=600.E0,
 Urelaxobcsbound=90.E0,
 Vrelaxobcsbound=90.E0,
 &
"""

DATA_PKG = """# Packages (lines beginning "#" are comments)
 &PACKAGES
 useMNC=.TRUE.,
 useDiagnostics=.TRUE.,
 useOBCS=.TRUE.,
 &
"""

DATA_MNC = """# NetCDF/MNC output
 &MNC_01
 mnc_use_outdir=.TRUE.,
 mnc_outdir_str='mnc_{run}_',
 monitor_mnc=.FALSE.,
 pickup_write_mnc=.FALSE.,
 pickup_read_mnc=.FALSE.,
 &
"""

DATA_EXCH2 = """# EXCH2 Package: Wrapper-2 User Choice
#--------------------
#  preDefTopol   :: pre-defined Topology selector:
#                :: = 0 : topology defined from processing "data.exch2";
#                :: = 1 : simple, single facet topology;
#                :: = 2 : customized topology (w2_set_myown_facets)
#                :: = 3 : 6-facet Cube (3 face-dims: nRed, nGreen, nBlue).
#  dimsFacets    :: facet pair of dimensions (n1x,n1y, n2x,n2y ...)
#  facetEdgeLink :: Face-Edge connectivity map:
#    facetEdgeLink(i,j)=XX.1 : face(j)-edge(i) (i=1,2,3,4 <==> N,S,E,W)
#    is connected to Northern edge of face "XX" ; similarly,
#    = XX.2 : to Southern.E, XX.3 = Eastern.E, XX.4 = Western.E of face "XX"
#  blankList     :: List of "blank" tiles
#  W2_mapIO      :: global map IO selector (-1 = old type ; 0 = 1 long line in X
#                :: 1 = compact, mostly in Y dir)
#  W2_printMsg   :: option for information messages printing
#                :: <0 : write to log file ; =0 : minimum print ;
#                :: =1 : no duplicated print ; =2 : all processes do print
#--------------------
 &W2_EXCH2_PARM01
# W2_printMsg= 1,
# W2_mapIO   = 1,
  dimsFacets = 512, 612,
 &
"""

DATA_DIAGNOSTICS = f"""# Diagnostic Package Choices
 &DIAGNOSTICS_LIST
  diag_mnc = .TRUE.,
  fields(1:5,1)  = 'UVEL    ','VVEL    ','WVEL    ','RHOAnoma',
                   'PHIHYD  ',
{LEVEL_LINES}
  fileName(1) = 'dyn3d_avg',
  frequency(1) = 1200.,

  fields(1:2,2) = 'ETAN    ','PHI_SURF',
  fileName(2) = 'surf_snap',
  frequency(2) = -1200.,

  fields(1:2,3) = 'ETAN    ','PHI_SURF',
  fileName(3) = 'surf_avg',
  frequency(3) = 1200.,
 &

 &DIAG_STATIS_PARMS
 &
"""


README = """# CTW/SVB production reruns

This directory contains four MITgcm input directories generated from
`../real_SVB` using `../scripts/setup_ctw_reruns.py`.

Runs:

- `Rplus`: realistic bathymetry with SVB
- `Rminus`: realistic bathymetry without SVB
- `Splus`: smooth30 bathymetry with SVB
- `Sminus`: smooth30 bathymetry without SVB

All four runs use the same code and runtime settings. The intentional
scientific difference is only `bathyFile`.

Runtime/forcing:

- `endTime=1296000.` for a 15-day run.
- `periodicExternalForcing=.TRUE.`, `externForcingPeriod=1800.`, and
  `externForcingCycle=1382400.` use a 16-day forcing cycle. The extra day is a
  zero-wind guard interval so interpolation at day 15 does not wrap back to the
  first nonzero wind record.
- Wind-stress files `upwelling_*_halfgauss2hr_16day_512x612.bin`
  preserve the original 120 stress records in N/m^2, then
  pad zeros through the 16-day cycle.
- Generate/rebuild these files with
  `python scripts/extend_forcing_for_15day_run.py`.

MPI/build settings:

- shared `../code/SIZE.h` is configured for 16 MPI ranks:
  `sNx=64, sNy=51, nSx=2, nSy=3, nPx=4, nPy=4`
- `pkg/mnc` is compiled/enabled so diagnostics and model output can be written
  as NetCDF.

Diagnostics:

- 3-D 20-minute averages: `UVEL`, `VVEL`, `WVEL`, `RHOAnoma`, `PHIHYD`
- 2-D 20-minute snapshots and averages: `ETAN`, `PHI_SURF`
- 75 saved vertical levels: all levels 1-65, every other level from 67 to 79,
  then 84, 90, and 100. The exact grid has 79 levels with centers shallower
  than 1500 m, so a strict 75-level cap cannot include every level above 1500 m.

Temperature/salinity handling:

- `hydrogThetaFile` and `hydrogSaltFile` point to the original 3-D initial
  hydrography files in `common_input`.
- `tRef=100*0.` and `sRef=100*35.` define a constant linear-EOS reference
  state; with `rhoNil=999.8`, this means `rhoNil` is the reference density at
  T=0 degC and S=35 in the model's linear salinity convention.
- Open-boundary T/S are prescribed from the same vertically varying,
  horizontally homogeneous initial profiles using:
  `OBWtFile`, `OBStFile`, `OBNtFile`, `OBWsFile`, `OBSsFile`, and `OBNsFile`.
- The prescribed T/S OBC files are repeated for the same 16-day forcing cycle
  because OBCS shares MITgcm's periodic external-forcing loader in this setup.
- No velocity or sea-level OBC files are prescribed. Orlanski settings are
  retained for the dynamical boundary variables.

Important audit note: the present bathymetry files have the expected sizes,
wet masks, and big-endian float64 interpretation, but their SHA256 checksums
do not match the checksums listed in `../README_CODEX_CTW_RERUNS.md`.

OBC mask audit:

- west boundary: 612 wet rows
- south boundary: 479 wet columns, configured as `479*1,33*0`
- north boundary: the README/build CSV says 71 active columns, but the actual
  final bathymetry north row has 73 wet columns in all four files. The rerun
  configs use the actual wet boundary mask: `73*-1,439*0`.

External-drive output:

- Prepare external execution directories with
  `python scripts/prepare_external_output_dirs.py`
- By default this creates `/Volumes/LaCie/CTW_SVB_reruns/<case>`.
- Run MITgcm from those per-case directories so diagnostics, NetCDF/MNC output,
  checkpoints, monitor files, and any binary dumps land on the external drive.
"""


def link_force(dst: Path, target: Path) -> None:
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    rel = os.path.relpath(target, dst.parent)
    dst.symlink_to(rel)


def write_obcs_ts() -> None:
    outdir = RERUNS / "common_input" / "obcs_ts_from_initial_hydrography"
    outdir.mkdir(exist_ok=True)
    theta = np.fromfile(ROOT / "input" / "feb_Tref_100zlev_512x612.bin", dtype=">f8").reshape(NR, NY, NX)
    salt = np.fromfile(ROOT / "input" / "feb_Sref_100zlev_512x612.bin", dtype=">f8").reshape(NR, NY, NX)
    tprof = theta[:, 0, 0]
    sprof = salt[:, 0, 0]
    if np.max(np.abs(theta - tprof[:, None, None])) > 0:
        raise ValueError("Theta initial condition is not horizontally homogeneous")
    if np.max(np.abs(salt - sprof[:, None, None])) > 0:
        raise ValueError("Salt initial condition is not horizontally homogeneous")
    arrays = {
        "OBWt.bin": np.repeat(tprof[:, None], NY, axis=1),
        "OBWs.bin": np.repeat(sprof[:, None], NY, axis=1),
        "OBSt.bin": np.repeat(tprof[:, None], NX, axis=1),
        "OBSs.bin": np.repeat(sprof[:, None], NX, axis=1),
        "OBNt.bin": np.repeat(tprof[:, None], NX, axis=1),
        "OBNs.bin": np.repeat(sprof[:, None], NX, axis=1),
    }
    for name, arr in arrays.items():
        arr.astype(">f8").tofile(outdir / name)


def main() -> None:
    RERUNS.mkdir(exist_ok=True)
    (RERUNS / "bathymetry").mkdir(exist_ok=True)
    (RERUNS / "common_input").mkdir(exist_ok=True)
    (RERUNS / "scripts").mkdir(exist_ok=True)
    (RERUNS / "build").mkdir(exist_ok=True)

    for name in [
        "upwelling_merid_halfgauss2hr_512x612.bin",
        "upwelling_zonal_halfgauss2hr_512x612.bin",
        "feb_Tref_100zlev_512x612.bin",
        "feb_Sref_100zlev_512x612.bin",
    ]:
        link_force(RERUNS / "common_input" / name, ROOT / "input" / name)

    for bathy in RUNS.values():
        link_force(RERUNS / "bathymetry" / bathy, ROOT / "input" / bathy)

    write_obcs_ts()

    for name in [
        "estimate_diagnostics_storage.py",
        "extend_forcing_for_15day_run.py",
        "make_obcs_ts_from_hydrography.py",
        "prepare_external_output_dirs.py",
        "validate_bathy.py",
        "validate_configs.py",
    ]:
        link_force(RERUNS / "scripts" / name, ROOT / "scripts" / name)

    for run, bathy in RUNS.items():
        run_dir = RERUNS / run
        run_dir.mkdir(exist_ok=True)
        shutil.copy2(SOURCE / "eedata", run_dir / "eedata")
        (run_dir / "data.exch2").write_text(DATA_EXCH2)
        (run_dir / "data.pkg").write_text(DATA_PKG)
        (run_dir / "data.mnc").write_text(DATA_MNC.format(run=run))
        (run_dir / "data").write_text(DATA_TEMPLATE.format(bathy=bathy))
        (run_dir / "data.obcs").write_text(DATA_OBCS)
        (run_dir / "data.diagnostics").write_text(DATA_DIAGNOSTICS)

    (RERUNS / "README.md").write_text(README)


if __name__ == "__main__":
    main()
