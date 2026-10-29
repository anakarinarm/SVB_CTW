#!/usr/bin/env python3
"""Create CTW/SVB rerun directories with a south-extended flat buffer."""

from __future__ import annotations

import shutil
from pathlib import Path

from setup_ctw_reruns import DATA_DIAGNOSTICS, DATA_PKG, link_force


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "real_SVB"
RERUNS = ROOT / "ctw_reruns_telescoped_south"
SOUTHEXT = ROOT / "input" / "telescoped_south_extended"

RUNS = {
    "Rplus": "bathy_realistic_plusSVB_SOUTHEXT.bin",
    "Rminus": "bathy_realistic_minusSVB_SOUTHEXT.bin",
    "Splus": "bathy_smooth30_plusSVB_SOUTHEXT.bin",
    "Sminus": "bathy_smooth30_minusSVB_SOUTHEXT.bin",
}

DATA_PKG_BINARY = DATA_PKG.replace(" useMNC=.TRUE.,\n", " useMNC=.FALSE.,\n")
DATA_DIAGNOSTICS_BINARY = DATA_DIAGNOSTICS.replace("  diag_mnc = .TRUE.,", "  diag_mnc = .FALSE.,")


def read_summary() -> dict[str, str]:
    out = {}
    for line in (SOUTHEXT / "telescoped_south_extended_summary.csv").read_text().splitlines()[1:]:
        key, value = line.split(",", 1)
        out[key] = value
    return out


def data_template(x_origin: str, y_origin: str, bathy: str) -> str:
    return f"""# Model parameters
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
 xgOrigin={x_origin},
 ygOrigin={y_origin},
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
 delXFile='../south_extended_input/delX_SOUTHEXT.bin',
 delYFile='../south_extended_input/delY_SOUTHEXT.bin',
 useMin4hFacEdges = .TRUE.,
 pCellMix_select  = 0,
 &

 &PARM05
 bathyFile='../south_extended_input/{bathy}',
 meridWindFile='../south_extended_input/upwelling_merid_halfgauss2hr_16day_512x648_SOUTHEXT.bin',
 zonalWindFile='../south_extended_input/upwelling_zonal_halfgauss2hr_16day_512x648_SOUTHEXT.bin',
 hydrogThetaFile='../south_extended_input/feb_Tref_100zlev_512x648_SOUTHEXT.bin',
 hydrogSaltFile='../south_extended_input/feb_Sref_100zlev_512x648_SOUTHEXT.bin',
 &
"""


DATA_OBCS_SOUTHEXT = """# Open-boundaries
 &OBCS_PARM01
 OB_Iwest = 648*1,
 OB_Jsouth = 479*1,33*0,
 OB_Jnorth = 73*-1,439*0,
 useOrlanskiWest=.TRUE.,
 useOrlanskiSouth=.TRUE.,
 useOrlanskiNorth=.TRUE.,
 useOBCSsponge=.TRUE.,
 useLinearSponge=.TRUE.,
 OBCSsponge_W=.TRUE.,
 OBCSsponge_S=.TRUE.,
 OBCSsponge_N=.TRUE.,
 OBCSsponge_UatEW=.TRUE.,
 OBCSsponge_VatEW=.TRUE.,
 OBCSsponge_UatNS=.TRUE.,
 OBCSsponge_VatNS=.TRUE.,
 OBCSsponge_Salt=.TRUE.,
 OBCSsponge_Theta=.TRUE.,
 useOBCSbalance=.TRUE.,
 useOBCSprescribe=.TRUE.,
 OBWtFile='../south_extended_input/obcs_ts_from_initial_hydrography_16day/OBWt.bin',
 OBStFile='../south_extended_input/obcs_ts_from_initial_hydrography_16day/OBSt.bin',
 OBNtFile='../south_extended_input/obcs_ts_from_initial_hydrography_16day/OBNt.bin',
 OBWsFile='../south_extended_input/obcs_ts_from_initial_hydrography_16day/OBWs.bin',
 OBSsFile='../south_extended_input/obcs_ts_from_initial_hydrography_16day/OBSs.bin',
 OBNsFile='../south_extended_input/obcs_ts_from_initial_hydrography_16day/OBNs.bin',
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
 spongeThickness=48,
 Urelaxobcsinner=43200.E0,
 Vrelaxobcsinner=43200.E0,
 Urelaxobcsbound=1800.E0,
 Vrelaxobcsbound=1800.E0,
 &
"""


README = """# CTW/SVB South-Extended Telescoped Reruns

This run family increases the domain from 512 x 612 x 100 to 512 x 648 x 100
by prepending 36 flat abyssal rows south of the previous southern boundary.
The original bay/interior rows are preserved; only the artificial flat buffer
is extended.

Key choices:

- 36 added southern rows, about 0.49 degrees or 54 km.
- Added rows are 4000 m deep over the same open-water columns as the previous
  south OBC and land over the remaining eastern columns.
- South wind taper is confined to the added flat rows; the original domain is
  not south-tapered, so wind stress over SVB is preserved.
- West and north flat-buffer/taper choices are retained from the telescoped
  setup.
- OBCS sponge remains active, linear, and 48 cells thick for this first
  south-extended test.
- Diagnostics are MITgcm MDS binary `.data/.meta`, not MNC NetCDF.
- Exch2 and MNC are not used.
- A 20-minute Rplus MPI smoke test with 16 ranks passed with `STOP NORMAL END`,
  OBCS checks OK, and binary diagnostics written.

Because Ny changed, use the matching executable built from
`code_telescoped_south_noexch2`.
"""


def main() -> None:
    summary = read_summary()
    RERUNS.mkdir(exist_ok=True)
    (RERUNS / "south_extended_input").mkdir(exist_ok=True)
    (RERUNS / "build").mkdir(exist_ok=True)
    for path in SOUTHEXT.glob("*"):
        if path.is_file() or path.is_dir():
            link_force(RERUNS / "south_extended_input" / path.name, path)
    for name in [
        "make_telescoped_south_extended_inputs.py",
        "setup_ctw_telescoped_south_extended_reruns.py",
        "prepare_external_output_dirs.py",
    ]:
        link_force(RERUNS / name, ROOT / "scripts" / name)

    for run, bathy in RUNS.items():
        run_dir = RERUNS / run
        run_dir.mkdir(exist_ok=True)
        shutil.copy2(SOURCE / "eedata", run_dir / "eedata")
        (run_dir / "data.pkg").write_text(DATA_PKG_BINARY)
        for stale in ["data.mnc", "data.exch2"]:
            path = run_dir / stale
            if path.exists() or path.is_symlink():
                path.unlink()
        (run_dir / "data.diagnostics").write_text(DATA_DIAGNOSTICS_BINARY)
        (run_dir / "data.obcs").write_text(DATA_OBCS_SOUTHEXT)
        (run_dir / "data").write_text(
            data_template(summary["xgOrigin_data_value"], summary["ygOrigin"], bathy)
        )
    (RERUNS / "README.md").write_text(README)
    print(f"wrote south-extended telescoped reruns to {RERUNS}")


if __name__ == "__main__":
    main()
