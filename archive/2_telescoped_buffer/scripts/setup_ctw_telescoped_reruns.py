#!/usr/bin/env python3
"""Create CTW/SVB rerun directories using telescoped flat abyssal buffers."""

from __future__ import annotations

import shutil
from pathlib import Path

from setup_ctw_reruns import DATA_DIAGNOSTICS, DATA_PKG, LEVEL_LINES, link_force


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "real_SVB"
BASE_RERUNS = ROOT / "ctw_reruns"
RERUNS = ROOT / "ctw_reruns_telescoped"
TELEBUF = ROOT / "input" / "telescoped_buffer"

RUNS = {
    "Rplus": "bathy_realistic_plusSVB_TELEBUF.bin",
    "Rminus": "bathy_realistic_minusSVB_TELEBUF.bin",
    "Splus": "bathy_smooth30_plusSVB_TELEBUF.bin",
    "Sminus": "bathy_smooth30_minusSVB_TELEBUF.bin",
}

DATA_PKG_BINARY = DATA_PKG.replace(" useMNC=.TRUE.,\n", " useMNC=.FALSE.,\n")
DATA_DIAGNOSTICS_BINARY = DATA_DIAGNOSTICS.replace("  diag_mnc = .TRUE.,", "  diag_mnc = .FALSE.,")


def read_summary() -> dict[str, str]:
    out = {}
    for line in (TELEBUF / "telescoped_buffer_summary.csv").read_text().splitlines()[1:]:
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
 delXFile='../telescoped_input/delX_TELEBUF.bin',
 delYFile='../telescoped_input/delY_TELEBUF.bin',
 useMin4hFacEdges = .TRUE.,
 pCellMix_select  = 0,
 &

 &PARM05
 bathyFile='../telescoped_input/{bathy}',
 meridWindFile='../telescoped_input/upwelling_merid_halfgauss2hr_16day_512x612_TELEBUF.bin',
 zonalWindFile='../telescoped_input/upwelling_zonal_halfgauss2hr_16day_512x612_TELEBUF.bin',
 hydrogThetaFile='../common_input/feb_Tref_100zlev_512x612.bin',
 hydrogSaltFile='../common_input/feb_Sref_100zlev_512x612.bin',
 &
"""


DATA_OBCS_TELE = """# Open-boundaries
 &OBCS_PARM01
 OB_Iwest = 612*1,
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
 spongeThickness=48,
 Urelaxobcsinner=43200.E0,
 Vrelaxobcsinner=43200.E0,
 Urelaxobcsbound=1800.E0,
 Vrelaxobcsbound=1800.E0,
 &
"""

README = """# CTW/SVB Telescoped-Buffer Reruns

This is a second-generation run family with the same 512 x 612 x 100 model
size and the same 16-rank executable, but with west/south/north boundary strips
replaced by a flat abyssal buffer.

Status note:

- The original all-wet rectangular buffer failed `OBCS_CHECK`. This version
  preserves the final wet/dry topology while flattening the wet sponge-region
  bathymetry and removes exch2 from the telescoped build.
- A 20-minute Rplus MPI smoke test with 16 ranks passed with `STOP NORMAL END`
  using `build_binary/mitgcmuv` and wrote MDS binary diagnostic files.
- Diagnostics are written as MITgcm MDS binary `.data/.meta` output, not MNC
  NetCDF, to speed up writing and postprocessing.

Key changes relative to `../ctw_reruns`:

- `delXFile` and `delYFile` define a telescoped grid.
- The old western 48 cells, southern 40 cells, and northern 40 cells are
  sacrificial buffer zones.
- Bathymetry in wet parts of the buffer is a flat 4000 m abyssal plane with a
  30-cell transition into the final approved bathymetry. Original land cells remain
  land so the OBCS wet/dry topology stays compatible with MITgcm.
- Wind stress is tapered to zero over the same buffer/transition mask.
- The OBC active spans preserve the original boundary topology: west all rows,
  south columns 1-479, north columns 1-73.
- Exch2 and MNC are not used by this run family.
- OBCS sponge is active and linear with 48-cell thickness.

The goal is to keep the artificial open boundaries off dynamically meaningful
shelf/slope topography so they stop behaving like continuations of the coastal
waveguide.
"""


def main() -> None:
    summary = read_summary()
    RERUNS.mkdir(exist_ok=True)
    (RERUNS / "telescoped_input").mkdir(exist_ok=True)
    (RERUNS / "build").mkdir(exist_ok=True)
    link_force(RERUNS / "common_input", BASE_RERUNS / "common_input")
    for path in TELEBUF.glob("*"):
        if path.is_file():
            link_force(RERUNS / "telescoped_input" / path.name, path)
    for name in [
        "make_telescoped_buffer_inputs.py",
        "setup_ctw_telescoped_reruns.py",
        "prepare_external_output_dirs.py",
    ]:
        link_force(RERUNS / name, ROOT / "scripts" / name)

    for run, bathy in RUNS.items():
        run_dir = RERUNS / run
        run_dir.mkdir(exist_ok=True)
        shutil.copy2(SOURCE / "eedata", run_dir / "eedata")
        (run_dir / "data.pkg").write_text(DATA_PKG_BINARY)
        mnc = run_dir / "data.mnc"
        if mnc.exists() or mnc.is_symlink():
            mnc.unlink()
        exch2 = run_dir / "data.exch2"
        if exch2.exists() or exch2.is_symlink():
            exch2.unlink()
        (run_dir / "data.diagnostics").write_text(DATA_DIAGNOSTICS_BINARY)
        (run_dir / "data.obcs").write_text(DATA_OBCS_TELE)
        (run_dir / "data").write_text(
            data_template(summary["xgOrigin_data_value"], summary["ygOrigin"], bathy)
        )
    (RERUNS / "README.md").write_text(README)
    print(f"wrote telescoped reruns to {RERUNS}")


if __name__ == "__main__":
    main()
