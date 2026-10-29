#!/usr/bin/env python3
"""Create CTW/SVB rerun directories with added west and south flat buffers."""

from __future__ import annotations

import shutil
from pathlib import Path

from setup_ctw_reruns import DATA_DIAGNOSTICS, DATA_PKG, link_force


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "real_SVB"
RERUNS = ROOT / "ctw_reruns_telescoped_west64_south72"
INPUTS = ROOT / "input" / "telescoped_west64_south72"

RUNS = {
    "Rplus": "bathy_realistic_plusSVB_W64S72.bin",
    "Rminus": "bathy_realistic_minusSVB_W64S72.bin",
    "Splus": "bathy_smooth30_plusSVB_W64S72.bin",
    "Sminus": "bathy_smooth30_minusSVB_W64S72.bin",
    "S50plus": "bathy_smooth50_plusSVB_W64S72.bin",
    "S50minus": "bathy_smooth50_minusSVB_W64S72.bin",
    "S70plus": "bathy_smooth70_plusSVB_W64S72.bin",
    "S70minus": "bathy_smooth70_minusSVB_W64S72.bin",
}

DATA_PKG_BINARY = DATA_PKG.replace(" useMNC=.TRUE.,\n", " useMNC=.FALSE.,\n")
DATA_DIAGNOSTICS_BINARY = DATA_DIAGNOSTICS.replace("  diag_mnc = .TRUE.,", "  diag_mnc = .FALSE.,")


def read_summary() -> dict[str, str]:
    out = {}
    for line in (INPUTS / "telescoped_west64_south72_summary.csv").read_text().splitlines()[1:]:
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
 delXFile='../west_south_input/delX_W64S72.bin',
 delYFile='../west_south_input/delY_W64S72.bin',
 useMin4hFacEdges = .TRUE.,
 pCellMix_select  = 0,
 &

 &PARM05
 bathyFile='../west_south_input/{bathy}',
 meridWindFile='../west_south_input/upwelling_merid_halfgauss2hr_16day_576x684_W64S72.bin',
 zonalWindFile='../west_south_input/upwelling_zonal_halfgauss2hr_16day_576x684_W64S72.bin',
 hydrogThetaFile='../west_south_input/feb_Tref_100zlev_576x684_W64S72.bin',
 hydrogSaltFile='../west_south_input/feb_Sref_100zlev_576x684_W64S72.bin',
 &
"""


DATA_OBCS = """# Open-boundaries
 &OBCS_PARM01
 OB_Iwest = 684*1,
 OB_Jsouth = 543*1,33*0,
 OB_Jnorth = 137*-1,439*0,
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
 OBCSsponge_Salt=.FALSE.,
 OBCSsponge_Theta=.FALSE.,
 useOBCSbalance=.TRUE.,
 useOBCSprescribe=.TRUE.,
 OBWtFile='../west_south_input/obcs_ts_from_initial_hydrography_16day/OBWt.bin',
 OBStFile='../west_south_input/obcs_ts_from_initial_hydrography_16day/OBSt.bin',
 OBNtFile='../west_south_input/obcs_ts_from_initial_hydrography_16day/OBNt.bin',
 OBWsFile='../west_south_input/obcs_ts_from_initial_hydrography_16day/OBWs.bin',
 OBSsFile='../west_south_input/obcs_ts_from_initial_hydrography_16day/OBSs.bin',
 OBNsFile='../west_south_input/obcs_ts_from_initial_hydrography_16day/OBNs.bin',
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


README = """# CTW/SVB West64/South72 Telescoped Reruns

This run family increases the domain from 512 x 612 x 100 to 576 x 684 x 100
by adding 64 western columns and 72 southern rows. The original bay/interior
grid is preserved and shifted away from both open boundaries.

Key choices:

- Added western columns and southern rows are flat 4400 m abyssal buffer at
  the open boundaries, with straight artificial coastlines, and are smoothly
  blended into the active-domain edge depths to avoid sponge-corner steps.
- The northern buffer/taper band locally continues existing wet-cell depths
  from the row immediately south of the band while preserving the land mask and
  north OBC wet span.
- The southwest added-buffer corner has a local 2-D wet-cell smoothing pass;
  the northwest added-buffer corner has a smoothed along-boundary profile with
  a gradual northward transition.
- The 48-cell OBCS sponge is entirely outside the original domain on the west
  and south sides; SVB stays outside the sponge and wind-taper regions.
- Added-grid spacing is telescoped by the same factor in X and Y so cell aspect
  ratios stay close to the original near-square physical grid.
- South OBC wet span is `543*1,33*0`; north OBC wet span is `137*-1,439*0`;
  west OBC is wet along all 684 rows.
- Boundary T/S are prescribed from the initial horizontally uniform profiles,
  but T/S sponge relaxation is disabled.
- Diagnostics are MITgcm MDS binary `.data/.meta`, not MNC NetCDF.
- Exch2 and MNC are not used.
- Additional smooth50 sensitivity runs are available as `S50plus` and
  `S50minus`; these are identical to `Splus` and `Sminus` except for the
  higher-smoothing bathymetry files
  `bathy_smooth50_plusSVB_W64S72.bin` and
  `bathy_smooth50_minusSVB_W64S72.bin`, generated from the original
  512 x 612 `*_FINAL.bin` source files.
- Additional smooth70 sensitivity runs are available as `S70plus` and
  `S70minus`; these were generated locally with
  `CTW_build_final_bathymetries.py --smooth-length-km 70 --smooth-label smooth70`
  and then extended into the same west64/south72 domain.
- A 20-minute Rplus MPI smoke test with 16 ranks passed with `STOP NORMAL END`,
  OBCS checks OK, T/S sponge off, and binary diagnostics written.

Because Nx and Ny changed, use the matching executable built from
`code_telescoped_west64_south72_noexch2`.
"""


def main() -> None:
    summary = read_summary()
    RERUNS.mkdir(exist_ok=True)
    (RERUNS / "west_south_input").mkdir(exist_ok=True)
    (RERUNS / "build").mkdir(exist_ok=True)
    for path in INPUTS.glob("*"):
        if path.is_file() or path.is_dir():
            link_force(RERUNS / "west_south_input" / path.name, path)
    for name in [
        "make_telescoped_west64_south72_inputs.py",
        "setup_ctw_telescoped_west64_south72_reruns.py",
        "prepare_external_output_dirs.py",
        "run_west64_south72_s50_15d_sequential.sh",
        "run_west64_south72_s70_15d_sequential.sh",
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
        (run_dir / "data.obcs").write_text(DATA_OBCS)
        (run_dir / "data").write_text(
            data_template(summary["xgOrigin_data_value"], summary["ygOrigin"], bathy)
        )
    (RERUNS / "README.md").write_text(README)
    print(f"wrote west64/south72 telescoped reruns to {RERUNS}")


if __name__ == "__main__":
    main()
