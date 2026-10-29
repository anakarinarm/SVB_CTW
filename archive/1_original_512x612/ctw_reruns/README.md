# CTW/SVB production reruns

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
