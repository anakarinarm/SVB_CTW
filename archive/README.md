# Archive: earlier CTW_SVB configurations

> Written by Claude, approved by KRM.

Moved here from the top level of `CTW_SVB` on 2026-10-29. Each stage replaced
the previous one; none of them produced production output (only smoke tests).
The final configuration is `../ctw_reruns_telescoped_west64_south72`.

| Group | Domain | What changed | Outcome |
|---|---|---|---|
| `1_original_512x612/` | `512 x 612 x 100` | `ctw_reruns/`: first reruns from `real_SVB` with the approved bathymetries, partial-cell settings and new diagnostics (MNC NetCDF, exch2); `code/` | Boundaries sat on shelf/slope topography and acted as a continuation of the waveguide |
| `2_telescoped_buffer/` | `512 x 612 x 100` | `ctw_reruns_telescoped/`: W/S/N strips replaced by a flat 4000 m telescoped buffer, wind tapered there, MDS output, no exch2; `code_telescoped*/` | The 48-cell sponge still overlapped the original domain |
| `3_telescoped_south/` | `512 x 648 x 100` | `ctw_reruns_telescoped_south/`: 36 flat rows added to the south; `code_telescoped_south_noexch2/` | Replaced by the west64/south72 extension |

Each group's `scripts/` holds the setup, input and QC scripts used only by
that stage. The links inside the run folders (for example
`ctw_reruns_telescoped/setup_ctw_telescoped_reruns.py`) point to them.

## What was kept

- Namelists, READMEs, validation files, figures, `STDOUT.0000` and
  `available_diagnostics.log` of the smoke tests.
- The smoke-test model output, the build folders and the stage-specific
  input binaries (`../input/telescoped_buffer/*.bin`,
  `../input/telescoped_south_extended/*.bin`) were approved for deletion on
  2026-10-29. The scripts in each group's `scripts/` can regenerate those
  inputs.
- `ctw_reruns/common_input` is now a link to
  `../../../input/common_input_512x612`. The final inputs are built from those
  files, so they stay in `input/`.

## Rerunning an archived stage

- **Paths.** All links were rewritten when the folders moved, so they resolve
  as before. Each group has an `input -> ../../input` link.
- **Scripts.** The scripts still compute paths relative to their old location
  (for example `../real_SVB` or `ctw_reruns/`). Check them before rerunning.
- **Executables.** Rebuild with `genmake2 -mods=<group>/code_*` (see the
  `Makefile` header of the final `build_binary` for the options), into the
  build folder that the run's `mitgcmuv` link expects.
