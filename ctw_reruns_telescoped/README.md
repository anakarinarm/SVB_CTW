# CTW/SVB Telescoped-Buffer Reruns

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
