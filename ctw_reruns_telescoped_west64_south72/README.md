# CTW/SVB West64/South72 Telescoped Reruns

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
