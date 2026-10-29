# CTW/SVB South-Extended Telescoped Reruns

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
