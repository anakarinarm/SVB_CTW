# CTW MITgcm reruns — Codex handoff

## Purpose

Prepare four controlled MITgcm reruns of the Sebastián Vizcaíno Bay (SVB) coastal-trapped-wave experiment using the **original MITgcm configuration as the physics/forcing baseline**, but replacing the bathymetry/partial-cell treatment and expanding diagnostics so that the reviewer questions can be answered quantitatively.

The scientific comparison must remain a **bathymetry-only experiment**. Do not modernize, retune, or otherwise change unrelated model physics unless a compile/runtime requirement makes it unavoidable and the change is explicitly documented.

---

# 1. Source-of-truth hierarchy

Use this order when deciding what to keep or change:

1. **Original MITgcm run files supplied by the user** (`data`, `data.pkg`, `data.obcs`, `data.diagnostics`, forcing files/scripts, `SIZE.h`, `packages.conf`, `CPP_OPTIONS.h`, etc.) are the source of truth for the original model physics and forcing.
2. The parameter changes in this README are the **approved rerun changes**.
3. `CTW_final_bathymetry_package/` is the source of truth for the new bathymetries.
4. The manuscript is useful as an audit of intent, but if manuscript text and the actual original configuration disagree, **do not silently choose one**. Report the discrepancy.

Do not reuse the historical bathymetry-generation notebooks or the historical 10×/1000× smoothing for the reruns. The new bathymetry has already been finalized.

---

# 2. Four production runs

Create all four runs from **one canonical copy of the original configuration**. Use the **same executable** for all four. Apart from `bathyFile`, run name/output directory, and any paths that necessarily differ, the four configurations should be identical.

| Short name | Geometry | Bathymetry file |
|---|---|---|
| `Rplus` | realistic + SVB | `bathy_realistic_plusSVB_FINAL.bin` |
| `Rminus` | realistic − SVB | `bathy_realistic_minusSVB_FINAL.bin` |
| `Splus` | smooth30 + SVB | `bathy_smooth30_plusSVB_FINAL.bin` |
| `Sminus` | smooth30 − SVB | `bathy_smooth30_minusSVB_FINAL.bin` |

These files are in `CTW_final_bathymetry_package/`.

The intended perturbation fields are

```text
delta_R(q) = q(Rplus) - q(Rminus)
delta_S(q) = q(Splus) - q(Sminus)
```

Do not compare only the raw `+SVB` runs when diagnosing the SVB-generated CTW.

### Bathymetry binary convention

The final binaries are:

- shape: `ny=612`, `nx=512`
- big-endian `float64`
- ocean depths **negative**
- land = `0`
- expected size of each file: `612*512*8 = 2,506,752 bytes`

SHA256 checksums:

```text
b9cd1468a972f6ce3b46b970a382bd79b5716acb27cd8279524e3a5a416ad0e1  bathy_realistic_minusSVB_FINAL.bin
5bd3f87847c6fa179953303795bf395a542dffc45f1e31b8249674e3e9293fa2  bathy_realistic_plusSVB_FINAL.bin
38fad368107da97a3cb1d0caca723040fe329cc417b41a1af984908447e5808e  bathy_smooth30_minusSVB_FINAL.bin
a55f58aea88a393c9bbf8d7c2b62fd923d8511156a853030733e1f938f888fed  bathy_smooth30_plusSVB_FINAL.bin
```

Ensure the original input precision/endian settings actually read these files as 64-bit big-endian. In a standard MITgcm setup this usually means checking `readBinaryPrec` and the compiler/runtime endian convention. Do not convert the files without documenting it.

---

# 3. Required partial-cell changes

These are approved changes and should be made in the rerun configuration.

## `data` / `PARM01`

```fortran
hFacMin   = 0.10,
hFacMinDr = 2.0,
```

## `data` / `PARM04`

```fortran
useMin4hFacEdges = .TRUE.,
pCellMix_select  = 0,
```

The final bathymetry already has a common minimum positive depth of **2 m**. Do not apply a second shallow-water clipping operation in the run scripts.

### Exact vertical grid

Keep the exact original 100-level vertical grid:

```fortran
delR =
 2, 2, 2, 2, 2, 2.2, 2.4, 2.6, 2.8, 3,
 3.2, 3.5, 3.8, 4.2, 4.5, 4.8, 5, 5, 5, 5,
 5, 5, 5.5, 5.5, 5.5, 5.5, 5.5, 6, 6.5, 7,
 7.5, 8, 8.5, 9.2, 10, 10, 10, 10, 10, 10,
 10, 11, 12, 13, 14, 15, 16, 17, 18, 20,
 22, 23, 25, 25, 25, 25, 25, 25, 25, 25,
 27, 28, 30, 30, 32, 34, 36, 38, 40, 43,
 46, 50, 55, 59, 65, 70, 75, 80, 85, 92,
 100, 100, 100, 100, 100, 100, 105, 115, 125, 132,
 145, 160, 170, 180, 190, 200, 210, 230, 245, 269.8,
```

This sums to exactly 4700 m.

### Why these settings are fixed

With `hFacMin=0.10`, `hFacMinDr=2 m`, and the 2-m minimum wet depth:

- no intended wet column becomes dry after hFac construction;
- minimum open bottom-cell thickness is 2 m;
- the final realistic/smooth pair masks remain identical within each source geometry;
- the tiny-island / one-cell-wet-bay topology audit remains clean.

The diagnostics are archived in:

```text
CTW_final_bathymetry_package/CTW_final_partial_cell_diagnostics.nc
CTW_final_bathymetry_package/final_partial_cell_summary.csv
CTW_final_bathymetry_package/final_topology_and_consistency_audit.csv
```

---

# 4. Important facts already baked into the new bathymetry

Do **not** reproduce these transformations in the MITgcm run scripts. They are already in the four final binary files.

1. Realistic and smooth members of a pair have the same land mask.
2. The immediate coastal depth is the same in realistic and smooth members of a pair.
3. The south open-boundary coastline is straightened through the OBC treatment zone to avoid a small awkward coastal bay at the boundary.
4. OBC bathymetry treatment uses a 15-cell flat normal-depth core plus a 15-cell smooth transition.
5. West, south, and north OBC bathymetry were treated; the northern land mask itself was left natural because that boundary is downstream of the wave.
6. The final straight southern coastline has last wet x-index `i=478` (Python zero-based) and longitude approximately `-114.5354004°`.
7. The no-SVB experiment preserves the continental slope seaward of the shelf break and only modifies the shelf/bay source geometry.
8. The no-SVB 250-m shelf-break spike near 28.6°N was deliberately regularized because following it literally created an awkward land spike.
9. The smooth experiment is the approved 30-km half-power **native-grid local-tangent** smoothing, not the historical repeated smoothing.

Useful build audit:

```text
grid: 512 x 612
west active OBC rows: 612
south active OBC columns: 479
north active OBC columns: 71
hFacMin: 0.10
hFacMinDr: 2 m
minimum positive depth: 2 m
```

Do not add another boundary flattening step after reading `bathyFile`.

---

# 5. OBCS checks Codex must make

The original physical boundary condition should remain **Orlanski radiation on north, south, and west**, with the eastern side closed, unless the actual original input files show otherwise.

Because the southern coastal mask was changed near the open boundary, inspect the original `data.obcs` / OBCS masks/files carefully:

- south: verify that boundary forcing/radiation is only active where the final bathymetry is wet;
- north: do not artificially straighten the coast; just verify that OBCS respects the final wet mask;
- west: all rows are wet through the treated boundary strip;
- do not prescribe boundary values into cells that are land in the new bathymetry.

If the old OBC implementation used hard-coded wet-index ranges derived from the old bathymetry, update those ranges from the **new wet mask**, not from the old constants.

---

# 6. Physics/forcing that should stay unchanged

The reviewer reruns are not the place to retune model physics. Start from the actual original configuration and preserve these values/options unless the original files prove they were different.

Known original/model-intent values to audit:

```text
deltaT = 60 s
run duration = 10 days = 864000 s = 14400 time steps
horizontal viscosity Ah ~= 100 m^2/s
horizontal tracer diffusivity Kh ~= 100 m^2/s
vertical viscosity Az ~= 1e-5 m^2/s
vertical tracer diffusivity Kz ~= 1e-5 m^2/s
```

Also preserve from the original configuration:

- momentum and tracer advection schemes;
- free-surface formulation;
- implicit/explicit viscosity/diffusion choices;
- equation of state;
- Coriolis treatment;
- bottom/lateral boundary conditions;
- initial T/S stratification files;
- initial velocities / pickup behavior;
- wind forcing implementation;
- OBCS numerical settings and relaxation/radiation parameters;
- all compile-time options unrelated to enabling the required diagnostics.

The manuscript intent was a horizontally uniform February-2008 stratification and a linear EOS with approximately

```text
alpha_T = 2.0e-4 1/degC
beta_S  = 7.4e-4 (in the original salinity convention)
```

but the **actual original namelist/files must win** if there is a discrepancy.

### Wind pulse

Preserve the exact original wind-forcing implementation. The documented target is:

```text
spatially homogeneous, upwelling-favorable alongshore stress
tau_max = 0.1 N/m^2
sigma = 1 hour
forcing rises to its maximum at t = 3.5 hr and then stops
orientation approximately -30 deg relative to north / toward the southeast
```

Do not rewrite the forcing formula from the manuscript if the original run already has a tested forcing file or routine. Reuse the original forcing exactly and verify its units/sign.

---

# 7. Compile-time/package policy

Build **one executable** and use it for all four production runs.

First inspect the original:

```text
packages.conf
CPP_OPTIONS.h
SIZE.h
DIAGNOSTICS_SIZE.h
data.pkg
```

### Required package capability

The rerun needs `pkg/diagnostics`. If the original build already has it, preserve the build and only modify `data.diagnostics` as needed.

If diagnostics is not compiled/enabled:

1. enable the diagnostics package in the same source tree;
2. set `useDiagnostics=.TRUE.` in `data.pkg`;
3. ensure `DIAGNOSTICS_SIZE.h` is large enough for the lists below;
4. recompile once;
5. use that executable for all four runs.

Do not change other packages while doing this.

### Pressure diagnostic compile check

After a short startup, inspect `available_diagnostics.log` and verify that at least these diagnostics exist:

```text
UVEL
VVEL
WVEL
THETA
SALT
RHOAnoma
PHIHYD
PHI_SURF
ETAN
```

The standard MITgcm diagnostics list defines `PHIHYD` as the hydrostatic pressure-potential anomaly and `PHI_SURF` as the surface dynamical pressure potential. If `PHIHYD` is absent because the corresponding calculation code was compiled out, enable the appropriate pressure-diagnostic code option and rebuild the **single common executable**.

Also look for these useful optional diagnostics:

```text
oceTAUX
oceTAUY
UVELPHI
VVELPHI
PHI_NH   # only relevant if the original configuration is nonhydrostatic
```

Do not assume a name is available merely because it exists in current MITgcm documentation; `available_diagnostics.log` from the actual executable is authoritative.

---

# 8. New diagnostics/output requirement

The old output cadence was 10 min for the first 5 days and 20 min for the last 5 days. For the reruns, use a **uniform 20-minute cadence for the entire 10-day integration** so all four runs can be differenced and Fourier-analyzed without cadence changes.

```text
output interval = 1200 s
```

Use positive `frequency=1200.` in `data.diagnostics` for time averages and negative `frequency=-1200.` for snapshots, following the diagnostics-package convention.

## A. Required 3-D time-averaged fields every 20 min

```text
UVEL
VVEL
WVEL
THETA
SALT
RHOAnoma
PHIHYD
```

Requested vertical coverage:

- all model levels 1–84 (level 84 center is about 1973 m);
- plus sparse deep-tail levels, suggested: 85, 88, 92, 96, 100.

This is intentionally generous because the rerun needs:

- perturbation kinetic and available-potential-energy diagnostics;
- complex harmonic structures near 0.84 cpd;
- mode-1/mode-2 structure comparison;
- pressure–velocity energy flux.

**Before production, Codex should estimate the disk volume of this diagnostic request using the actual I/O precision and number of tiles. Do not silently reduce the vertical levels.** If it is too large, report the estimate and let the user choose a reduced level set.

## B. Required 2-D fields every 20 min

Save at least:

```text
ETAN
PHI_SURF
```

Prefer **snapshots** every 20 min for ETAN (phase propagation / phase speed). A 20-min averaged ETAN stream can also be kept if convenient.

Save the actual applied wind stress if available:

```text
oceTAUX
oceTAUY
```

or the exact diagnostic names used by the forcing package in this build.

## C. Optional pressure-flux cross-check

If `UVELPHI` and `VVELPHI` are available and storage is acceptable, save them as a separate diagnostic stream for a consistency check. They are not a substitute for the perturbation-energy flux described below.

## D. Static grid fields

Make sure each run has, once, the geometry needed for offline integration:

```text
XC, YC, XG, YG
Depth
drF / delR
hFacC, hFacW, hFacS
dxG, dyG (or equivalent face metrics)
```

If the original NetCDF/MNC grid output already provides these, retain it.

---

# 9. Why pressure must be saved

The old `Fx` / `Fy` products were generated **offline**; they were not native MITgcm diagnostics. Do not try to reproduce them by looking for `Fx` or `Fy` fields in MITgcm.

The new analysis should reconstruct the SVB perturbation first:

```text
delta_u   = u_plus   - u_minus
delta_v   = v_plus   - v_minus
delta_w   = w_plus   - w_minus
delta_rho = rho_plus - rho_minus
delta_p   = p_plus   - p_minus
```

and then form perturbation wave-energy flux from the perturbation fields, schematically

```text
F = delta_p * delta_u_vector
```

Do **not** use

```text
F_plus - F_minus
```

as though it were equal to the flux of the perturbation; the operation is quadratic and those quantities are not equivalent.

Save both `PHIHYD` and `PHI_SURF` so the pressure convention can be reconstructed correctly offline. Do not blindly add `g*ETAN` to `PHIHYD`; first verify the exact MITgcm pressure-potential definitions in this build so that surface pressure is not double-counted.

---

# 10. Reviewer-driven analysis that the outputs must support

The reruns are specifically intended to separate several mechanisms that the old outputs could not distinguish cleanly.

The new output must support:

1. **Section-integrated alongshore energy flux** across the manuscript transects.
2. **Cross-section-integrated wave energy density** so that

   ```text
   c_E = F_section / E_section
   ```

   can be compared with Brink group velocity.
3. **Cross-shore redistribution** of flux without confusing it with alongshore transport.
4. **Mode-character changes near 0.84 cpd**, especially around 31.69°N and 32.39°N, where new Brink mode-1/mode-2 eigenfunctions are available at the actual CTW frequency.
5. Comparison of realistic vs smooth runs without the old mask/coastal-depth confounder.

For mode analysis, use the complex harmonic coefficient of the perturbation at approximately 0.84 cpd and compare/joint-fit it to the local Brink modes. Do not classify modes solely from peak/trough snapshots.

---

# 11. Known analysis/code pitfalls — do not repeat them

These were found while auditing the original runs and postprocessing.

## A. Old realistic/smooth masks were not actually identical

The old 10× and 1000× realized grids differed in hundreds of wet/dry cells. The new bathymetries explicitly fix this. Do not regenerate masks from the old bathymetry or old smoother.

## B. Old coastal depths differed between smoothing cases

The new R/S members have the same coast/mask within a source geometry. Verify the generated `hFacC` wet-column mask after model startup before launching production.

## C. Historical OBC bathymetry code did not match the manuscript description

The old bathymetry notebook flattened west + south but not north even though the manuscript described north/south/west open boundaries. The approved final bathymetries already contain the corrected common OBC bathymetry treatment. Do not reapply the historical loop.

## D. Historical no-SVB documentation was inconsistent

The old no-SVB construction actually used the 250-m shelf-break logic even though one manuscript description referred to 150 m. The final no-SVB bathymetry is already approved and should not be reconstructed inside the run setup.

## E. Do not reuse the old repeated smoother

The historical smoother had both scientific/provenance problems and an implementation with directional bias. The final 30-km local-tangent bathymetry is already generated.

## F. Figure-7 flux-vector projection bug

In the old `EFLUX_fig7.py`, a geographic-distance helper expected `(lon1,lat1,lon2,lat2)` but was called with latitude/longitude swapped. The resulting blue “normal” projection was not a trustworthy cross-shore flux. Do not copy that geometry code into new analysis.

Use transect tangent/normal vectors derived consistently from the actual section endpoints/grid geometry.

## G. Old “period” windows were actually 86400-s daily windows

The CTW frequency is about 0.84 cpd, so the physical period is about 1.19 days (~102,900 s), not 86400 s. Label old 86400-s blocks as daily windows, not CTW periods. For the reruns, prefer harmonic analysis at the observed frequency.

---

# 12. Run reproducibility / pairing requirements

All four runs must have:

- the same executable;
- the same start time and end time;
- the same `deltaT`;
- the same diagnostic timestamps;
- the same forcing files and forcing timing;
- the same initial stratification;
- the same package settings;
- the same output precision;
- the same OBC numerical options;
- the same viscosity/diffusivity/advection settings.

Only the approved bathymetry differs.

Within `Rplus/Rminus` and `Splus/Sminus`, exact timestamp alignment is essential because perturbation fields are formed by direct subtraction.

---

# 13. Pre-production validation Codex should automate

Before submitting the four 10-day jobs, create a small validation script and perform these checks.

## Configuration diff

Generate a machine-readable diff of the four input directories. The only scientifically intentional differences should be:

```text
bathyFile
run/output names and paths
```

plus unavoidable bookkeeping.

## Bathymetry file check

For each run:

- file exists;
- size = 2,506,752 bytes;
- SHA256 matches this README;
- startup log reports sensible depths and no binary/endian failure.

## hFac/grid check

After initialization, verify:

```text
Rplus mask == Splus mask
Rminus mask == Sminus mask
```

and compare resulting `Depth` / wet-column counts with the supplied audit:

```text
Rplus/Splus wet columns  : 211096
Rminus/Sminus wet columns: 204908
```

There should be no hFac-induced input-mask changes.

## Short smoke tests

Before production, run a short integration for each distinct geometry (at least Rplus and Splus; ideally all four):

- check for NaNs/Infs;
- check `max|U|`, `max|V|`, `max|W|`, `max|ETAN|`;
- inspect solver convergence;
- inspect OBC behavior;
- confirm diagnostic files are produced at exactly 1200-s cadence;
- confirm required diagnostic names and vertical levels are present.

A short zero-wind/rest test is useful for detecting bathymetry/OBC startup currents before testing the wind pulse. Do not replace the production forcing with this; it is only a numerical smoke test.

---

# 14. Suggested directory layout

```text
ctw_reruns/
├── build/                     # one common executable
├── common_input/              # original forcing / initial T,S / OBC files
├── bathymetry/
│   ├── bathy_realistic_plusSVB_FINAL.bin
│   ├── bathy_realistic_minusSVB_FINAL.bin
│   ├── bathy_smooth30_plusSVB_FINAL.bin
│   └── bathy_smooth30_minusSVB_FINAL.bin
├── Rplus/
├── Rminus/
├── Splus/
├── Sminus/
└── scripts/
    ├── validate_configs.py
    ├── validate_bathy.py
    ├── estimate_diagnostics_storage.py
    └── submit_all.sh
```

Prefer symlinks from run directories to immutable common forcing/initial files and the bathymetry directory so that the four experiments cannot silently drift apart.

---

# 15. Codex task sequence

Codex should proceed in this order:

1. Inspect the original MITgcm configuration and summarize all runtime/compile settings before changing anything.
2. Compare actual original settings with the expected values in this README and flag discrepancies.
3. Make one clean rerun branch/directory; do not edit the only copy of the original config.
4. Enable/resize diagnostics only if required.
5. Insert the approved partial-cell settings and exact `delR`.
6. Create the four run directories from one canonical template.
7. Set only the corresponding `bathyFile` in each.
8. Audit OBC masks/ranges against the new bathymetry, especially the straightened southern coast.
9. Configure the new 20-min diagnostic streams.
10. Estimate output volume before production.
11. Compile one executable, if recompilation is necessary.
12. Run initialization/smoke tests and inspect `available_diagnostics.log`.
13. Produce configuration/bathymetry/hFac validation reports.
14. Only then prepare the four 10-day production submissions.

At every step, prefer an explicit diff over an implicit assumption.

---

# 16. Files to give Codex alongside this README

At minimum:

```text
CTW_final_bathymetry_package/
    bathy_realistic_plusSVB_FINAL.bin
    bathy_realistic_minusSVB_FINAL.bin
    bathy_smooth30_plusSVB_FINAL.bin
    bathy_smooth30_minusSVB_FINAL.bin
    CTW_final_bathymetries.nc
    CTW_final_partial_cell_diagnostics.nc
    FINAL_data_partial_cells_and_delR.txt
    final_build_summary.csv
    final_partial_cell_summary.csv
    final_topology_and_consistency_audit.csv

CTW_build_final_bathymetries.py   # provenance/reproduction, not needed at runtime
```

Also give Codex the **original MITgcm configuration tree** and all original forcing/initial-condition files.

Optional scientific-context files useful for postprocessing planning:

```text
CTW_mode1_mode2_084cpd/CTW_mode1_mode2_projection_template.ipynb
CTW_new_bathymetry_spectra_review.ipynb
CTW_extract_cross_section_integrated_flux.ipynb
CTW_T5_T6_Brink_detailed_analysis.ipynb
```

---

# 17. Final non-negotiable principle

The reruns are designed to answer whether short-scale bathymetric complexity modifies CTW propagation and energetics. Therefore **do not change the dynamical model at the same time as the bathymetry**.

The clean experiment is:

```text
same code
same forcing
same stratification
same numerics
same diagnostics
same timestamps
new controlled bathymetry only
```

with the approved partial-cell settings applied identically to all four runs.
