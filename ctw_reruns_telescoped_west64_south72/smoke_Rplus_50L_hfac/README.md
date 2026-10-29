# Rplus_50L hFac Smoke Test

This isolated copy runs five 60-s timesteps. Only endTime (300 s), monitorFreq
(60 s), dumpFreq (0) and chkptFreq (0) differ from Rplus_50L. Grid output,
including hFacC/hFacW/hFacS, is written on LaCie. The full production case
remains unlaunched. The configuration description below refers to the template.

Smaller-SVB experiment copied from the current production Rplus configuration.
The only namelist change is:

```diff
- bathyFile='../west_south_input/bathy_realistic_plusSVB_W64S72.bin',
+ bathyFile='../west_south_input/bathy_realistic_SVB50L_mouth_W64S72.bin',
```

data.obcs, data.diagnostics, data.pkg and eedata are byte-identical to Rplus.
The grid, stratification, initial conditions, boundary conditions, sponge,
mixing and wind forcing are unchanged. The run uses the existing production
executable, exactly 16 MPI ranks, a 60-s timestep, 1200-s diagnostic output,
and the original 15-day duration (endTime=1296000 s).

The supplied bathymetry is used as-is. The no-bay counterpart is the existing
Rminus run; no additional minus case is required.

External execution directory:
/Volumes/LaCie/CTW_SVB_reruns_telescoped_west64_south72/Rplus_50L

Prepared only; no run has been launched. On 2026-10-13 LaCie had approximately
327 GiB free, insufficient for this case's approximately 1.3 TB of output.
Make sufficient space available before launching. Production configuration
and output files were not modified.
