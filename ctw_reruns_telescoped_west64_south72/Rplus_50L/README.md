# Rplus_50L

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

The full run was launched in the detached screen session CTW_Rplus_50L after
space was freed on LaCie (approximately 2.7 TiB available at launch). The
five-timestep test in smoke_Rplus_50L_hfac completed normally before launch.
The production configuration remains unchanged.

Reattach from a terminal:

```bash
screen -r CTW_Rplus_50L
```

Detach again with Ctrl-a followed by d. Output progress is recorded in
STDOUT.0000 and the MPI launcher log is Rplus_50L.mpirun.log in the external
execution directory. run_screen.sh records the launch command.
