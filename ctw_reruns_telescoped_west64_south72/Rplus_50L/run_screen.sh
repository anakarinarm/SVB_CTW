#!/usr/bin/env bash
set -euo pipefail
set -o noclobber
cd /Volumes/LaCie/CTW_SVB_reruns_telescoped_west64_south72/Rplus_50L
echo "Starting Rplus_50L: 15 days, 16 MPI ranks"
echo "Progress: $(pwd)/STDOUT.0000"
/opt/homebrew/bin/mpirun -np 16 ./mitgcmuv > Rplus_50L.mpirun.log 2>&1
grep -q 'STOP NORMAL END\|Execution ended Normally' STDOUT.0000
echo "Rplus_50L completed normally"
