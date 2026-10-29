#!/usr/bin/env bash
set -euo pipefail
set -o noclobber
ROOT="${1:-/Volumes/LaCie/CTW_SVB_reruns_telescoped_west64_south72}"
NPROC="${NPROC:-16}"
[[ "$NPROC" == 16 ]] || { echo 'Executable requires exactly 16 MPI ranks'; exit 1; }
for run in Rplus_G5d Rminus_G5d; do
  dir="${ROOT}/${run}"
  [[ -x "${dir}/mitgcmuv" ]] || { echo "Missing executable in ${dir}"; exit 1; }
  grep -q 'endTime=2592000' "${dir}/data"
  grep -q 'externForcingCycle=2678400' "${dir}/data"
  if [[ -n "$(find "${dir}" -maxdepth 1 -type f \( -name '*.data' -o -name '*.meta' -o -name 'STDOUT.*' -o -name 'STDERR.*' -o -name '*.mpirun.log' \) -print -quit)" ]]; then
    echo "Refusing to overwrite outputs in ${dir}"
    exit 1
  fi
done
if ! df -Pk "$ROOT" | awk 'NR == 2 { if ($4 * 1024 < 5.3e12) exit 1 }'; then
  echo 'The pair needs approximately 5.3 TB free with the unchanged production output settings'
  exit 1
fi
for run in Rplus_G5d Rminus_G5d; do
  echo "[$(date)] Starting ${run}, 30 days, 16 MPI ranks"
  (cd "${ROOT}/${run}"; mpirun -np "$NPROC" ./mitgcmuv) >"${ROOT}/${run}/${run}.mpirun.log" 2>&1
  grep -q 'STOP NORMAL END' "${ROOT}/${run}/STDOUT.0000"
  echo "[$(date)] Finished ${run}"
done
