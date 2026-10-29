#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/Volumes/LaCie/CTW_SVB_reruns_telescoped_west64_south72}"
NPROC="${NPROC:-16}"
RUNS=(S70plus S70minus)
STAMP="$(date +%Y%m%d_%H%M%S)"
LOGDIR="${ROOT}/batch_logs_s70_${STAMP}"

mkdir -p "${LOGDIR}"

clean_outputs() {
  local run_dir="$1"
  find "${run_dir}" -maxdepth 1 -type f \( \
    -name '*.data' -o \
    -name '*.meta' -o \
    -name 'STDOUT.*' -o \
    -name 'STDERR.*' -o \
    -name 'monitor.*' -o \
    -name 'pickup*' -o \
    -name 'available_diagnostics.log' -o \
    -name 'output.txt' \
  \) -delete
}

for run in "${RUNS[@]}"; do
  run_dir="${ROOT}/${run}"
  log="${LOGDIR}/${run}.mpirun.log"

  if [[ ! -d "${run_dir}" ]]; then
    echo "missing run directory: ${run_dir}" | tee -a "${LOGDIR}/batch.log"
    exit 1
  fi
  if [[ ! -x "${run_dir}/mitgcmuv" ]]; then
    echo "missing executable: ${run_dir}/mitgcmuv" | tee -a "${LOGDIR}/batch.log"
    exit 1
  fi
  if ! grep -q "endTime=1296000" "${run_dir}/data"; then
    echo "unexpected endTime in ${run_dir}/data; expected 1296000 seconds" | tee -a "${LOGDIR}/batch.log"
    exit 1
  fi

  echo "[$(date)] cleaning generated outputs for ${run}" | tee -a "${LOGDIR}/batch.log"
  clean_outputs "${run_dir}"

  echo "[$(date)] starting ${run} with ${NPROC} MPI ranks" | tee -a "${LOGDIR}/batch.log"
  (
    cd "${run_dir}"
    mpirun -np "${NPROC}" ./mitgcmuv
  ) >"${log}" 2>&1

  if ! grep -q "STOP NORMAL END" "${run_dir}/STDOUT.0000"; then
    echo "[$(date)] ${run} finished without STOP NORMAL END; see ${log}" | tee -a "${LOGDIR}/batch.log"
    exit 1
  fi
  echo "[$(date)] finished ${run}" | tee -a "${LOGDIR}/batch.log"
done

echo "[$(date)] all S70 runs finished" | tee -a "${LOGDIR}/batch.log"
