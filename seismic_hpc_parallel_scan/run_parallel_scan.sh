#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"

if ! command -v parallel >/dev/null 2>&1; then
  echo "GNU parallel is not installed."
  echo "Install it with: sudo apt install parallel"
  exit 1
fi

detect_threads() {
  if [[ -n "${THREADS:-}" ]]; then
    echo "${THREADS}"
    return
  fi

  local t=""
  t=$(python3 - <<'PY'
import os
n = os.cpu_count()
print(n if n is not None else 1)
PY
)
  if [[ "$t" =~ ^[0-9]+$ ]] && (( t >= 1 )); then
    echo "$t"
    return
  fi

  t=$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)
  if [[ "$t" =~ ^[0-9]+$ ]] && (( t >= 1 )); then
    echo "$t"
    return
  fi

  t=$(nproc 2>/dev/null || true)
  if [[ "$t" =~ ^[0-9]+$ ]] && (( t >= 1 )); then
    echo "$t"
    return
  fi

  echo 1
}

THREADS="$(detect_threads)"
JOBS_STANDARD="${JOBS_STANDARD:-${JOBS:-$(( THREADS > 2 ? (3 * THREADS) / 4 : 1 ))}}"
REBUILD_JOBS="${REBUILD_JOBS:-1}"
RETRY_FAILED="${RETRY_FAILED:-1}"
RUN_MERGE="${RUN_MERGE:-1}"
JOBLOG_DIR="${JOBLOG_DIR:-${HERE}/joblogs}"
mkdir -p "${JOBLOG_DIR}"
mkdir -p "${HERE}/scan_outputs"

if [[ "${REBUILD_JOBS}" == "1" ]]; then
  bash "${HERE}/generate_jobs.sh"
fi

if [[ ! -f "${HERE}/jobs.txt" ]]; then
  echo "jobs.txt was not found. Generate jobs first."
  exit 1
fi

rm -f "${HERE}/jobs_joblevel.txt"
rm -f "${HERE}/jobs_hybrid.txt"
rm -f "${HERE}"/jobs_hybrid_step_*.txt
: > "${HERE}/jobs_joblevel.txt"
: > "${HERE}/jobs_hybrid.txt"

extract_step() {
  python3 - "$1" <<'PY'
import re, sys
cmd = sys.argv[1]
m = re.search(r'--window_step_bins\s+([0-9]+)', cmd)
print(m.group(1) if m else "")
PY
}

while IFS= read -r cmd; do
  [[ -z "$cmd" ]] && continue
  step="$(extract_step "$cmd")"
  if [[ "$step" =~ ^[0-9]+$ ]] && (( step <= 30 )); then
    echo "$cmd" >> "${HERE}/jobs_hybrid.txt"
    echo "$cmd" >> "${HERE}/jobs_hybrid_step_${step}.txt"
  else
    echo "$cmd" >> "${HERE}/jobs_joblevel.txt"
  fi
done < "${HERE}/jobs.txt"

count_lines() {
  local f="$1"
  [[ -s "$f" ]] && wc -l < "$f" || echo 0
}

JOBLEVEL_COUNT="$(count_lines "${HERE}/jobs_joblevel.txt")"
HYBRID_COUNT="$(count_lines "${HERE}/jobs_hybrid.txt")"

echo "Detected THREADS=${THREADS}"
echo "Using JOBS_STANDARD=${JOBS_STANDARD}"
echo "  job-level queue : ${JOBLEVEL_COUNT}"
echo "  hybrid queue    : ${HYBRID_COUNT}"

python3 - <<'PY' "${HERE}" "${THREADS}" "${JOBS_STANDARD}"
import json, os, platform, socket, sys
from importlib import metadata as md
from pathlib import Path
from datetime import datetime, timezone

here, threads, jobs_standard = sys.argv[1:]

def ver(name):
    try:
        return md.version(name)
    except Exception:
        return "unknown"

payload = {
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "hostname": socket.gethostname(),
    "platform": platform.platform(),
    "python_version": sys.version.replace("\n", " "),
    "cpu_count_logical": os.cpu_count(),
    "threads_visible": int(threads),
    "jobs_standard": int(jobs_standard),
    "packages": {k: ver(k) for k in ["numpy", "pandas", "matplotlib", "numba"]},
    "thread_env": {k: os.environ.get(k) for k in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]},
}
with open(Path(here, "scan_outputs", "pipeline_environment.json"), "w", encoding="utf-8") as f:
    json.dump(payload, f, indent=2)
PY

calc_hybrid_policy() {
  local step="$1"
  local threads="$2"
  local outer=1
  local inner=1

  if (( step <= 1 )); then
    if (( threads <= 8 )); then
      outer=1; inner=$(( threads > 1 ? threads - 1 : 1 ))
    elif (( threads <= 16 )); then
      outer=2; inner=6
    elif (( threads <= 32 )); then
      outer=4; inner=6
    else
      outer=$(( threads / 8 )); (( outer < 4 )) && outer=4
      inner=6
    fi
  elif (( step <= 7 )); then
    if (( threads <= 8 )); then
      outer=1; inner=$(( threads > 1 ? threads - 1 : 1 ))
    elif (( threads <= 16 )); then
      outer=3; inner=4
    elif (( threads <= 32 )); then
      outer=4; inner=4
    else
      outer=$(( threads / 8 )); (( outer < 4 )) && outer=4
      inner=4
    fi
  else
    if (( threads <= 8 )); then
      outer=2; inner=2
    elif (( threads <= 16 )); then
      outer=6; inner=2
    elif (( threads <= 32 )); then
      outer=8; inner=2
    else
      outer=$(( (3 * threads) / 8 )); (( outer < 8 )) && outer=8
      inner=2
    fi
  fi

  (( outer < 1 )) && outer=1
  (( inner < 1 )) && inner=1
  echo "${outer} ${inner}"
}

run_queue() {
  local label="$1"
  local jobfile="$2"
  local jobs="$3"
  local window_workers="$4"
  local joblog
  local rc=0
  local failed=0

  joblog="${JOBLOG_DIR}/${label}_joblog.tsv"

  if [[ ! -s "${jobfile}" ]]; then
    echo "Queue ${label}: no jobs to run."
    return 0
  fi

  echo "Running ${label} queue with -j ${jobs} (WINDOW_WORKERS=${window_workers})"
  set +e
  WINDOW_WORKERS="${window_workers}" MFDFA_WINDOW_WORKERS="${window_workers}" \
    parallel --bar --joblog "${joblog}" -j "${jobs}" < "${jobfile}"
  rc=$?
  set -e

  if [[ "${RETRY_FAILED}" == "1" && -s "${joblog}" ]]; then
    failed=$(awk 'NR>1 && $7 != 0 {c++} END {print c+0}' "${joblog}")
    if (( failed > 0 )); then
      echo "Retrying ${failed} failed jobs in ${label} queue via --resume-failed"
      set +e
      WINDOW_WORKERS="${window_workers}" MFDFA_WINDOW_WORKERS="${window_workers}" \
        parallel --resume-failed --joblog "${joblog}" -j "${jobs}" < "${jobfile}"
      set -e
    fi
  fi

  if [[ -s "${joblog}" ]]; then
    awk 'NR>1 && $7 != 0 {print $9}' "${joblog}" > "${JOBLOG_DIR}/${label}_failed_commands.txt" || true
  fi

  return ${rc}
}

run_queue "joblevel" "${HERE}/jobs_joblevel.txt" "${JOBS_STANDARD}" 1

hybrid_rc=0
shopt -s nullglob
for stepfile in "${HERE}"/jobs_hybrid_step_*.txt; do
  [[ ! -s "$stepfile" ]] && continue
  step="${stepfile##*_step_}"
  step="${step%.txt}"
  if [[ ! "$step" =~ ^[0-9]+$ ]]; then
    continue
  fi
  read -r outer_jobs inner_workers <<< "$(calc_hybrid_policy "$step" "$THREADS")"
  run_queue "hybrid_step_${step}" "$stepfile" "$outer_jobs" "$inner_workers" || hybrid_rc=$?
done
shopt -u nullglob

if [[ "${RUN_MERGE}" == "1" ]]; then
  bash "${HERE}/merge_summaries.sh"
fi

echo "Parallel scan finished."
echo "Logs saved under ${JOBLOG_DIR}"
exit ${hybrid_rc}
