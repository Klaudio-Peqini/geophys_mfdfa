#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
THREADS=${THREADS:-$(nproc)}
JOB_COUNTS=(${JOB_COUNTS:-2 4 6 8 10 12 14 16})
HYBRID_COUNTS=(${HYBRID_COUNTS:-1 1 2 2 3 3 4 4})
BENCH_DIR="${BENCH_DIR:-${HERE}/benchmarks}"
LOAD_INTERVAL="${LOAD_INTERVAL:-1}"
REBUILD_JOBS_ONCE=${REBUILD_JOBS_ONCE:-1}

mkdir -p "${BENCH_DIR}"

if [[ "${REBUILD_JOBS_ONCE}" == "1" ]]; then
  bash "${HERE}/generate_jobs.sh"
fi

idx=0
for jobs in "${JOB_COUNTS[@]}"; do
  hybrid="${HYBRID_COUNTS[$idx]}"
  idx=$((idx + 1))
  run_dir="${BENCH_DIR}/jobs_${jobs}"
  mkdir -p "${run_dir}"

  rm -rf "${HERE}/scan_outputs"
  mkdir -p "${HERE}/scan_outputs"

  load_csv="${run_dir}/cpu_load_samples.csv"
  INTERVAL="${LOAD_INTERVAL}" DURATION=0 bash "${HERE}/sample_cpu_load.sh" "${load_csv}" &
  sampler_pid=$!

  start=$(date +%s)
  set +e
  /usr/bin/time -f '%e' -o "${run_dir}/elapsed_seconds.txt" \
    env THREADS="${THREADS}" JOBS_STANDARD="${jobs}" JOBS_HYBRID="${hybrid}" REBUILD_JOBS=0 JOBLOG_DIR="${run_dir}/joblogs" \
    bash "${HERE}/run_parallel_scan.sh" > "${run_dir}/stdout.log" 2> "${run_dir}/stderr.log"
  rc=$?
  set -e
  end=$(date +%s)

  kill "${sampler_pid}" >/dev/null 2>&1 || true
  wait "${sampler_pid}" 2>/dev/null || true

  printf 'jobs_standard,jobs_hybrid,wall_clock_sec,return_code\n%s,%s,%s,%s\n' "${jobs}" "${hybrid}" "$((end - start))" "${rc}" > "${run_dir}/run_meta.csv"
  if [[ "${rc}" -ne 0 ]]; then
    echo "Benchmark failed for JOBS_STANDARD=${jobs}, JOBS_HYBRID=${hybrid}. See ${run_dir}/stderr.log"
    exit "${rc}"
  fi

  cp "${HERE}/scan_outputs/all_parameter_scan_results.csv" "${run_dir}/all_parameter_scan_results.csv"
  echo "Completed benchmark for JOBS_STANDARD=${jobs}, JOBS_HYBRID=${hybrid}"
done
