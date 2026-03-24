#!/usr/bin/env bash
set -euo pipefail

# Keep each Python process single-threaded internally.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

if ! command -v parallel >/dev/null 2>&1; then
  echo "GNU parallel is not installed."
  echo "Install it with: sudo apt install parallel"
  exit 1
fi

THREADS=${THREADS:-$(nproc)}
JOBS=${JOBS:-$(( THREADS > 2 ? THREADS / 2 : 1 ))}

bash ./generate_jobs.sh
parallel --bar --joblog parallel_joblog.txt -j "${JOBS}" < jobs.txt

echo "Parallel scan finished with ${JOBS} concurrent jobs on ${THREADS} visible CPU threads."
