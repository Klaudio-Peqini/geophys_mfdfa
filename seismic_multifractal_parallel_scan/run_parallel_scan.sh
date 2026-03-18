#!/usr/bin/env bash
#set -euo pipefail
set -uo pipefall

# For a 16-thread machine, a good starting point is 8 parallel jobs.
JOBS=12

if ! command -v parallel >/dev/null 2>&1; then
  echo "GNU parallel is not installed."
  echo "Install it with: sudo apt install parallel"
  exit 1
fi

source generate_jobs.sh
parallel -j "${JOBS}" --joblog parallel_joblog.txt < jobs.txt

echo "Parallel scan finished."
