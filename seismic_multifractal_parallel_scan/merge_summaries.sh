#!/usr/bin/env bash
set -euo pipefail

OUT="scan_outputs/all_parameter_scan_results.csv"

find scan_outputs -name run_summary.csv | sort | xargs awk 'FNR==1 && NR!=1 {next} {print}' > "${OUT}"

echo "Merged summary CSV written to ${OUT}"
