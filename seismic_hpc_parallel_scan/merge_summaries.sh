#!/usr/bin/env bash
set -euo pipefail

OUT="scan_outputs/all_parameter_scan_results.csv"
mkdir -p "$(dirname "${OUT}")"

mapfile -t files < <(find scan_outputs -name run_summary.csv -print | sort)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "No run_summary.csv files were found under scan_outputs/."
  exit 1
fi

awk 'FNR==1 && NR!=1 {next} {print}' "${files[@]}" > "${OUT}"
echo "Merged summary CSV written to ${OUT}"
