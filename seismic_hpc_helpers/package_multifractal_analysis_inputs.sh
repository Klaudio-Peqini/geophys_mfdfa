#!/usr/bin/env bash
set -euo pipefail

# Package the most relevant outputs for detailed multifractal interpretation.
# Default focus: runs with small window shifts step=1 and step=7.
#
# Usage examples:
#   bash package_multifractal_analysis_inputs_fixed.sh
#   bash package_multifractal_analysis_inputs_fixed.sh scan_outputs
#   bash package_multifractal_analysis_inputs_fixed.sh scan_outputs multifractal_analysis_bundle
#   bash package_multifractal_analysis_inputs_fixed.sh scan_outputs multifractal_analysis_bundle "run_counts_1D_m5.0_*"
#
# Result:
#   Creates <bundle_name>.zip containing, for each selected run folder:
#     - summary.json
#     - run_summary.csv
#     - acf.csv
#     - binned_series.csv
#     - delta_alpha_windowed.csv
#     - delta_h_h0_h5.csv
#     - multifractal_spectrum.csv
#     - shuffle_test_summary.csv
#     - windowed_hq.csv
#     - mfdfa_Fq.csv (if present)
#     - plots_weekly_meeting/* (if present)
#
# Notes:
#   - By default, only step=1 and step=7 folders are included.
#   - You can override the run-name pattern with a custom glob as 3rd argument.

SCAN_DIR="${1:-scan_outputs}"
BUNDLE_NAME="${2:-multifractal_analysis_bundle}"
RUN_GLOB="${3:-run_*}"
CUSTOM_GLOB_PROVIDED=0
if [[ $# -ge 3 ]]; then
  CUSTOM_GLOB_PROVIDED=1
fi

if [[ ! -d "$SCAN_DIR" ]]; then
  echo "ERROR: scan directory '$SCAN_DIR' does not exist."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
STAGING_DIR="${BUNDLE_NAME}_${STAMP}"
ZIP_FILE="${BUNDLE_NAME}_${STAMP}.zip"
MANIFEST="${STAGING_DIR}/MANIFEST.txt"

mkdir -p "$STAGING_DIR"

echo "Packaging from: $SCAN_DIR" | tee "$MANIFEST"
echo "Run glob: $RUN_GLOB" | tee -a "$MANIFEST"
echo "Created: $(date)" | tee -a "$MANIFEST"
echo >> "$MANIFEST"

included=0
skipped=0
seen_any=0

shopt -s nullglob
for run_dir in "$SCAN_DIR"/$RUN_GLOB; do
  [[ -d "$run_dir" ]] || continue
  seen_any=1

  run_base="$(basename "$run_dir")"

  # Default focus: only small-shift runs s1 or s7.
  if [[ "$CUSTOM_GLOB_PROVIDED" -eq 0 ]]; then
    if [[ ! "$run_base" =~ _s(1|7)$ ]]; then
      ((skipped+=1))
      continue
    fi
  fi

  dest_dir="$STAGING_DIR/$run_base"
  mkdir -p "$dest_dir"

  echo "[$run_base]" >> "$MANIFEST"

  for f in \
    summary.json \
    run_summary.csv \
    acf.csv \
    binned_series.csv \
    delta_alpha_windowed.csv \
    delta_h_h0_h5.csv \
    multifractal_spectrum.csv \
    shuffle_test_summary.csv \
    windowed_hq.csv \
    mfdfa_Fq.csv
  do
    if [[ -f "$run_dir/$f" ]]; then
      cp "$run_dir/$f" "$dest_dir/"
      echo "  + $f" >> "$MANIFEST"
    else
      echo "  - $f (missing)" >> "$MANIFEST"
    fi
  done

  if [[ -d "$run_dir/plots_weekly_meeting" ]]; then
    mkdir -p "$dest_dir/plots_weekly_meeting"
    cp -r "$run_dir/plots_weekly_meeting/." "$dest_dir/plots_weekly_meeting/"
    echo "  + plots_weekly_meeting/" >> "$MANIFEST"
  else
    echo "  - plots_weekly_meeting/ (missing)" >> "$MANIFEST"
  fi

  echo >> "$MANIFEST"
  ((included+=1))
done
shopt -u nullglob

if [[ "$seen_any" -eq 0 ]]; then
  echo "ERROR: no run folders matched '$RUN_GLOB' under '$SCAN_DIR'."
  rm -rf "$STAGING_DIR"
  exit 1
fi

if [[ "$included" -eq 0 ]]; then
  echo "ERROR: matching folders were found, but none satisfied the selection rules."
  echo "Hint: by default only _s1 and _s7 runs are included."
  rm -rf "$STAGING_DIR"
  exit 1
fi

zip -qr "$ZIP_FILE" "$STAGING_DIR"

echo "Included runs : $included"
echo "Skipped runs  : $skipped"
echo "Staging dir   : $STAGING_DIR"
echo "ZIP archive   : $ZIP_FILE"
echo
echo "Done. You can upload: $ZIP_FILE"
