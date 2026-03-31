#!/usr/bin/env bash
set -euo pipefail

# Package only delta_h_h0_h5.csv for selected run folders.
#
# Usage:
#   bash package_delta_h_only.sh
#   bash package_delta_h_only.sh scan_outputs
#   bash package_delta_h_only.sh scan_outputs delta_h_bundle
#   bash package_delta_h_only.sh scan_outputs delta_h_bundle "run_counts_1D_m5.0_*"
#
# Result:
#   Creates <bundle_name>_<timestamp>.zip containing:
#     <bundle>/<run_folder>/delta_h_h0_h5.csv
#
# Notes:
#   - By default, only step=1 and step=7 folders are included.
#   - You can override the run-name pattern with a custom glob as 3rd argument.

SCAN_DIR="${1:-scan_outputs}"
BUNDLE_NAME="${2:-delta_h_bundle}"
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

  if [[ -f "$run_dir/delta_h_h0_h5.csv" ]]; then
    cp "$run_dir/delta_h_h0_h5.csv" "$dest_dir/"
    echo "  + delta_h_h0_h5.csv" >> "$MANIFEST"
    ((included+=1))
  else
    echo "  - delta_h_h0_h5.csv (missing)" >> "$MANIFEST"
    rmdir "$dest_dir" 2>/dev/null || true
    ((skipped+=1))
  fi

  echo >> "$MANIFEST"
done
shopt -u nullglob

if [[ "$seen_any" -eq 0 ]]; then
  echo "ERROR: no run folders matched '$RUN_GLOB' under '$SCAN_DIR'."
  rm -rf "$STAGING_DIR"
  exit 1
fi

if [[ "$included" -eq 0 ]]; then
  echo "ERROR: matching folders were found, but none contained delta_h_h0_h5.csv."
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
