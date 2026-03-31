#!/usr/bin/env bash
set -euo pipefail

# Batch plot generator for all run_* folders inside scan_outputs
# Example usage:
#   bash run_plots_all.sh
#   bash run_plots_all.sh /path/to/scan_outputs
#   bash run_plots_all.sh /path/to/scan_outputs /path/to/plot_weekly_meeting_figures.py

SCAN_OUTPUTS_DIR="${1:-scan_outputs}"
PLOT_SCRIPT="${2:-$(dirname "$0")/plot_figures.py}"
LOG_FILE="${3:-plots_batch.log}"

if [[ ! -d "$SCAN_OUTPUTS_DIR" ]]; then
    echo "Error: scan_outputs directory not found: $SCAN_OUTPUTS_DIR" >&2
    exit 1
fi

if [[ ! -f "$PLOT_SCRIPT" ]]; then
    echo "Error: plotting script not found: $PLOT_SCRIPT" >&2
    exit 1
fi

if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
elif command -v python >/dev/null 2>&1; then
    PYTHON_BIN="python"
else
    echo "Error: neither python3 nor python was found in PATH." >&2
    exit 1
fi

TOTAL=0
OK=0
FAIL=0

: > "$LOG_FILE"

echo "Batch plotting started at $(date)" | tee -a "$LOG_FILE"
echo "scan_outputs directory: $SCAN_OUTPUTS_DIR" | tee -a "$LOG_FILE"
echo "plotting script: $PLOT_SCRIPT" | tee -a "$LOG_FILE"
echo "python executable: $PYTHON_BIN" | tee -a "$LOG_FILE"
echo | tee -a "$LOG_FILE"

shopt -s nullglob
RUN_DIRS=("$SCAN_OUTPUTS_DIR"/run_*)
shopt -u nullglob

if [[ ${#RUN_DIRS[@]} -eq 0 ]]; then
    echo "No run_* folders found inside $SCAN_OUTPUTS_DIR" | tee -a "$LOG_FILE"
    exit 1
fi

for RUN_DIR in "${RUN_DIRS[@]}"; do
    [[ -d "$RUN_DIR" ]] || continue
    TOTAL=$((TOTAL + 1))
    BASENAME="$(basename "$RUN_DIR")"
    OUT_DIR="$RUN_DIR/plots_weekly_meeting"

    echo "[$TOTAL] Processing $BASENAME ..." | tee -a "$LOG_FILE"

    if "$PYTHON_BIN" "$PLOT_SCRIPT" --input-dir "$RUN_DIR" --output-dir "$OUT_DIR" >> "$LOG_FILE" 2>&1; then
        echo "[$TOTAL] OK   -> $OUT_DIR" | tee -a "$LOG_FILE"
        OK=$((OK + 1))
    else
        echo "[$TOTAL] FAIL -> $BASENAME" | tee -a "$LOG_FILE"
        FAIL=$((FAIL + 1))
    fi

    echo | tee -a "$LOG_FILE"
done

echo "Batch plotting finished at $(date)" | tee -a "$LOG_FILE"
echo "Total folders : $TOTAL" | tee -a "$LOG_FILE"
echo "Successful    : $OK" | tee -a "$LOG_FILE"
echo "Failed        : $FAIL" | tee -a "$LOG_FILE"

if [[ $FAIL -gt 0 ]]; then
    exit 2
fi
