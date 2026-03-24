#!/usr/bin/env bash
set -euo pipefail

SCRIPT="./analyze_seismic_multifractal_fast.py"
AGG="./aggregate_seismic_results.py"
CATALOG="./eq_data_earthquake_reviewed_mag4.csv"
CACHE_DIR="./cache_series"
BASE_OUT="./scan_outputs"

BINS=("1D" "3D" "7D")
SERIES=("energy" "counts")
MAGMINS=("none" "5.0")

QMIN=-10
QMAX=10
QSTEP=1
WINDOW_SIZES=(3650 1825 1670 1095)
WINDOW_STEPS=(365 182 30 7)

fit_smin_for_window() {
  case "$1" in
    3650) echo 121 ;;
    1825) echo 60 ;;
    1670) echo 60 ;;
    1095) echo 30 ;;
    *) echo 60 ;;
  esac
}

fit_smax_for_window() {
  case "$1" in
    3650) echo 1460 ;;
    1825) echo 730 ;;
    1670) echo 730 ;;
    1095) echo 365 ;;
    *) echo 365 ;;
  esac
}

mkdir -p "${BASE_OUT}"
rm -f jobs.txt
: > jobs.txt

python3 ./prepare_binned_series_cache.py \
  --catalog "${CATALOG}" \
  --cache-dir "${CACHE_DIR}" \
  --bins "${BINS[@]}" \
  --series "${SERIES[@]}" \
  --magmins "${MAGMINS[@]}"

run_counter=0
for bin in "${BINS[@]}"; do
  for series in "${SERIES[@]}"; do
    for magmin in "${MAGMINS[@]}"; do
      cache_csv="${CACHE_DIR}/series_${series}_${bin}_m${magmin}.csv"
      for win in "${WINDOW_SIZES[@]}"; do
        for step in "${WINDOW_STEPS[@]}"; do
          if [[ "${win}" -eq 1095 && "${step}" -gt 182 ]]; then
            continue
          fi

          fit_smin=$(fit_smin_for_window "${win}")
          fit_smax=$(fit_smax_for_window "${win}")
          run_id=$(printf "run_%04d" "${run_counter}")
          outdir="${BASE_OUT}/${run_id}_${series}_${bin}_w${win}_s${step}_m${magmin}"

          cmd="python3 ${SCRIPT} \
            --input-series-csv ${cache_csv} \
            --out ${outdir} \
            --bin ${bin} \
            --series ${series} \
            --log1p \
            --skip-series-export \
            --qmin ${QMIN} --qmax ${QMAX} --qstep ${QSTEP} \
            --fit_smin ${fit_smin} \
            --fit_smax ${fit_smax} \
            --window_size_bins ${win} \
            --window_step_bins ${step} \
            --window_q_values 0 1 2 3 4 5 \
            --n_shuffles 20"

          cmd+=" && python3 ${AGG} \
            --run_dir ${outdir} \
            --run_id ${run_id} \
            --bin ${bin} \
            --series ${series} \
            --magmin ${magmin} \
            --window_size_bins ${win} \
            --window_step_bins ${step} \
            --qmin ${QMIN} \
            --qmax ${QMAX} \
            --qstep ${QSTEP} \
            --fit_smin ${fit_smin} \
            --fit_smax ${fit_smax}"

          echo "${cmd}" >> jobs.txt
          run_counter=$((run_counter + 1))
        done
      done
    done
  done
done

echo "Created jobs.txt with ${run_counter} jobs."
