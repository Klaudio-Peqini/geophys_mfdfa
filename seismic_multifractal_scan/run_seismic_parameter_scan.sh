#!/usr/bin/env bash
#set -euo pipefail
set -uo pipefail

SCRIPT="analyze_seismic_multifractal.py"
CATALOG="eq_data_earthquake_reviewed_mag4.csv"
BASE_OUT="scan_outputs"

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

MASTER_CSV="${BASE_OUT}/all_parameter_scan_results.csv"
echo "run_id,bin,series,magmin,window_size_bins,window_step_bins,qmin,qmax,qstep,fit_smin,fit_smax,h2,h0,spectrum_width_alpha,alpha_at_fmax,falpha_max,delta_h_mean,delta_h_std,delta_h_min,delta_h_max,n_windows,summary_json,windowed_hq_csv,delta_h_csv,delta_h_png,delta_alpha_csv,delta_alpha_png" > "${MASTER_CSV}"

run_counter=0

for bin in "${BINS[@]}"; do
  for series in "${SERIES[@]}"; do
    for magmin in "${MAGMINS[@]}"; do
      for win in "${WINDOW_SIZES[@]}"; do
        for step in "${WINDOW_STEPS[@]}"; do
          if [[ "${win}" -eq 1095 && "${step}" -gt 182 ]]; then
            continue
          fi

          fit_smin=$(fit_smin_for_window "${win}")
          fit_smax=$(fit_smax_for_window "${win}")

          run_id=$(printf "run_%04d" "${run_counter}")
          outdir="${BASE_OUT}/${run_id}_${series}_${bin}_w${win}_s${step}_m${magmin}"
          mkdir -p "${outdir}"

          echo "============================================================"
          echo "Running ${run_id}"
          echo "bin=${bin}, series=${series}, magmin=${magmin}, window=${win}, step=${step}"
          echo "============================================================"

          cmd=(
            python3 "${SCRIPT}"
            --catalog "${CATALOG}"
            --out "${outdir}"
            --bin "${bin}"
            --series "${series}"
            --log1p
            --qmin "${QMIN}" --qmax "${QMAX}" --qstep "${QSTEP}"
            --fit_smin "${fit_smin}"
            --fit_smax "${fit_smax}"
            --window_size_bins "${win}"
            --window_step_bins "${step}"
            --window_q_values 0 1 2 3 4 5
            --n_shuffles 20
          )

          if [[ "${magmin}" != "none" ]]; then
            cmd+=(--magmin "${magmin}")
          fi

          "${cmd[@]}"

          python aggregate_seismic_results.py             --run_dir "${outdir}"             --master_csv "${MASTER_CSV}"             --run_id "${run_id}"             --bin "${bin}"             --series "${series}"             --magmin "${magmin}"             --window_size_bins "${win}"             --window_step_bins "${step}"             --qmin "${QMIN}"             --qmax "${QMAX}"             --qstep "${QSTEP}"             --fit_smin "${fit_smin}"             --fit_smax "${fit_smax}"

          run_counter=$((run_counter + 1))
        done
      done
    done
  done
done

echo
echo "Finished."
echo "Master CSV: ${MASTER_CSV}"
