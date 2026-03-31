#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"

OUTDIR="${HERE}/scan_outputs"
CACHEDIR="${HERE}/cache_series"
ANALYZER="${HERE}/analyze_seismic_multifractal_fast.py"
AGGREGATOR="${HERE}/aggregate_seismic_results.py"

mkdir -p "${OUTDIR}"

WINDOW_SIZES=(3650 1825 1670 1095)
WINDOW_STEPS=(365 182 30 7 1)

BINS=("1D" "3D" "7D")
SERIES=("energy" "counts")
MAGMINS=("none" "5.0")

QMIN=-10
QMAX=10
QSTEP=1

rm -f "${HERE}/jobs.txt"

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

for bin in "${BINS[@]}"; do
  for series in "${SERIES[@]}"; do
    for mag in "${MAGMINS[@]}"; do

      csv_file="${CACHEDIR}/series_${series}_${bin}_m${mag}.csv"
      npz_file="${CACHEDIR}/series_${series}_${bin}_m${mag}.npz"

      if [[ -f "${csv_file}" ]]; then
        INPUT_ARG="--input-series-csv ${csv_file}"
      elif [[ -f "${npz_file}" ]]; then
        INPUT_ARG="--input-series-npz ${npz_file}"
      else
        echo "ERROR: missing cache for ${series} ${bin} m=${mag}"
        echo "Expected one of:"
        echo "  ${csv_file}"
        echo "  ${npz_file}"
        exit 1
      fi

      for win in "${WINDOW_SIZES[@]}"; do
        for step in "${WINDOW_STEPS[@]}"; do

          smin=$(fit_smin_for_window "$win")
          smax=$(fit_smax_for_window "$win")

          run_id="run_${series}_${bin}_m${mag}_w${win}_s${step}"
          out="${OUTDIR}/${run_id}"

          cmd="python3 ${ANALYZER} \
            ${INPUT_ARG} \
            --out ${out} \
            --bin ${bin} \
            --series ${series} \
            --window-size-bins ${win} \
            --window-step-bins ${step} \
            --fit-smin ${smin} \
            --fit-smax ${smax} \
            --window-workers 1 \
            --window-q-values 0 1 2 3 4 5 \
            --n-shuffles 20 \
            --skip-plots \
            && python3 ${AGGREGATOR} \
            --run_dir ${out} \
            --run_id ${run_id} \
            --bin ${bin} \
            --series ${series} \
            --magmin ${mag} \
            --window_size_bins ${win} \
            --window_step_bins ${step} \
            --qmin ${QMIN} \
            --qmax ${QMAX} \
            --qstep ${QSTEP} \
            --fit_smin ${smin} \
            --fit_smax ${smax}"

          echo "$cmd" >> "${HERE}/jobs.txt"

        done
      done

    done
  done
done

echo "Created jobs.txt with $(wc -l < "${HERE}/jobs.txt") runnable jobs."
