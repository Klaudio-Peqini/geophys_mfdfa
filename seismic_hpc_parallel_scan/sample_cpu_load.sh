#!/usr/bin/env bash
set -euo pipefail

OUT="${1:-cpu_load_samples.csv}"
INTERVAL="${INTERVAL:-1}"
DURATION="${DURATION:-0}"

cpu_count=$(grep -c '^cpu[0-9]' /proc/stat)
printf 'timestamp_epoch,cpu_id,usage_percent\n' > "${OUT}"

read_stat() {
  awk '/^cpu[0-9]+ / {print $1,$2+$3+$4+$5+$6+$7+$8,$5+$6}' /proc/stat
}

declare -A prev_total prev_idle
while read -r cpu total idle; do
  prev_total["$cpu"]="$total"
  prev_idle["$cpu"]="$idle"
done < <(read_stat)

start_ts=$(date +%s)
while true; do
  sleep "${INTERVAL}"
  ts=$(date +%s)
  while read -r cpu total idle; do
    old_total=${prev_total["$cpu"]}
    old_idle=${prev_idle["$cpu"]}
    dt=$((total - old_total))
    di=$((idle - old_idle))
    usage=0
    if [[ ${dt} -gt 0 ]]; then
      usage=$(awk -v dt="$dt" -v di="$di" 'BEGIN {printf "%.4f", 100.0*(dt-di)/dt}')
    fi
    cpu_id=${cpu#cpu}
    printf '%s,%s,%s\n' "$ts" "$cpu_id" "$usage" >> "${OUT}"
    prev_total["$cpu"]="$total"
    prev_idle["$cpu"]="$idle"
  done < <(read_stat)

  if [[ "${DURATION}" != "0" ]]; then
    now=$(date +%s)
    if (( now - start_ts >= DURATION )); then
      break
    fi
  fi
done
