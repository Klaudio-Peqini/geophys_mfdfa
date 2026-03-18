# seismic_multifractal_parallel_scan

Parallel parameter-space scan utilities for seismic multifractal analysis.

## Files

- `generate_jobs.sh`  
  Builds `jobs.txt`, one command per line.
- `run_parallel_scan.sh`  
  Runs all jobs with GNU Parallel.
- `merge_summaries.sh`  
  Merges all per-run summaries into one final CSV.
- `aggregate_seismic_results.py`  
  Creates:
  - `delta_h_h0_h5.csv`
  - `delta_h_h0_h5.png`
  - one `run_summary.csv` per run

## Expected companion files

Keep these scripts in the same folder as:
- `analyze_seismic_multifractal.py`
- `eq_data_earthquake_reviewed_mag4.csv`

## Recommended usage

```bash
chmod +x generate_jobs.sh
chmod +x run_parallel_scan.sh
chmod +x merge_summaries.sh

./run_parallel_scan.sh
./merge_summaries.sh
```

## Final result

```bash
scan_outputs/all_parameter_scan_results.csv
```

## Parallel jobs

For a 16-thread machine, start with:

```bash
JOBS=8
```

inside `run_parallel_scan.sh`.

Then, if the machine remains comfortable in RAM and I/O usage, try 10 or 12.
