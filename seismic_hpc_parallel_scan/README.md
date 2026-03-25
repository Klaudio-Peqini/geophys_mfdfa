# seismic_hpc_parallel_scan

A mature HPC/scientific pipeline for multifractal analysis of seismic catalogs with:

- cached binned-series reuse
- duplicate-free windowed metrics
- adaptive hybrid parallelism driven by `window_step_bins`
- full timing instrumentation
- resumable scans with skipped completed runs
- failed-job retry support
- environment capture for reproducibility
- benchmark and load-analysis tooling

## Folder contents

```text
seismic_hpc_parallel_scan/
├── analyze_seismic_multifractal_fast.py
├── prepare_binned_series_cache.py
├── aggregate_seismic_results.py
├── generate_jobs.sh
├── run_parallel_scan.sh
├── merge_summaries.sh
├── sample_cpu_load.sh
├── sweep_job_counts.sh
├── analyze_benchmark_and_loads.py
├── requirements.txt
└── README.md
```

## What changed in this mature version

### 1. No duplicated scientific work in windowed analysis
All windowed quantities are computed in a single pass and then reused for:

- $h(q)$ per window
- $Δh(t)$ via $h(0,t) - h(5,t)$
- $Δα(t)$
- `windowed_hq.csv`
- plots
- run summaries

There is no second pass that recomputes the same windows just to obtain derived quantities.

### 2. Adaptive hybrid parallel policy
The pipeline distinguishes between two execution modes:

- **job-level mode** for larger `window_step_bins`
- **hybrid mode** for denser windows (`window_step_bins <= ADAPTIVE_WINDOW_THRESHOLD`)

By default:

- large steps stay in the GNU Parallel outer-job queue
- dense steps are assigned to a hybrid queue and also receive a nontrivial `--window-workers` value inside Python

This prevents wasting inner-worker overhead on coarse jobs while accelerating dense-window runs.

### 3. Full timing instrumentation
Each run writes detailed timing information into `summary.json`, including:

- cache or catalog loading
- catalog filtering
- series building
- series export
- series plotting
- preprocessing
- ACF computation / writing / plotting
- full-series MFDFA computation / writing / plotting
- spectrum computation / writing / plotting
- windowed computation / writing / plotting
- summary writing
- total runtime

### 4. Resume and failure recovery
Completed runs are skipped automatically when these files exist:

- `summary.json`
- `run_summary.csv`
- `run_complete.ok`

The launcher also writes job logs and retries failed jobs once with GNU Parallel `--resume-failed`.

### 5. Reproducibility capture
Each run stores:

- `environment.json`
- `summary.json["environment"]`

The launcher also writes:

- `scan_outputs/pipeline_environment.json`

These capture machine and software details such as Python version, package versions, hostname, and thread-related environment variables.

### 6. Production mode includes step = 1 day
The default job generator now includes:

```bash
WINDOW_STEPS=(365 182 30 7 1)
```

so dense one-bin sliding-window scans are part of the normal production workflow.

## Quick start

### 1. Create the environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

Install GNU Parallel if needed:

```bash
sudo apt install parallel
```

### 2. Put the catalog in the folder

Expected default name:

```text
eq_data_earthquake_reviewed_mag4.csv
```

### 3. Generate jobs

```bash
bash generate_jobs.sh
```

### 4. Run the scan

```bash
bash run_parallel_scan.sh
```

### 5. Merge summaries

```bash
bash merge_summaries.sh
```

The merged table is written to:

```text
scan_outputs/all_parameter_scan_results.csv
```

## Cluster usage

A typical cluster launch is:

```bash
THREADS=48 JOBS_STANDARD=20 JOBS_HYBRID=6 bash run_parallel_scan.sh
```

Important environment defaults inside the launcher:

```bash
OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
```

This prevents hidden oversubscription.

### Recommended strategy
- For coarse scans, use a larger `JOBS_STANDARD`
- For dense scans, keep `JOBS_HYBRID` smaller because each hybrid job can also use inner window workers

## Benchmark instructions

To benchmark different outer-job settings:

```bash
JOB_COUNTS="4 6 8 10 12 14" \
HYBRID_COUNTS="1 1 2 2 3 3" \
bash sweep_job_counts.sh
```

Then analyze the results:

```bash
python3 analyze_benchmark_and_loads.py --bench-dir benchmarks
```

This produces:
- scaling summaries
- speedup and efficiency plots
- CPU-load PDF plots

## Important configurable variables

### In `generate_jobs.sh`

- `WINDOW_STEPS`
- `WINDOW_SIZES`
- `ADAPTIVE_WINDOW_THRESHOLD`
- `WINDOW_WORKER_CAP`
- `SKIP_PLOTS`
- `LIGHT_OUTPUT`
- `SKIP_COMPLETED`

### In `run_parallel_scan.sh`

- `THREADS`
- `JOBS_STANDARD`
- `JOBS_HYBRID`
- `RETRY_FAILED`
- `REBUILD_JOBS`
- `RUN_MERGE`

## Outputs

Each run directory contains the main scientific products:

- `summary.json`
- `environment.json`
- `status.json`
- `multifractal_spectrum.csv`
- `windowed_hq.csv`
- `delta_h_h0_h5.csv`
- `delta_alpha_windowed.csv`
- `run_summary.csv`
- `run_complete.ok`

Global files:

- `scan_outputs/all_parameter_scan_results.csv`
- `scan_outputs/pipeline_environment.json`

Launcher logs:

- `joblogs/joblevel_joblog.tsv`
- `joblogs/hybrid_joblog.tsv`
- `joblogs/joblevel_failed_commands.txt`
- `joblogs/hybrid_failed_commands.txt`

## Cache files

The cache builder writes both:

- CSV cache files
- compressed NPZ cache files

Examples:

```text
cache_series/series_energy_1D_mnone.csv
cache_series/series_energy_1D_mnone.npz
```

The NPZ cache is preferred because it avoids repeated CSV parsing overhead.

## Troubleshooting

### `EmptyDataError: No columns to parse from file`
The catalog exists but is empty or invalid. Check:

```bash
ls -lh eq_data_earthquake_reviewed_mag4.csv
head -5 eq_data_earthquake_reviewed_mag4.csv
```

### `TypeError: len() of unsized object` when reading NPZ
Use the provided mature scripts. The NPZ loader already handles scalar and one-element metadata formats safely.

### GNU Parallel is missing
Install it:

```bash
sudo apt install parallel
```

### The scan reruns too much after a crash
Set:

```bash
SKIP_COMPLETED=1
```

and rerun. Completed jobs will be skipped automatically.

### Hybrid jobs overload the machine
Lower:

```bash
JOBS_HYBRID
```

or reduce:

```bash
WINDOW_WORKER_CAP
```

## Recommended production pattern

1. benchmark machine-specific `JOBS_STANDARD` and `JOBS_HYBRID`
2. keep `SKIP_PLOTS=1` for broad sweeps
3. use cached NPZ inputs
4. enable `SKIP_COMPLETED=1`
5. inspect merged summaries and rerun only specific subsets if needed
