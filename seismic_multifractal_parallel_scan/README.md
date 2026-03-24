# seismic-parallel-refactor

A refactored workflow for running large multifractal seismic parameter scans efficiently on a laptop, workstation, or HPC node.

This folder is designed to be copied directly into a GitHub repository. Its purpose is not only to run the scan, but also to make the **parallelization logic explicit and reproducible**.

---

## 1. Why this refactor was needed

In the original workflow, each parallel job was doing much more than the scientifically essential work:

- re-reading the same earthquake catalog CSV,
- re-parsing timestamps,
- re-filtering and re-binning the same series,
- recomputing some windowed quantities twice,
- producing many PNG files per run,
- and then re-reading intermediate CSV files again during aggregation.

When this is repeated across a large parameter grid, the wall time starts to be dominated not only by the MFDFA itself, but also by repeated **I/O**, repeated **preprocessing**, and repeated **postprocessing**.

The refactor therefore applies three main principles:

1. **Cache once, reuse many times**.
2. **Keep one analysis job single-threaded internally**.
3. **Parallelize at the job level, not inside every numerical kernel**.

---

## 2. Folder contents

```text
seismic_parallel_refactor/
├── analyze_seismic_multifractal_fast.py
├── prepare_binned_series_cache.py
├── aggregate_seismic_results.py
├── generate_jobs.sh
├── run_parallel_scan.sh
├── merge_summaries.sh
├── requirements.txt
├── .gitignore
└── README.md
```

### Main role of each file

- `analyze_seismic_multifractal_fast.py`  
  Main optimized analysis script. It can work either from the raw earthquake catalog or from a precomputed cached series.

- `prepare_binned_series_cache.py`  
  Builds the reusable cached time series for each `(bin, series, magmin)` combination.

- `aggregate_seismic_results.py`  
  Produces a compact one-row summary for each run without doing any redundant plotting.

- `generate_jobs.sh`  
  Generates `jobs.txt` after first building the reusable cache.

- `run_parallel_scan.sh`  
  Executes the jobs with GNU parallel, while forcing NumPy/BLAS backends to stay single-threaded per job.

- `merge_summaries.sh`  
  Merges all run summaries into one CSV file.

---

## 3. What changed relative to the original workflow

### 3.1 Cached input series

The parameter grid varies mainly in:

- MFDFA fitting range,
- window size,
- window step,
- number of shuffles,
- and output directory.

However, the binned seismic series depends only on:

- `bin`,
- `series`,
- `magmin`.

So instead of rebuilding the same binned series for every run, this refactor computes it once and stores it in:

```text
cache_series/series_<series>_<bin>_m<magmin>.csv
```

This turns repeated catalog parsing and repeated resampling into a one-time cost.

### 3.2 Python interpreter consistency

All shell scripts now use `python3`, not `python`.

This matters on clusters because:

- `python` may not exist,
- it may point to Python 2,
- or it may refer to a different environment than the one used interactively.

### 3.3 Safe shell behavior

The old launcher had a typo in the shell flags. The refactor uses:

```bash
set -euo pipefail
```

consistently.

### 3.4 Prevent hidden oversubscription

In job-level parallel scans, NumPy, OpenBLAS, MKL, or NumExpr may try to use multiple threads **inside each job**. That creates hidden oversubscription.

For example, on a 16-thread machine:

- 12 GNU parallel jobs × 4 BLAS threads/job = 48 runnable threads,
- which causes context switching, cache pressure, and slower runtimes.

The new launcher forces:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
```

This is one of the most important practical changes.

### 3.5 No duplicate windowed recomputation

The original workflow computed windowed `h(q)` once, then later computed windowed `Δα` by running another full window loop. In the refactored code, both are extracted from the **same window pass**.

### 3.6 Lighter aggregation

The old aggregation script reopened `windowed_hq.csv`, rebuilt `Δh`, and regenerated a plot. In the refactor, the main analysis script already writes the relevant summary CSVs, and the aggregator only collects numerical summaries.

### 3.7 Timing instrumentation

Each run now records timing blocks in `summary.json`, for example:

- catalog loading,
- series building,
- preprocessing,
- ACF block,
- main MFDFA block,
- windowed block,
- shuffle block,
- total runtime.

This makes future bottleneck studies much easier.

---

## 4. Installation

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

Install GNU parallel if needed:

```bash
sudo apt install parallel
```

---

## 5. Expected external prerequisites

Before running the scan, make sure the machine already has:

- Python 3
- `bash`
- GNU `parallel`
- enough free disk space for `scan_outputs/`
- enough memory for several concurrent Python jobs

For cluster runs, also make sure that:

- the same Python environment is visible from the working directory,
- `python3` is on `PATH`,
- and the node-local or shared filesystem can handle many concurrent writes.

---

## 6. Input file expected by default

The default scripts assume the earthquake catalog is named:

```text
eq_data_earthquake_reviewed_mag4.csv
```

in the same directory as the scripts.

The catalog must contain at least:

- `time`
- `mag`

and optionally:

- `latitude`
- `longitude`

The current timestamp parser expects the format:

```text
YYYY-MM-DD HH:MM:SS
```

in UTC-style processing.

---

## 7. Running the full workflow

### Step 1: generate the cached series and jobs

```bash
bash generate_jobs.sh
```

This will:

1. create `cache_series/`
2. generate the reusable binned series
3. create `jobs.txt`

### Step 2: execute in parallel

```bash
bash run_parallel_scan.sh
```

By default:

- `THREADS=$(nproc)`
- `JOBS=THREADS/2`

You can override this manually:

```bash
JOBS=8 THREADS=16 bash run_parallel_scan.sh
```

### Step 3: merge all summaries

```bash
bash merge_summaries.sh
```

Final merged output:

```text
scan_outputs/all_parameter_scan_results.csv
```

---

## 8. Recommended job counts

These are starting points, not rigid rules.

### 16-thread laptop or workstation

Start with:

```bash
JOBS=6 bash run_parallel_scan.sh
```

or

```bash
JOBS=8 bash run_parallel_scan.sh
```

### 32-thread machine

Start with:

```bash
JOBS=12 bash run_parallel_scan.sh
```

### Why not use all threads?

Because each job also performs:

- plotting,
- CSV writing,
- JSON writing,
- pandas parsing,
- memory allocation,
- and repeated short NumPy kernels.

These workloads are not ideal for full saturation. In practice, using about **50% to 70%** of visible threads often gives better wall-clock efficiency.

---

## 9. Performance logic behind the redesign

### 9.1 Where the runtime really goes

For this type of workflow, the expensive parts are usually:

1. windowed MFDFA,
2. shuffled multifractal tests,
3. repeated detrending across many scales and windows,
4. large numbers of PNG writes,
5. repeated catalog loading and resampling.

### 9.2 Why caching helps so much

Suppose many runs share the same `(bin, series, magmin)` configuration. Then the catalog-to-series transformation is identical. Recomputing it in every job wastes time and filesystem bandwidth.

Caching makes the scan more similar to this logic:

- **once**: raw catalog → binned signal
- **many times**: binned signal → MFDFA parameter study

This is the right separation of concerns.

### 9.3 Why single-threaded jobs are better here

The scan is already embarrassingly parallel at the outer level. If each job also becomes internally multi-threaded, the machine gets oversubscribed and loses efficiency.

So the best strategy is usually:

- **many single-threaded jobs**,
- not **fewer multi-threaded jobs**.

---

## 10. Output structure

Each run directory inside `scan_outputs/` contains the normal scientific outputs plus timing metadata.

Common files include:

- `acf.csv`
- `mfdfa_Fq.csv`
- `multifractal_spectrum.csv`
- `windowed_hq.csv`
- `delta_h_h0_h5.csv`
- `delta_alpha_windowed.csv`
- `summary.json`
- `run_summary.csv`

The key new file for performance diagnosis is:

- `summary.json`

because it stores `timings_sec`.

---

## 11. Optional further accelerations

This refactor deliberately stays close to your original scientific code. If later you want even more speed, the next candidates are:

### A. Disable plots for the large scan

Add `--skip-plots` to large production scans, and regenerate plots only for selected final runs.

This often gives a visible speedup and reduces filesystem load.

### B. Reduce `n_shuffles`

Shuffles are scientifically useful, but they are expensive because each shuffle repeats a large fraction of the MFDFA pipeline.

A good workflow is:

- coarse scan with smaller `n_shuffles`,
- focused reruns with larger `n_shuffles` on the best parameter regions.

### C. Profile the windowed block separately

Because `windowed_block_sec` is reported explicitly, you can now identify whether runtime is dominated by:

- the main full-series MFDFA,
- the moving-window study,
- or the surrogate tests.

### D. Consider algorithmic acceleration later

If you need still more speed, the next serious step would be to optimize the detrending kernel itself, for example by:

- reducing repeated `polyfit` calls,
- exploiting closed-form linear detrending for `poly=1`,
- or moving the inner loops to Numba/Cython/C++.

That is a second-stage optimization. The current refactor addresses the **workflow-level bottlenecks first**, which is the correct first move.

---

## 12. Suggested GitHub usage

A good repository structure would be:

```text
your_repo/
├── seismic_parallel_refactor/
│   ├── analyze_seismic_multifractal_fast.py
│   ├── prepare_binned_series_cache.py
│   ├── aggregate_seismic_results.py
│   ├── generate_jobs.sh
│   ├── run_parallel_scan.sh
│   ├── merge_summaries.sh
│   ├── requirements.txt
│   └── README.md
└── eq_data_earthquake_reviewed_mag4.csv
```

Then from inside `seismic_parallel_refactor/` you can run the workflow directly.

---

## 13. Final recommendation

For your very first serious parallelization attempt, the most important lesson is this:

> Do not start by making every numerical kernel parallel. First separate the workflow into what is **truly variable** and what is **repeated unnecessarily**.

That is exactly what this refactor does.

