# geophysics-multifractal                               

A research-oriented Python toolkit for **multifractal analysis**, **long-range correlation analysis**, and **scaling diagnostics** applied to:

- 🌍 Geomagnetic time series (paleointensity, VADM, RPI, secular variation)
- 🌋 Seismic catalogs (event counts, magnitudes, energy release, inter-event times)
- <img width="20" height="20" alt="image" src="https://github.com/user-attachments/assets/27de1aba-4deb-4a60-8ac7-1ca89627d0c7" /> Geophysical catalogs of various types

The package implements:

- Multifractal Detrended Fluctuation Analysis (MFDFA)
- Structure functions and scaling exponents
- Multifractal spectrum via Legendre transform
- Cross-MFDFA (MF-DCCA style)
- Autocorrelation and cross-correlation (FFT + direct)
- Sliding-window multifractal analysis
- Surrogate testing (phase randomization & shuffle)
- Publication-ready plotting utilities

---

# Scientific Motivation

Many natural systems exhibit:

- Long-range correlations
- Scale-free clustering
- Intermittent bursts
- Heavy-tailed fluctuations
- Multiscale variability

Such behavior is observed in:

- Geodynamo dynamics (outer core convection, reversals)
- Earthquake triggering and stress redistribution
- Turbulent and cascading systems

This package provides tools to quantify:

- Memory effects
- Scaling laws
- Intermittency
- Multifractal complexity

---

# Theoretical Background

## 1. Autocorrelation Function

The autocorrelation function:

$C(τ)$

quantifies memory across time lags.

Decay behavior distinguishes:

- Exponential decay → short memory
- Power-law decay → long-range correlation

Relation to Hurst exponent:

$C(τ) \sim τ^{2H - 2}$

---

## 2. MFDFA (Multifractal Detrended Fluctuation Analysis)

The algorithm computes:

$F_q(s)$

and extracts generalized Hurst exponents:

$h(q)$

Special cases:

- $h(2)$ → classical Hurst exponent
- Constant $h(q)$ → monofractal
- Nonlinear $h(q)$ → multifractal

---

## 3. Multifractal Spectrum

From:

$τ(q) = qh(q) - 1$

we derive:

$α = \frac{dτ}{dq}$ \
$f(α) = qα - τ(q)$

Key measurable quantities:

- Spectrum width $Δα$ → strength of multifractality
- Peak location $α_0$ → dominant scaling exponent
- Asymmetry → dominance of large vs small fluctuations

---

# Why pandas?

Geomagnetic and seismic datasets often contain:

- Irregular sampling
- Event catalogs
- Time stamps or age axes
- Missing values

pandas ensures:

- Robust parsing
- Resampling of event catalogs
- Time alignment
- Safe handling of NaNs

Core numerical algorithms use numpy for performance and transparency.

---

# Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .
```

Optional:

```bash
pip install -r requirements.txt
```

---

# Repository Structure

```
geomag-seismic-multifractal/
├── data/
│   ├── geomagnetic/
│   └── seismic/
├── src/gsmf/
│   ├── io/
│   ├── preprocessing/
│   ├── multifractal/
│   ├── analysis/
│   └── plotting/
├── notebooks/
├── examples/
├── results/
│   ├── figures/
│   └── tables/
├── tests/
└── scripts/
```

---

# Quick Start – Seismic Example

```python
import numpy as np
from gsmf.io import load_earthquake_catalog, events_to_counts
from gsmf.analysis import autocorrelation
from gsmf.multifractal import mfdfa, legendre_spectrum

cat = load_earthquake_catalog("data/seismic/eq_data_earthquake_reviewed_mag4.csv")

# Convert event catalog to daily counts
x = events_to_counts(cat, freq="1D").values

# Autocorrelation
acf = autocorrelation(x, max_lag=365, method="fft")

# Multifractal analysis
scales = np.unique(np.logspace(np.log10(16), np.log10(2048), 20).astype(int))
res = mfdfa(x, scales=scales, qs=np.linspace(-5,5,21), poly_order=2)

# Multifractal spectrum
spec = legendre_spectrum(res.qs, res.hq)
```

---

# Geomagnetic Workflow Example

Geomagnetic series are typically provided as:

- Age (kyr)
- RPI or VADM

The recommended workflow:

1. Interpolate onto uniform age grid (if needed)
2. Compute autocorrelation
3. Run MFDFA
4. Compute spectrum
5. Perform surrogate testing
6. Run windowed MFDFA for nonstationarity

Example script provided in:

scripts/run_geomag_analysis.py

---

# CLI Usage

```bash
gsmf-quick   --catalog data/seismic/eq_data_earthquake_reviewed_mag4.csv   --freq 1D   --max-lag 365   --out-prefix results/tables/quick
```

Outputs:

- ACF table
- h(q) table
- Multifractal spectrum table

---

# Plotting

All plotting helpers use matplotlib.

Example ACF plot:

```python
import matplotlib.pyplot as plt
from gsmf.analysis import autocorrelation
from gsmf.plotting import plot_acf

acf = autocorrelation(x, max_lag=18000, method="fft")

plot_acf(
    acf,
    style="m-",
    marker="+",
    markersize=4,
    linewidth=1.0
)

plt.show()
```

Available plotting functions:

- plot_acf()
- plot_ccf()
- plot_mfdfa_Fq()
- plot_hq()
- plot_spectrum()
- plot_structure_functions()

---

# Meaningful Scientific Outputs

## From Autocorrelation

- Correlation time
- Power-law decay exponent
- Presence of long-range memory

## From MFDFA

- Hurst exponent h(2)
- Multifractal strength Δα
- Scaling regime consistency
- Intermittency diagnostics

## From Surrogate Testing

Distinguish:

- Correlation-driven multifractality
- Distribution-driven multifractality

---

# Applications

## Geomagnetic Data

- Paleointensity variability
- Reversal precursors
- Core turbulence signatures
- Long-memory geodynamo behavior

## Seismic Data

- Aftershock clustering
- Stress transfer cascades
- SOC-like scaling
- Pre-large-event intermittency

## Examples
In the folder ```examples```, several Python scripts are shown that compute the relevant quantities and enable drawing meaningful conclusions.
In the same folder is a mock file. It allows the generation of mock data to apply the various functions in the repo.
The mock data have a structure that mimics the real data from geophysical sources.

---

# How to run the analysis Python file in the main branch?

The Python script in the main branch constructs the multifractal spectrum,
as well as calculates other quantities, for multiple cases:

1) Counts per day
```
python analyze_seismic_multifractal.py \
  --catalog /path/to/file/eq_data_earthquake_reviewed_mag4.csv \
  --out results_counts_1D \
  --bin 1D \
  --series counts \
  --log1p \
  --qmin -30 --qmax 30 --qstep 1 \
  --fit_smin 121 \
  --acf_shuffle_band \
  --n_shuffles 20
```
2) Temporal evolution of selected h(q)
```   
python analyze_seismic_multifractal.py \
  --catalog /path/to/file/eq_data_earthquake_reviewed_mag4.csv \
  --out results_windowed \
  --bin 1D \
  --series counts \
  --log1p \
  --window_size_bins 3650 \
  --window_step_bins 365 \
  --window_q_values 0 2 5
```
3) For a given subperiod
```
python analyze_seismic_multifractal.py \
  --catalog /file/to/path/eq_data_earthquake_reviewed_mag4.csv \
  --out results_subperiod \
  --bin 1D \
  --series energy \
  --log1p \
  --part_start 1980-01-01 \
  --part_end 2000-12-31
```
4) Energy per 1 day
```
python analyze_seismic_multifractal.py \
  --catalog /path/to/file/eq_data_earthquake_reviewed_mag4.csv \
  --out results_energy_1D \
  --bin 1D \
  --series energy \
  --log1p \
  --qmin -30 --qmax 30 --qstep 1 \
  --acf_shuffle_band \
  --n_shuffles 20
```
5) Energy per 3 days
```
python analyze_seismic_multifractal.py \
  --catalog /path/to/file/eq_data_earthquake_reviewed_mag4.csv \
  --out results_energy_3D \
  --bin 3D \
  --series energy \
  --log1p \
  --qmin -30 --qmax 30 --qstep 1 \
  --acf_shuffle_band \
  --n_shuffles 20
```
6) Energy per 7 days
```
python analyze_seismic_multifractal.py \
  --catalog /path/to/file/eq_data_earthquake_reviewed_mag4.csv \
  --out results_energy_7D \
  --bin 7D \
  --series energy \
  --log1p \
  --qmin -30 --qmax 30 --qstep 1 \
  --acf_shuffle_band \
  --n_shuffles 20
```
7) Energy for magnitudes greater than 5.0
```
python analyze_seismic_multifractal.py \
  --catalog /path/to/file/eq_data_earthquake_reviewed_mag4.csv \
  --out results_energy_maggt5_1D \
  --bin 1D \
  --series energy \
  --magmin 5.0 \
  --log1p \
  --qmin -30 --qmax 30 --qstep 1 \
  --acf_shuffle_band \
  --n_shuffles 20
```
8) Energy for magnitudes greater than 5.5, binned every 7 days
```
python analyze_seismic_multifractal.py \
  --catalog /path/to/file/eq_data_earthquake_reviewed_mag4.csv \
  --out results_energy_maggt5p5_7D \
  --bin 7D \
  --series energy \
  --magmin 5.5 \
  --log1p \
  --qmin -30 --qmax 30 --qstep 1 \
  --acf_shuffle_band \
  --n_shuffles 20
```
9) Counts for magnitudes greater than a threshold
```
python analyze_seismic_multifractal.py \
  --catalog /path/to/file/eq_data_earthquake_reviewed_mag4.csv \
  --out results_counts_maggt5_1D \
  --bin 1D \
  --series counts \
  --magmin 5.0 \
  --log1p \
  --qmin -30 --qmax 30 --qstep 1 \
  --acf_shuffle_band \
  --n_shuffles 20
```
10) Mean magnitude per 3 days for events above 4.5
```
python analyze_seismic_multifractal.py \
  --catalog /path/to/file/eq_data_earthquake_reviewed_mag4.csv \
  --out results_meanmag_maggt4p5_3D \
  --bin 3D \
  --series meanmag \
  --magmin 4.5 \
  --qmin -30 --qmax 30 --qstep 1 \
  --acf_shuffle_band \
  --n_shuffles 20
```
11) Restrict to a specific time interval as well
```
python analyze_seismic_multifractal.py \
  --catalog /path/to/file/eq_data_earthquake_reviewed_mag4.csv \
  --out results_energy_1D_1980_2020_maggt5 \
  --bin 1D \
  --series energy \
  --magmin 5.0 \
  --start 1980-01-01 \
  --end 2020-12-31 \
  --log1p \
  --qmin -30 --qmax 30 --qstep 1 \
  --acf_shuffle_band \
  --n_shuffles 20
```
12) Add a selected fitting range for the MFDFA lines, if you already know the physically meaningful scaling region:
```
python analyze_seismic_multifractal.py \
  --catalog /path/to/file/eq_data_earthquake_reviewed_mag4.csv \
  --out results_energy_1D_fit \
  --bin 1D \
  --series energy \
  --log1p \
  --qmin -30 --qmax 30 --qstep 1 \
  --fit_smin 121 \
  --fit_smax 4113 \
  --acf_shuffle_band \
  --n_shuffles 20
```
13) Temporal evolution of $h(q)$ for energy. Example: 10-year windows stepped by 1 year for daily bins:
```
python analyze_seismic_multifractal.py \
  --catalog /path/to/file/eq_data_earthquake_reviewed_mag4.csv \
  --out results_energy_windowed \
  --bin 1D \
  --series energy \
  --log1p \
  --qmin -30 --qmax 30 --qstep 1 \
  --window_size_bins 3650 \
  --window_step_bins 365 \
  --window_q_values 0 2 5 10
```
14) Full analysis of all relevant multifractal quantities, like $Δh(t)$, or $Δα(t)$
```
python analyze_seismic_multifractal.py \
  --catalog /file/to/path/eq_data_earthquake_reviewed_mag4.csv \
  --out results_seismic_full \
  --bin 1D \
  --series energy \
  --log1p \
  --qmin -10 --qmax 10 --qstep 1 \
  --fit_smin 121 \
  --fit_smax 4113 \
  --maxlag 365 \
  --acf_shuffle_band \
  --acf_shuffle_count 200 \
  --n_shuffles 50 \
  --window_size_bins 3650 \
  --window_step_bins 365 \
  --window_q_values 0 1 2 3 4 5
```

### Take notice:
- It is mathematically possible to run for $q∈[−30,30]$, but in seismic count series the negative side can become very unstable because tiny local variances dominate.
- For publication-quality results, a conservative range such as $[−5,5]$, or $[−10,10]$ is reccomended.
- Use --log1p almost always for --series energy.
- For sparse subsets such as ```--magmin 5.5``` or higher, try ```--bin 3D``` or ```--bin 7D``` as well, because daily bins can become too sparse.

---

## Scanning parameters space

A GitHub-ready subfolder for systematic seismic multifractal analysis.

```
geomag-seismic-multifractal/
└── seismic_multifractal_scan/
```

### Contents

- `analyze_seismic_multifractal.py`
  Main seismic multifractal analysis script.
- `run_seismic_parameter_scan.sh`
  Bash driver for parameter-space analysis.
- `aggregate_seismic_results.py`
  A helper that extracts quantitative outputs from each run and appends them to a master CSV.
- `requirements.txt`
  Minimal Python dependencies.

### What the pipeline does

For each selected parameter combination, it:
1. builds a binned seismic series from the catalog
2. runs ACF + MFDFA
3. computes windowed `h(q,t)`
4. saves:
   - `multifractality_strength_h0_h5.png`
   - `delta_alpha_windowed.png`
5. derives:
   - `delta_h_h0_h5.csv`
   - `delta_h_h0_h5.png`
6. appends a single summary row to:
   - `scan_outputs/all_parameter_scan_results.csv`

### Required input

Place your catalog CSV in the same folder, for example:

- `eq_data_earthquake_reviewed_mag4.csv`

### Quick start

```bash
source run_seismic_parameter_scan.sh
```

---

## Parallel scanning parameters space

Parallel parameter-space scan utilities for seismic multifractal analysis.

```
geomag-seismic-multifractal/
└── seismic_multifractal_parallel_scan/
```

### Files

In the respective README.md file are shown details regarding
the content of the folder and the mode of use of each.

### Output:

- delta_h_h0_h5.csv
- delta_h_h0_h5.png
- one run_summary.csv per run
- merged .csv file

### Recommended usage

```
source run_parallel_scan.sh # run all the created jobs
source merge_summaries.sh # merges all the output files into one .csv
```

### recommended for parallel run

For a 16-thread machine, start with:

```
JOBS=8
```

inside run_parallel_scan.sh. Then, if the machine remains comfortable
in RAM and I/O usage, try 10 or 12. If the scripts are run in a cluster,
make use of the threads appropriately.

---

## What you’ll get in the output folder

```
- binned_series.csv, binned_series.png
- acf.csv, acf.png
- mfdfa_Fq.csv, mfdfa_Fq.png
- multifractal_spectrum.csv, spectrum.png
- summary.json (parameters + key spectrum summaries)
```

---
## seismic-parallel-hpc

A refactored workflow for running large multifractal seismic parameter scans efficiently on a laptop, workstation, or HPC node.
This folder is designed to be copied directly into a GitHub repository. Its purpose is not only to run the scan, but also to make the **parallelization logic explicit and reproducible**.

### Why this refactor was needed

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

For more details, please refer to the README.md file in the folder.

---

# Reproducibility & Philosophy

This project emphasizes:

- Transparent algorithms
- Minimal dependencies
- Explicit scaling diagnostics
- Clear separation between I/O and analysis
- Reproducible outputs (CSV tables + figures)

The goal is not only computation, but interpretation and scientific clarity.

---

# Testing

```bash
pytest
```

Basic smoke tests validate:

- ACF correctness
- MFDFA regression stability
- Spectrum consistency

---

# Citation

See CITATION.cff.

If this package contributes to your research, please cite appropriately.

---

# Future Roadmap

- Automatic scaling-region detection
- Bootstrap confidence intervals for spectrum width
- ETAS-based synthetic seismic benchmarks
- Geodynamo simulation comparisons
- Advanced statistical hypothesis testing

---

# License

MIT License.
