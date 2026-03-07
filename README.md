# geomag-seismic-multifractal

A research-oriented Python toolkit for **multifractal analysis**, **long-range correlation analysis**, and **scaling diagnostics** applied to:

- 🌍 Geomagnetic time series (paleointensity, VADM, RPI, secular variation)
- 🌋 Seismic catalogs (event counts, magnitudes, energy release, inter-event times)
- Geophysical catalogs of various types

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

The Python script in the main branch constructs the multifractal spectrum for three case:
- Counts for M > 4 (the whole catalogue)
- Energy for M > 4 (the whole catalogue)
- M > 5 earthquakes

1) Counts per day (good first pass)
```
python analyze_seismic_multifractal.py \
  --catalog /mnt/data/eq_data_earthquake_reviewed_mag4.csv \
  --bin 1D \
  --series counts \
  --log1p \
  --maxlag 400 \
  --out results_counts_1D \
  --precomputed_corr /mnt/data/eq_self_singlepair_19680831_W3350_perwin_selflag_singlepair.csv
```
2) Energy proxy per day (often more “physical” than counts)
```   
python analyze_seismic_multifractal.py \
  --catalog /mnt/data/eq_data_earthquake_reviewed_mag4.csv \
  --bin 1D \
  --series energy \
  --log1p \
  --maxlag 400 \
  --out results_energy_1D
```
3) Restrict to a time interval or magnitude band
```
python analyze_seismic_multifractal.py \
  --catalog /mnt/data/eq_data_earthquake_reviewed_mag4.csv \
  --bin 1D \
  --series counts \
  --magmin 5.0 \
  --start 1970-01-01 \
  --end 2020-12-31 \
  --log1p \
  --out results_mag5_1970_2020
```

## What you’ll get in the output folder

```
- binned_series.csv, binned_series.png
- acf.csv, acf.png
- mfdfa_Fq.csv, mfdfa_Fq.png
- multifractal_spectrum.csv, spectrum.png
- summary.json (parameters + key spectrum summaries)
```

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
