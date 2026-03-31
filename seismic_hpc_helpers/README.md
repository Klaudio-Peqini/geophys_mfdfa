# HPC Helper Scripts for Multifractal Analysis            

This folder contains utility scripts to automate plotting and packaging of results from large-scale multifractal analyses performed on HPC systems.

## Contents

### 1. `plot_figures.py`
Generates all key plots for a single run directory:
- ACF
- Binned series
- Δα (windowed)
- Δh = h(0) − h(5)
- Multifractal spectrum
- Shuffle test summary
- Windowed h(q)

**Usage**
```bash
python3 plot_figures.py --input-dir path/to/run_folder
```

---
### 2. `run_plots_all.sh`
Batch script to execute `plot_figures.py` across all run folders.

**Usage**
```bash
bash run_plots_all.sh scan_outputs plot_figures.py
```

---
### 3. `package_multifractal_analysis_inputs.sh`
Packages all relevant outputs for selected runs into a ZIP archive for detailed analysis.

Includes:
- CSV/JSON outputs
- Generated plots

Default behavior:
- Selects only `_s1` and `_s7` runs (small-shift windows)

**Usage**
```bash
bash package_multifractal_analysis_inputs.sh
```

---

### 4. `package_delta_h_only.sh`
Creates a lightweight ZIP containing only:
```
delta_h_h0_h5.csv
```
for each selected run, preserving folder structure.

Useful for:
- fast comparative analysis
- statistical post-processing

**Usage**
```bash
bash package_delta_h_only.sh
```

---

## Typical Workflow

1. Run analysis on HPC → generate `scan_outputs/`
2. Produce plots:
```bash
bash run_plots_all.sh
```
3. Package results:
```bash
bash package_multifractal_analysis_inputs.sh
```
or (lightweight):
```bash
bash package_delta_h_only.sh
```

---

## Notes

- Scripts are designed for large batch processing (100+ runs).
- Folder naming convention `run_*` is assumed.
- Small-shift analysis focuses on `step = 1` and `step = 7`.
