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
