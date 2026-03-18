# seismic_multifractal_scan

A GitHub-ready subfolder for systematic seismic multifractal analysis.

## Contents

- `analyze_seismic_multifractal.py`  
  Main seismic multifractal analysis script.
- `run_seismic_parameter_scan.sh`  
  Bash driver for parameter-space analysis.
- `aggregate_seismic_results.py`  
  Helper that extracts quantitative outputs from each run and appends them to a master CSV.
- `requirements.txt`  
  Minimal Python dependencies.

## What the pipeline does

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

## Required input

Place your catalog CSV in the same folder, for example:

- `eq_data_earthquake_reviewed_mag4.csv`

## Quick start

```bash
chmod +x run_seismic_parameter_scan.sh
./run_seismic_parameter_scan.sh
```

## Expected final output

Main aggregated CSV:

```bash
scan_outputs/all_parameter_scan_results.csv
```

This CSV includes:
- binning choice
- series type
- magnitude threshold
- window length and step
- fitting range
- `h0`, `h2`
- spectrum width `Δα`
- summary statistics of `Δh(t)=h(0,t)-h(5,t)`
- paths to the main output files for each run

## Suggested GitHub location

```text
geomag-seismic-multifractal/
└── seismic_multifractal_scan/
```
