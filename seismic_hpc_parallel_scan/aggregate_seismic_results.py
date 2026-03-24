#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    p = argparse.ArgumentParser(description='Collect one-row summaries from each analysis run.')
    p.add_argument('--run_dir', required=True)
    p.add_argument('--run_id', required=True)
    p.add_argument('--bin', required=True)
    p.add_argument('--series', required=True)
    p.add_argument('--magmin', required=True)
    p.add_argument('--window_size_bins', required=True, type=int)
    p.add_argument('--window_step_bins', required=True, type=int)
    p.add_argument('--qmin', required=True, type=float)
    p.add_argument('--qmax', required=True, type=float)
    p.add_argument('--qstep', required=True, type=float)
    p.add_argument('--fit_smin', required=True, type=int)
    p.add_argument('--fit_smax', required=True, type=int)
    args = p.parse_args()

    run_dir = Path(args.run_dir)
    summary_path = run_dir / 'summary.json'
    delta_h_csv = run_dir / 'delta_h_h0_h5.csv'
    delta_alpha_csv = run_dir / 'delta_alpha_windowed.csv'
    run_summary_csv = run_dir / 'run_summary.csv'

    h2 = h0 = width = alpha_at_fmax = falpha_max = float('nan')
    timings = {}
    if summary_path.exists():
        with open(summary_path, 'r', encoding='utf-8') as f:
            summary = json.load(f)
        spec = summary.get('spectrum_summary', {})
        h2 = spec.get('h2', float('nan'))
        h0 = spec.get('h0', float('nan'))
        width = spec.get('spectrum_width_alpha', float('nan'))
        alpha_at_fmax = spec.get('alpha_at_fmax', float('nan'))
        falpha_max = spec.get('falpha_max', float('nan'))
        timings = summary.get('timings_sec', {})

    delta_h_mean = delta_h_std = delta_h_min = delta_h_max = float('nan')
    delta_alpha_mean = delta_alpha_std = delta_alpha_min = delta_alpha_max = float('nan')
    n_windows = 0

    if delta_h_csv.exists():
        df = pd.read_csv(delta_h_csv)
        if 'delta_h' in df.columns:
            delta_h_mean = float(df['delta_h'].mean())
            delta_h_std = float(df['delta_h'].std())
            delta_h_min = float(df['delta_h'].min())
            delta_h_max = float(df['delta_h'].max())
            n_windows = int(len(df))

    if delta_alpha_csv.exists():
        da = pd.read_csv(delta_alpha_csv)
        if 'delta_alpha' in da.columns:
            delta_alpha_mean = float(da['delta_alpha'].mean())
            delta_alpha_std = float(da['delta_alpha'].std())
            delta_alpha_min = float(da['delta_alpha'].min())
            delta_alpha_max = float(da['delta_alpha'].max())

    row = {
        'run_id': args.run_id,
        'bin': args.bin,
        'series': args.series,
        'magmin': args.magmin,
        'window_size_bins': args.window_size_bins,
        'window_step_bins': args.window_step_bins,
        'qmin': args.qmin,
        'qmax': args.qmax,
        'qstep': args.qstep,
        'fit_smin': args.fit_smin,
        'fit_smax': args.fit_smax,
        'h2': h2,
        'h0': h0,
        'spectrum_width_alpha': width,
        'alpha_at_fmax': alpha_at_fmax,
        'falpha_max': falpha_max,
        'delta_h_mean': delta_h_mean,
        'delta_h_std': delta_h_std,
        'delta_h_min': delta_h_min,
        'delta_h_max': delta_h_max,
        'delta_alpha_mean': delta_alpha_mean,
        'delta_alpha_std': delta_alpha_std,
        'delta_alpha_min': delta_alpha_min,
        'delta_alpha_max': delta_alpha_max,
        'n_windows': n_windows,
        'summary_json': str(summary_path),
        'delta_h_csv': str(delta_h_csv),
        'delta_alpha_csv': str(delta_alpha_csv),
        'load_cached_series_sec': timings.get('load_cached_series_sec', float('nan')),
        'load_catalog_sec': timings.get('load_catalog_sec', float('nan')),
        'filter_catalog_sec': timings.get('filter_catalog_sec', float('nan')),
        'build_series_sec': timings.get('build_series_sec', float('nan')),
        'preprocess_series_sec': timings.get('preprocess_series_sec', float('nan')),
        'acf_block_sec': timings.get('acf_block_sec', float('nan')),
        'mfdfa_main_sec': timings.get('mfdfa_main_sec', float('nan')),
        'spectrum_block_sec': timings.get('spectrum_block_sec', float('nan')),
        'windowed_block_sec': timings.get('windowed_block_sec', float('nan')),
        'shuffle_test_sec': timings.get('shuffle_test_sec', float('nan')),
        'phase_surrogate_test_sec': timings.get('phase_surrogate_test_sec', float('nan')),
        'total_runtime_sec': timings.get('total_runtime_sec', float('nan')),
    }
    pd.DataFrame([row]).to_csv(run_summary_csv, index=False)


if __name__ == '__main__':
    main()
