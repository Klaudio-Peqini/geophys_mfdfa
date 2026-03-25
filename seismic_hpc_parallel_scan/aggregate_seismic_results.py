#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)


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
    summary = read_json(run_dir / 'summary.json')
    env = read_json(run_dir / 'environment.json')
    status = read_json(run_dir / 'status.json')
    spec = summary.get('spectrum_summary', {})
    timings = summary.get('timings_sec', {})
    runtime = summary.get('runtime_features', {})
    window_meta = summary.get('windowed_hq', {})

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
        'h2': spec.get('h2', float('nan')),
        'h0': spec.get('h0', float('nan')),
        'spectrum_width_alpha': spec.get('spectrum_width_alpha', float('nan')),
        'alpha_at_fmax': spec.get('alpha_at_fmax', float('nan')),
        'falpha_max': spec.get('falpha_max', float('nan')),
        'window_count': window_meta.get('window_count', 0),
        'window_workers': runtime.get('window_workers', 1),
        'window_parallel_enabled': runtime.get('window_parallel_enabled', False),
        'total_runtime_sec': timings.get('total_runtime_sec', float('nan')),
        'load_cached_series_sec': timings.get('load_cached_series_sec', float('nan')),
        'load_catalog_sec': timings.get('load_catalog_sec', float('nan')),
        'build_series_sec': timings.get('build_series_sec', float('nan')),
        'preprocess_series_sec': timings.get('preprocess_series_sec', float('nan')),
        'mfdfa_main_compute_sec': timings.get('mfdfa_main_compute_sec', float('nan')),
        'windowed_compute_sec': timings.get('windowed_compute_sec', float('nan')),
        'summary_write_sec': timings.get('summary_write_sec', float('nan')),
        'hostname': env.get('hostname'),
        'python_version': env.get('python_version'),
        'cpu_count_logical': env.get('cpu_count_logical'),
        'status_success': status.get('success', False),
        'summary_json': str((run_dir / 'summary.json').resolve()),
        'environment_json': str((run_dir / 'environment.json').resolve()),
    }

    for name, filename in {
        'windowed_hq_csv': 'windowed_hq.csv',
        'delta_h_csv': 'delta_h_h0_h5.csv',
        'delta_alpha_csv': 'delta_alpha_windowed.csv',
        'multifractal_spectrum_csv': 'multifractal_spectrum.csv',
    }.items():
        row[name] = str((run_dir / filename).resolve()) if (run_dir / filename).exists() else ''

    pd.DataFrame([row]).to_csv(run_dir / 'run_summary.csv', index=False)
    with open(run_dir / 'run_complete.ok', 'w', encoding='utf-8') as f:
        f.write('ok\n')


if __name__ == '__main__':
    main()
