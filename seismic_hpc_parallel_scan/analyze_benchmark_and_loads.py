#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def summarize_joblog(path: Path) -> dict:
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep='\t')
    out = {'joblog_rows': len(df)}
    if 'Runtime' in df.columns:
        out['job_runtime_mean_sec'] = float(df['Runtime'].mean())
        out['job_runtime_std_sec'] = float(df['Runtime'].std())
        out['job_runtime_max_sec'] = float(df['Runtime'].max())
    if 'Exitval' in df.columns:
        out['failed_jobs'] = int((df['Exitval'] != 0).sum())
    return out


def summarize_load(path: Path) -> tuple[dict, pd.DataFrame | None]:
    if not path.exists():
        return {}, None
    df = pd.read_csv(path)
    out = {}
    if not df.empty:
        out['cpu_usage_mean_percent'] = float(df['usage_percent'].mean())
        out['cpu_usage_std_percent'] = float(df['usage_percent'].std())
        out['cpu_usage_p95_percent'] = float(df['usage_percent'].quantile(0.95))
        per_cpu = df.groupby('cpu_id')['usage_percent'].agg(['mean', 'std', 'median']).reset_index()
        return out, per_cpu
    return out, None


def make_scaling_plot(df: pd.DataFrame, out_png: Path) -> None:
    plt.figure()
    plt.plot(df['jobs_standard'], df['wall_clock_sec'], marker='o')
    plt.xlabel('Standard concurrent jobs')
    plt.ylabel('Wall-clock time (s)')
    plt.title('Parallel scaling benchmark')
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def make_speedup_plot(df: pd.DataFrame, out_png: Path) -> None:
    base = float(df['wall_clock_sec'].iloc[0])
    speedup = base / df['wall_clock_sec']
    efficiency = speedup / df['jobs_standard']
    plt.figure()
    plt.plot(df['jobs_standard'], speedup, marker='o', label='Speedup')
    plt.plot(df['jobs_standard'], efficiency, marker='s', label='Efficiency')
    plt.xlabel('Standard concurrent jobs')
    plt.ylabel('Relative metric')
    plt.title('Speedup and efficiency')
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def make_cpu_pdf_plot(load_df: pd.DataFrame, out_png: Path) -> None:
    plt.figure(figsize=(10, 6))
    for cpu_id, grp in load_df.groupby('cpu_id'):
        vals = grp['usage_percent'].dropna().to_numpy()
        if vals.size == 0:
            continue
        bins = np.linspace(0, 100, 31)
        hist, edges = np.histogram(vals, bins=bins, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        plt.plot(centers, hist, lw=1.0, alpha=0.6, label=f'cpu{cpu_id}')
    plt.xlabel('CPU utilization (%)')
    plt.ylabel('Probability density')
    plt.title('Per-thread CPU-load probability density')
    if load_df['cpu_id'].nunique() <= 16:
        plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def main() -> None:
    p = argparse.ArgumentParser(description='Analyze benchmark runs and sampled CPU loads.')
    p.add_argument('--bench-dir', required=True)
    args = p.parse_args()

    bench_dir = Path(args.bench_dir)
    rows = []
    combined_loads = []

    for run_dir in sorted(bench_dir.glob('jobs_*')):
        jobs = int(run_dir.name.split('_')[-1])
        meta_path = run_dir / 'run_meta.csv'
        if not meta_path.exists():
            continue
        meta = pd.read_csv(meta_path).iloc[0].to_dict()
        row = {'jobs': jobs, 'wall_clock_sec': float(meta['wall_clock_sec'])}
        elapsed_path = run_dir / 'elapsed_seconds.txt'
        if elapsed_path.exists():
            row['elapsed_time_cmd_sec'] = float(elapsed_path.read_text().strip())
        row.update(summarize_joblog(run_dir / 'parallel_joblog.txt'))
        load_summary, per_cpu = summarize_load(run_dir / 'cpu_load_samples.csv')
        row.update(load_summary)
        rows.append(row)

        if per_cpu is not None:
            per_cpu['jobs'] = jobs
            per_cpu.to_csv(run_dir / 'cpu_load_per_thread_summary.csv', index=False)
        load_path = run_dir / 'cpu_load_samples.csv'
        if load_path.exists():
            df_load = pd.read_csv(load_path)
            df_load['jobs'] = jobs
            combined_loads.append(df_load)

    if not rows:
        raise SystemExit('No benchmark runs found.')

    summary = pd.DataFrame(rows).sort_values('jobs').reset_index(drop=True)
    summary['speedup_vs_first'] = summary['wall_clock_sec'].iloc[0] / summary['wall_clock_sec']
    summary['efficiency_vs_first'] = summary['speedup_vs_first'] / summary['jobs']
    summary.to_csv(bench_dir / 'benchmark_summary.csv', index=False)

    make_scaling_plot(summary, bench_dir / 'benchmark_scaling.png')
    make_speedup_plot(summary, bench_dir / 'benchmark_speedup_efficiency.png')

    if combined_loads:
        loads = pd.concat(combined_loads, ignore_index=True)
        loads.to_csv(bench_dir / 'combined_cpu_load_samples.csv', index=False)
        best_jobs = int(summary.loc[summary['wall_clock_sec'].idxmin(), 'jobs'])
        best_loads = loads[loads['jobs'] == best_jobs].copy()
        make_cpu_pdf_plot(best_loads, bench_dir / f'cpu_load_pdf_jobs_{best_jobs}.png')

    print(f'Wrote {bench_dir / "benchmark_summary.csv"}')


if __name__ == '__main__':
    main()
