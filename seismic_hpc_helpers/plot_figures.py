
#!/usr/bin/env python3
"""
plot_figures.py

Create the standard plots for a single MFDFA seismic run.

Expected input files in --input-dir:
    acf.csv
    binned_series.csv
    delta_alpha_windowed.csv
    delta_h_h0_h5.csv
    multifractal_spectrum.csv
    shuffle_test_summary.csv
    windowed_hq.csv
Optional:
    summary.json
    run_summary.csv
    mfdfa_Fq.csv

Outputs:
    PNG files in --output-dir
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import AutoDateLocator, ConciseDateFormatter


def _read_csv(path: Path, required: bool = True) -> Optional[pd.DataFrame]:
    if path.exists():
        return pd.read_csv(path)
    if required:
        raise FileNotFoundError(f"Missing required file: {path}")
    return None


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _parse_time_column(df: pd.DataFrame, column: str) -> pd.Series:
    if column not in df.columns:
        raise KeyError(f"Column '{column}' not found in dataframe.")
    return pd.to_datetime(df[column], utc=True, errors="coerce").dt.tz_convert(None)


def _setup_time_axis(ax):
    locator = AutoDateLocator()
    formatter = ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)


def _apply_common_style(ax, title: str, xlabel: str, ylabel: str):
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)


def _safe_run_label(summary: dict, run_summary: Optional[pd.DataFrame]) -> str:
    if run_summary is not None and "run_id" in run_summary.columns and len(run_summary) > 0:
        return str(run_summary["run_id"].iloc[0])
    if summary.get("out"):
        return Path(summary["out"]).name
    parts = []
    if summary.get("series"):
        parts.append(str(summary["series"]))
    if summary.get("bin"):
        parts.append(str(summary["bin"]))
    return "_".join(parts) if parts else "run"


def _savefig(fig, outpath: Path):
    fig.tight_layout()
    fig.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {outpath}")


def plot_acf(acf: pd.DataFrame, outdir: Path, run_label: str):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(acf["lag_bins"], acf["acf"], linewidth=1.6, label="ACF")

    if "white_noise_95_band" in acf.columns:
        band = pd.to_numeric(acf["white_noise_95_band"], errors="coerce")
        if band.notna().any():
            band_value = float(band.dropna().iloc[0])
            ax.axhline(+band_value, linestyle="--", linewidth=1.2, label="95% white-noise band")
            ax.axhline(-band_value, linestyle="--", linewidth=1.2)

    if "acf_physically_relevant_analytic" in acf.columns:
        mask = acf["acf_physically_relevant_analytic"].astype(str).str.lower().isin(["true", "1", "yes"])
        if mask.any():
            ax.scatter(
                acf.loc[mask, "lag_bins"],
                acf.loc[mask, "acf"],
                s=18,
                marker="o",
                label="physically relevant"
            )

    _apply_common_style(ax, f"Autocorrelation function — {run_label}", "Lag (bins)", "ACF")
    ax.legend()
    _savefig(fig, outdir / "plot_acf.png")


def plot_binned_series(series_df: pd.DataFrame, outdir: Path, run_label: str):
    t = _parse_time_column(series_df, "time")
    y = pd.to_numeric(series_df["value"], errors="coerce")

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(t, y, linewidth=0.9)
    _setup_time_axis(ax)
    _apply_common_style(ax, f"Binned series — {run_label}", "Time", "Value")
    _savefig(fig, outdir / "plot_binned_series.png")


def plot_delta_alpha(delta_alpha: pd.DataFrame, outdir: Path, run_label: str):
    t = _parse_time_column(delta_alpha, "window_mid_time")
    y = pd.to_numeric(delta_alpha["delta_alpha"], errors="coerce")

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(t, y, linewidth=1.0, label=r"$\Delta \alpha$")
    ax.axhline(np.nanmean(y), linestyle="--", linewidth=1.2, label=f"mean = {np.nanmean(y):.3f}")
    _setup_time_axis(ax)
    _apply_common_style(ax, f"Windowed multifractal width Δα — {run_label}", "Window midpoint", r"$\Delta \alpha$")
    ax.legend()
    _savefig(fig, outdir / "plot_delta_alpha_windowed.png")


def plot_delta_h(delta_h_df: pd.DataFrame, outdir: Path, run_label: str):
    t = _parse_time_column(delta_h_df, "window_mid_time")
    y = pd.to_numeric(delta_h_df["delta_h"], errors="coerce")

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(t, y, linewidth=1.0, label=r"$\Delta h$")
    ax.axhline(np.nanmean(y), linestyle="--", linewidth=1.2, label=f"mean = {np.nanmean(y):.3f}")
    _setup_time_axis(ax)
    _apply_common_style(ax, f"Windowed Δh (h(0) − h(5)) — {run_label}", "Window midpoint", r"$\Delta h$")
    ax.legend()
    _savefig(fig, outdir / "plot_delta_h_h0_h5.png")


def plot_multifractal_spectrum(spec: pd.DataFrame, outdir: Path, run_label: str):
    alpha = pd.to_numeric(spec["alpha"], errors="coerce")
    falpha = pd.to_numeric(spec["falpha"], errors="coerce")
    q = pd.to_numeric(spec["q"], errors="coerce") if "q" in spec.columns else None

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(alpha, falpha, marker="o", linewidth=1.4)
    if q is not None:
        idx_max = np.nanargmax(falpha.to_numpy())
        ax.scatter([alpha.iloc[idx_max]], [falpha.iloc[idx_max]], s=40, label=f"max at q={q.iloc[idx_max]:g}")
    _apply_common_style(ax, f"Multifractal spectrum — {run_label}", r"$\alpha$", r"$f(\alpha)$")
    if q is not None:
        ax.legend()
    _savefig(fig, outdir / "plot_multifractal_spectrum.png")


def plot_shuffle_summary(shuffle: pd.DataFrame, summary: dict, outdir: Path, run_label: str):
    width = pd.to_numeric(shuffle["spectrum_width_alpha"], errors="coerce").dropna()
    real_width = None
    if "spectrum_summary" in summary and "spectrum_width_alpha" in summary["spectrum_summary"]:
        real_width = float(summary["spectrum_summary"]["spectrum_width_alpha"])
    elif "shuffle_test" in summary and "mean_shuffled_width" in summary["shuffle_test"]:
        # real width not present there, so leave as None if absent
        pass

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(width, bins=min(15, max(5, int(math.sqrt(max(len(width), 1))))), alpha=0.8, label="Shuffled widths")
    if real_width is not None:
        ax.axvline(real_width, linestyle="--", linewidth=2.0, label=f"Real width = {real_width:.3f}")
    if "shuffle_test" in summary and "mean_shuffled_width" in summary["shuffle_test"]:
        ax.axvline(summary["shuffle_test"]["mean_shuffled_width"], linestyle=":", linewidth=1.5,
                   label=f"Shuffle mean = {summary['shuffle_test']['mean_shuffled_width']:.3f}")
    p_value = summary.get("shuffle_test", {}).get("p_value_width_ge_real", None)
    title = f"Shuffle-test summary — {run_label}"
    if p_value is not None:
        title += f" (p = {p_value:.3g})"
    _apply_common_style(ax, title, r"Shuffled spectrum width $\Delta \alpha$", "Count")
    ax.legend()
    _savefig(fig, outdir / "plot_shuffle_test_summary.png")


def _extract_hq_columns(windowed_hq: pd.DataFrame) -> list[tuple[float, str]]:
    cols = []
    for col in windowed_hq.columns:
        if not col.startswith("h_q"):
            continue
        suffix = col[3:]  # after h_q
        try:
            qval = float(suffix)
        except ValueError:
            continue
        cols.append((qval, col))
    cols.sort(key=lambda x: x[0])
    return cols


def plot_windowed_hq(windowed_hq: pd.DataFrame, outdir: Path, run_label: str,
                     q_values: Optional[Iterable[float]] = None, max_lines: int = 7):
    t = _parse_time_column(windowed_hq, "window_mid_time")
    all_cols = _extract_hq_columns(windowed_hq)
    if not all_cols:
        raise ValueError("No h_q* columns found in windowed_hq.csv")

    available_qs = [q for q, _ in all_cols]
    if q_values is None:
        # sensible default subset across the q range
        preferred = [-10, -5, 0, 2, 5, 10]
        selected = []
        for pq in preferred:
            nearest = min(all_cols, key=lambda x: abs(x[0] - pq))
            if nearest not in selected:
                selected.append(nearest)
        selected = selected[:max_lines]
    else:
        selected = []
        for q_target in q_values:
            nearest = min(all_cols, key=lambda x: abs(x[0] - q_target))
            if nearest not in selected:
                selected.append(nearest)
        selected = selected[:max_lines]

    fig, ax = plt.subplots(figsize=(12, 6))
    for qval, col in selected:
        y = pd.to_numeric(windowed_hq[col], errors="coerce")
        ax.plot(t, y, linewidth=1.0, label=f"h(q={qval:g})")

    _setup_time_axis(ax)
    _apply_common_style(ax, f"Windowed generalized Hurst exponents — {run_label}",
                        "Window midpoint", "h(q)")
    ax.legend(ncol=2)
    _savefig(fig, outdir / "plot_windowed_hq.png")


def plot_mfdfa_scaling(mfdfa_fq: Optional[pd.DataFrame], outdir: Path, run_label: str,
                       q_values: Optional[Iterable[float]] = None, max_lines: int = 6):
    if mfdfa_fq is None:
        return

    x = pd.to_numeric(mfdfa_fq["scale"], errors="coerce")
    cols = []
    for col in mfdfa_fq.columns:
        if col.startswith("Fq_q"):
            try:
                q = float(col.replace("Fq_q", ""))
            except ValueError:
                continue
            cols.append((q, col))
    cols.sort(key=lambda x: x[0])

    if q_values is None:
        preferred = [-10, -5, 0, 2, 5, 10]
        selected = []
        for pq in preferred:
            nearest = min(cols, key=lambda x: abs(x[0] - pq))
            if nearest not in selected:
                selected.append(nearest)
        selected = selected[:max_lines]
    else:
        selected = []
        for q_target in q_values:
            nearest = min(cols, key=lambda x: abs(x[0] - q_target))
            if nearest not in selected:
                selected.append(nearest)
        selected = selected[:max_lines]

    fig, ax = plt.subplots(figsize=(8, 6))
    for qval, col in selected:
        y = pd.to_numeric(mfdfa_fq[col], errors="coerce")
        mask = (x > 0) & (y > 0)
        ax.plot(np.log10(x[mask]), np.log10(y[mask]), marker="o", linewidth=1.2, label=f"q={qval:g}")

    _apply_common_style(ax, f"MFDFA fluctuation scaling — {run_label}", r"$\log_{10}(s)$", r"$\log_{10}(F_q(s))$")
    ax.legend(ncol=2)
    _savefig(fig, outdir / "plot_mfdfa_scaling.png")


def write_manifest(outdir: Path, run_label: str, summary: dict):
    manifest_path = outdir / "plot_manifest.txt"
    lines = [
        f"Run label: {run_label}",
        "",
        "Generated figures:",
        "  - plot_acf.png",
        "  - plot_binned_series.png",
        "  - plot_delta_alpha_windowed.png",
        "  - plot_delta_h_h0_h5.png",
        "  - plot_multifractal_spectrum.png",
        "  - plot_shuffle_test_summary.png",
        "  - plot_windowed_hq.png",
        "  - plot_mfdfa_scaling.png (optional, only if mfdfa_Fq.csv exists)",
        "",
    ]
    if summary:
        lines.append("Summary snapshot:")
        lines.append(json.dumps(summary, indent=2))
    manifest_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"[saved] {manifest_path}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Create weekly-meeting MFDFA plots from one run directory.")
    p.add_argument("--input-dir", type=Path, default=Path("."), help="Directory containing the CSV outputs.")
    p.add_argument("--output-dir", type=Path, default=None, help="Directory where PNG files will be written.")
    p.add_argument(
        "--hq-q-values", type=float, nargs="*", default=None,
        help="Optional list of q values to show in the windowed h(q) plot."
    )
    p.add_argument(
        "--scaling-q-values", type=float, nargs="*", default=None,
        help="Optional list of q values to show in the MFDFA scaling plot."
    )
    return p


def main():
    args = build_parser().parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = (args.output_dir or input_dir / "plots_weekly_meeting").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    acf = _read_csv(input_dir / "acf.csv")
    binned_series = _read_csv(input_dir / "binned_series.csv")
    delta_alpha = _read_csv(input_dir / "delta_alpha_windowed.csv")
    delta_h = _read_csv(input_dir / "delta_h_h0_h5.csv")
    spectrum = _read_csv(input_dir / "multifractal_spectrum.csv")
    shuffle = _read_csv(input_dir / "shuffle_test_summary.csv")
    windowed_hq = _read_csv(input_dir / "windowed_hq.csv")
    mfdfa_fq = _read_csv(input_dir / "mfdfa_Fq.csv", required=False)
    run_summary = _read_csv(input_dir / "run_summary.csv", required=False)
    summary = _read_json(input_dir / "summary.json")

    run_label = _safe_run_label(summary, run_summary)

    plot_acf(acf, output_dir, run_label)
    plot_binned_series(binned_series, output_dir, run_label)
    plot_delta_alpha(delta_alpha, output_dir, run_label)
    plot_delta_h(delta_h, output_dir, run_label)
    plot_multifractal_spectrum(spectrum, output_dir, run_label)
    plot_shuffle_summary(shuffle, summary, output_dir, run_label)
    plot_windowed_hq(windowed_hq, output_dir, run_label, q_values=args.hq_q_values)
    plot_mfdfa_scaling(mfdfa_fq, output_dir, run_label, q_values=args.scaling_q_values)
    write_manifest(output_dir, run_label, summary)

    print("\nDone. Figures written to:", output_dir)


if __name__ == "__main__":
    main()
