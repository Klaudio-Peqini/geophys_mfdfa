#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run_dir", required=True)
    p.add_argument("--run_id", required=True)
    p.add_argument("--bin", required=True)
    p.add_argument("--series", required=True)
    p.add_argument("--magmin", required=True)
    p.add_argument("--window_size_bins", required=True, type=int)
    p.add_argument("--window_step_bins", required=True, type=int)
    p.add_argument("--qmin", required=True, type=float)
    p.add_argument("--qmax", required=True, type=float)
    p.add_argument("--qstep", required=True, type=float)
    p.add_argument("--fit_smin", required=True, type=int)
    p.add_argument("--fit_smax", required=True, type=int)
    args = p.parse_args()

    run_dir = Path(args.run_dir)
    summary_path = run_dir / "summary.json"
    windowed_path = run_dir / "windowed_hq.csv"
    delta_h_csv = run_dir / "delta_h_h0_h5.csv"
    delta_h_png = run_dir / "delta_h_h0_h5.png"
    delta_alpha_csv = run_dir / "delta_alpha_windowed.csv"
    delta_alpha_png = run_dir / "delta_alpha_windowed.png"
    run_summary_csv = run_dir / "run_summary.csv"

    h2 = h0 = width = alpha_at_fmax = falpha_max = float("nan")
    if summary_path.exists():
        with open(summary_path, "r", encoding="utf-8") as f:
            summary = json.load(f)
        spec = summary.get("spectrum_summary", {})
        h2 = spec.get("h2", float("nan"))
        h0 = spec.get("h0", float("nan"))
        width = spec.get("spectrum_width_alpha", float("nan"))
        alpha_at_fmax = spec.get("alpha_at_fmax", float("nan"))
        falpha_max = spec.get("falpha_max", float("nan"))

    delta_h_mean = delta_h_std = delta_h_min = delta_h_max = float("nan")
    n_windows = 0

    if windowed_path.exists():
        df = pd.read_csv(windowed_path)
        if "window_mid_time" in df.columns and "h_q0" in df.columns and "h_q5" in df.columns:
            df["window_mid_time"] = pd.to_datetime(df["window_mid_time"])
            df["delta_h_h0_h5"] = df["h_q0"] - df["h_q5"]
            df.to_csv(delta_h_csv, index=False)

            delta_h_mean = float(df["delta_h_h0_h5"].mean())
            delta_h_std = float(df["delta_h_h0_h5"].std())
            delta_h_min = float(df["delta_h_h0_h5"].min())
            delta_h_max = float(df["delta_h_h0_h5"].max())
            n_windows = int(len(df))

            plt.figure(figsize=(11, 6))
            plt.plot(df["window_mid_time"], df["delta_h_h0_h5"], marker="o", linewidth=1.1)
            plt.xlabel("Window mid-time")
            plt.ylabel(r"$\Delta h(t)=h(0,t)-h(5,t)$")
            plt.title("Temporal evolution of multifractality strength")
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(delta_h_png, dpi=180)
            plt.close()

    row = {
        "run_id": args.run_id,
        "bin": args.bin,
        "series": args.series,
        "magmin": args.magmin,
        "window_size_bins": args.window_size_bins,
        "window_step_bins": args.window_step_bins,
        "qmin": args.qmin,
        "qmax": args.qmax,
        "qstep": args.qstep,
        "fit_smin": args.fit_smin,
        "fit_smax": args.fit_smax,
        "h2": h2,
        "h0": h0,
        "spectrum_width_alpha": width,
        "alpha_at_fmax": alpha_at_fmax,
        "falpha_max": falpha_max,
        "delta_h_mean": delta_h_mean,
        "delta_h_std": delta_h_std,
        "delta_h_min": delta_h_min,
        "delta_h_max": delta_h_max,
        "n_windows": n_windows,
        "summary_json": str(summary_path),
        "windowed_hq_csv": str(windowed_path),
        "delta_h_csv": str(delta_h_csv),
        "delta_h_png": str(delta_h_png),
        "delta_alpha_csv": str(delta_alpha_csv),
        "delta_alpha_png": str(delta_alpha_png),
    }

    pd.DataFrame([row]).to_csv(run_summary_csv, index=False)

if __name__ == "__main__":
    main()
