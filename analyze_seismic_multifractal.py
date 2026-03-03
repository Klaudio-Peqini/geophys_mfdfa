#!/usr/bin/env python3
"""
analyze_seismic_multifractal.py

A lightweight, research-friendly script (in the style of the
`geomag-seismic-multifractal` package) to:

1) Load an earthquake event catalog from CSV (time, lat, lon, mag)
2) Build an evenly sampled time series (event counts, mean magnitude, or energy proxy)
3) Compute:
   - Autocorrelation function (ACF) (FFT-based, fast for long series)
   - Multifractal Detrended Fluctuation Analysis (MFDFA):
       F_q(s), h(q), tau(q), alpha(q), f(alpha)
   - A few useful summary quantities from the spectrum (width, peak, etc.)
4) Save outputs (CSVs + plots) into an output folder.

Example:
  python analyze_seismic_multifractal.py \
      --catalog /path/to/eq_data.csv \
      --bin 1D \
      --series counts \
      --qmin -5 --qmax 5 --qstep 0.5 \
      --poly 1 \
      --scales 16 24 36 54 81 122 183 274 411 616 924 \
      --maxlag 365 \
      --out results_mag4_counts_1D

Notes:
- For catalogs, the main scientific choice is how to map events -> an evenly sampled
  time series. Common options are:
    * counts  : number of events per bin
    * meanmag : mean magnitude per bin (NaNs -> 0 by default)
    * energy  : sum(10**(1.5*mag)) per bin (energy proxy, Gutenberg-Richter style)
- If your series has strong nonstationarity, consider --log1p and/or --detrend linear.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass
from typing import Iterable, Tuple, Dict, Optional

import numpy as np
import pandas as pd

# Matplotlib is used only for saving figures (no GUI needed).
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------
# I/O and time series building
# -----------------------------

def load_catalog(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "time" not in df.columns:
        raise ValueError("Catalog CSV must contain a 'time' column.")
    df["time"] = pd.to_datetime(df["time"], errors="coerce", utc=True)
    if df["time"].isna().any():
        bad = int(df["time"].isna().sum())
        raise ValueError(f"Failed to parse {bad} timestamps in 'time' column.")
    # Keep only expected columns if present
    for col in ["latitude", "longitude", "mag"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.sort_values("time").reset_index(drop=True)
    return df


def filter_catalog(
    df: pd.DataFrame,
    mag_min: float | None = None,
    mag_max: float | None = None,
    start: str | None = None,
    end: str | None = None,
    lat_min: float | None = None,
    lat_max: float | None = None,
    lon_min: float | None = None,
    lon_max: float | None = None,
) -> pd.DataFrame:
    out = df.copy()
    if mag_min is not None and "mag" in out.columns:
        out = out[out["mag"] >= mag_min]
    if mag_max is not None and "mag" in out.columns:
        out = out[out["mag"] <= mag_max]
    if start is not None:
        out = out[out["time"] >= pd.to_datetime(start, utc=True)]
    if end is not None:
        out = out[out["time"] <= pd.to_datetime(end, utc=True)]
    if lat_min is not None and "latitude" in out.columns:
        out = out[out["latitude"] >= lat_min]
    if lat_max is not None and "latitude" in out.columns:
        out = out[out["latitude"] <= lat_max]
    if lon_min is not None and "longitude" in out.columns:
        out = out[out["longitude"] >= lon_min]
    if lon_max is not None and "longitude" in out.columns:
        out = out[out["longitude"] <= lon_max]
    return out.reset_index(drop=True)


def build_series(
    df: pd.DataFrame,
    bin_rule: str = "1D",
    series_kind: str = "counts",
    fillna: float = 0.0,
) -> pd.Series:
    """
    Convert event catalog -> evenly sampled series.

    Returns a pandas Series indexed by bin start time (UTC).
    """
    if series_kind not in {"counts", "meanmag", "energy"}:
        raise ValueError("series_kind must be one of: counts, meanmag, energy")

    g = df.set_index("time")

    if series_kind == "counts":
        s = g["mag"].resample(bin_rule).size() if "mag" in g.columns else g.resample(bin_rule).size()
        s = s.astype(float)

    elif series_kind == "meanmag":
        if "mag" not in g.columns:
            raise ValueError("meanmag requires a 'mag' column.")
        s = g["mag"].resample(bin_rule).mean()

    else:  # energy proxy
        if "mag" not in g.columns:
            raise ValueError("energy requires a 'mag' column.")
        # proxy ~ 10^(1.5M); absolute scaling doesn't matter for MFDFA
        e = np.power(10.0, 1.5 * g["mag"].astype(float))
        s = e.resample(bin_rule).sum()

    # Fill missing bins
    s = s.asfreq(bin_rule)
    s = s.fillna(fillna)
    s.name = f"{series_kind}_{bin_rule}"
    return s


def preprocess_series(
    s: pd.Series,
    log1p: bool = False,
    standardize: bool = True,
    detrend: str = "none",  # none|linear
) -> np.ndarray:
    x = s.astype(float).to_numpy()

    if log1p:
        x = np.log1p(np.maximum(x, 0.0))

    if detrend == "linear":
        t = np.arange(len(x), dtype=float)
        # least squares fit: x ~ a t + b
        A = np.vstack([t, np.ones_like(t)]).T
        a, b = np.linalg.lstsq(A, x, rcond=None)[0]
        x = x - (a * t + b)

    if standardize:
        mu = float(np.mean(x))
        sig = float(np.std(x))
        if sig > 0:
            x = (x - mu) / sig
        else:
            x = x - mu

    return x


# -----------------------------
# Autocorrelation (FFT-based)
# -----------------------------

def autocorr_fft(x: np.ndarray, max_lag: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Unbiased autocorrelation estimate using FFT.
    Returns (lags, acf) where acf[0] = 1.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 2:
        raise ValueError("Need at least 2 samples for autocorrelation.")

    x = x - np.mean(x)
    # next pow2 for speed
    nfft = 1 << (2 * n - 1).bit_length()
    fx = np.fft.rfft(x, n=nfft)
    acf_full = np.fft.irfft(fx * np.conj(fx), n=nfft)[:n]
    # unbiased normalization by (n-k)
    denom = np.arange(n, 0, -1, dtype=float)
    acf_full = acf_full / denom
    acf_full = acf_full / acf_full[0]

    max_lag = int(min(max_lag, n - 1))
    lags = np.arange(max_lag + 1, dtype=int)
    return lags, acf_full[: max_lag + 1]


# -----------------------------
# MFDFA implementation
# -----------------------------

def _poly_detrend(y: np.ndarray, order: int) -> np.ndarray:
    """Return residuals after polynomial fit of given order."""
    x = np.arange(len(y), dtype=float)
    coeffs = np.polyfit(x, y, deg=order)
    trend = np.polyval(coeffs, x)
    return y - trend


def mfdfa(
    x: np.ndarray,
    scales: Iterable[int],
    qs: np.ndarray,
    poly_order: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute MFDFA fluctuation functions F_q(s) for each scale s and moment q.

    Returns:
      scales_arr: shape (S,)
      Fq: shape (len(qs), S) with Fq[q_i, s_j]
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 16:
        raise ValueError("Series too short for MFDFA (need at least ~16 points).")

    # profile (integrated signal)
    y = np.cumsum(x - np.mean(x))

    scales_arr = np.array(sorted(set(int(s) for s in scales if s >= 4)), dtype=int)
    if scales_arr.size == 0:
        raise ValueError("No valid scales provided (need integers >= 4).")

    qs = np.asarray(qs, dtype=float)
    Fq = np.full((qs.size, scales_arr.size), np.nan, dtype=float)

    for j, s in enumerate(scales_arr):
        ns = n // s
        if ns < 2:
            continue

        # segments: forward and backward for better statistics
        seg_vars = []

        # forward
        for v in range(ns):
            seg = y[v * s : (v + 1) * s]
            res = _poly_detrend(seg, poly_order)
            seg_vars.append(np.mean(res**2))

        # backward
        for v in range(ns):
            start = n - (v + 1) * s
            seg = y[start : start + s]
            res = _poly_detrend(seg, poly_order)
            seg_vars.append(np.mean(res**2))

        seg_vars = np.array(seg_vars, dtype=float)
        # avoid zeros
        seg_vars = np.maximum(seg_vars, 1e-30)

        for i, q in enumerate(qs):
            if abs(q) < 1e-12:
                # q=0 -> geometric mean
                Fq[i, j] = np.exp(0.5 * np.mean(np.log(seg_vars)))
            else:
                Fq[i, j] = (np.mean(seg_vars ** (q / 2.0))) ** (1.0 / q)

    return scales_arr, Fq


def estimate_hq(
    scales: np.ndarray,
    Fq: np.ndarray,
    qs: np.ndarray,
    fit_range: Optional[Tuple[int, int]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit log10(F_q(s)) ~ h(q) log10(s) + c for each q.

    fit_range:
      - None: use all valid scales
      - (smin, smax): only scales within [smin, smax]
    Returns:
      hq: shape (len(qs),)
      r2: shape (len(qs),) coefficient of determination for the fit
    """
    scales = np.asarray(scales, dtype=float)
    qs = np.asarray(qs, dtype=float)

    if fit_range is not None:
        smin, smax = fit_range
        mask = (scales >= smin) & (scales <= smax)
    else:
        mask = np.ones_like(scales, dtype=bool)

    log_s = np.log10(scales[mask])
    hq = np.full(qs.shape, np.nan, dtype=float)
    r2 = np.full(qs.shape, np.nan, dtype=float)

    for i in range(qs.size):
        y = Fq[i, mask]
        ok = np.isfinite(y) & (y > 0)
        if ok.sum() < 3:
            continue
        log_F = np.log10(y[ok])
        x = log_s[ok]

        A = np.vstack([x, np.ones_like(x)]).T
        m, c = np.linalg.lstsq(A, log_F, rcond=None)[0]
        hq[i] = m

        # r^2
        pred = m * x + c
        ss_res = float(np.sum((log_F - pred) ** 2))
        ss_tot = float(np.sum((log_F - np.mean(log_F)) ** 2))
        r2[i] = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return hq, r2


def multifractal_spectrum(qs: np.ndarray, hq: np.ndarray) -> Dict[str, np.ndarray]:
    """
    Compute tau(q), alpha(q), f(alpha) from h(q).
    Uses numerical derivative of tau(q) wrt q.
    """
    qs = np.asarray(qs, dtype=float)
    hq = np.asarray(hq, dtype=float)

    tau = qs * hq - 1.0

    # numerical derivative d tau / d q
    # handle non-uniform q spacing safely:
    alpha = np.gradient(tau, qs)
    falpha = qs * alpha - tau

    return {"q": qs, "h": hq, "tau": tau, "alpha": alpha, "falpha": falpha}


def spectrum_summaries(spec: Dict[str, np.ndarray]) -> Dict[str, float]:
    """
    A few convenient scalars:
      - H (≈ h(2))
      - width of spectrum: alpha_max - alpha_min (finite points)
      - alpha at maximum f(alpha)
      - f(alpha) max
      - h(0), h(2)
    """
    q = spec["q"]
    h = spec["h"]
    alpha = spec["alpha"]
    falpha = spec["falpha"]

    def _pick(q0: float) -> float:
        idx = int(np.argmin(np.abs(q - q0)))
        return float(h[idx])

    ok = np.isfinite(alpha) & np.isfinite(falpha)
    alpha_ok = alpha[ok]
    f_ok = falpha[ok]

    if alpha_ok.size:
        width = float(np.max(alpha_ok) - np.min(alpha_ok))
        imax = int(np.argmax(f_ok))
        alpha_peak = float(alpha_ok[imax])
        f_peak = float(f_ok[imax])
    else:
        width = float("nan")
        alpha_peak = float("nan")
        f_peak = float("nan")

    return {
        "H_h2": _pick(2.0),
        "h0": _pick(0.0),
        "h2": _pick(2.0),
        "spectrum_width_alpha": width,
        "alpha_at_fmax": alpha_peak,
        "falpha_max": f_peak,
    }


# -----------------------------
# Plotting helpers
# -----------------------------

def save_plot_series(s: pd.Series, outpng: str) -> None:
    plt.figure()
    plt.plot(s.index, s.values)
    plt.xlabel("time (UTC)")
    plt.ylabel(s.name or "series")
    plt.title("Binned earthquake series")
    plt.tight_layout()
    plt.savefig(outpng, dpi=180)
    plt.close()


def save_plot_acf(lags: np.ndarray, acf: np.ndarray, outpng: str, bin_rule: str) -> None:
    plt.figure()
    plt.plot(lags, acf)
    plt.xlabel(f"lag (bins of {bin_rule})")
    plt.ylabel("ACF")
    plt.title("Autocorrelation function")
    plt.tight_layout()
    plt.savefig(outpng, dpi=180)
    plt.close()


def save_plot_fqs(scales: np.ndarray, qs: np.ndarray, Fq: np.ndarray, outpng: str) -> None:
    plt.figure()
    for i, q in enumerate(qs):
        y = Fq[i]
        ok = np.isfinite(y) & (y > 0)
        if ok.sum() < 3:
            continue
        plt.plot(np.log10(scales[ok]), np.log10(y[ok]), label=f"q={q:g}")
    plt.xlabel("log10(scale)")
    plt.ylabel("log10(Fq)")
    plt.title("MFDFA fluctuation functions")
    plt.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    plt.savefig(outpng, dpi=180)
    plt.close()


def save_plot_spectrum(alpha: np.ndarray, falpha: np.ndarray, outpng: str) -> None:
    ok = np.isfinite(alpha) & np.isfinite(falpha)
    plt.figure()
    plt.plot(alpha[ok], falpha[ok], marker="o", linestyle="-")
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$f(\alpha)$")
    plt.title("Multifractal spectrum")
    plt.tight_layout()
    plt.savefig(outpng, dpi=180)
    plt.close()


# -----------------------------
# Main
# -----------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="MFDFA + ACF analysis for earthquake catalog CSVs.")
    p.add_argument("--catalog", required=True, help="Path to earthquake catalog CSV (must contain 'time').")
    p.add_argument("--out", required=True, help="Output folder.")
    p.add_argument("--bin", default="1D", help="Pandas resample rule (e.g. 1H, 6H, 1D, 7D). Default: 1D")
    p.add_argument("--series", default="counts", choices=["counts", "meanmag", "energy"],
                   help="Which binned series to analyze. Default: counts")
    p.add_argument("--fillna", type=float, default=0.0, help="Fill value for empty bins. Default: 0.0")

    # filters
    p.add_argument("--magmin", type=float, default=None, help="Minimum magnitude filter.")
    p.add_argument("--magmax", type=float, default=None, help="Maximum magnitude filter.")
    p.add_argument("--start", type=str, default=None, help="Start time (e.g. 1960-01-01).")
    p.add_argument("--end", type=str, default=None, help="End time (e.g. 2020-12-31).")
    p.add_argument("--latmin", type=float, default=None)
    p.add_argument("--latmax", type=float, default=None)
    p.add_argument("--lonmin", type=float, default=None)
    p.add_argument("--lonmax", type=float, default=None)

    # preprocessing
    p.add_argument("--log1p", action="store_true", help="Apply log1p to series (recommended for counts/energy).")
    p.add_argument("--nostd", action="store_true", help="Do NOT standardize (z-score) the series.")
    p.add_argument("--detrend", choices=["none", "linear"], default="none", help="Optional detrending before analysis.")

    # ACF
    p.add_argument("--maxlag", type=int, default=365, help="Max lag in bins for ACF. Default 365.")

    # MFDFA params
    p.add_argument("--poly", type=int, default=1, help="Polynomial order for DFA detrending (m). Default 1.")
    p.add_argument("--qmin", type=float, default=-5.0)
    p.add_argument("--qmax", type=float, default=5.0)
    p.add_argument("--qstep", type=float, default=0.5)
    p.add_argument("--scales", nargs="*", type=int, default=None,
                   help="List of integer scales (window sizes) in bins. If omitted, use log-spaced defaults.")
    p.add_argument("--fit_smin", type=int, default=None, help="Min scale used in h(q) fit.")
    p.add_argument("--fit_smax", type=int, default=None, help="Max scale used in h(q) fit.")

    # Optional: a precomputed correlation-vs-lag CSV (like your *_selflag_singlepair.csv)
    p.add_argument("--precomputed_corr", type=str, default=None,
                   help="Optional path to CSV with columns like lag_hours/lag_days/corr_sign to plot alongside ACF.")

    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    # 1) load + filter catalog
    cat = load_catalog(args.catalog)
    cat_f = filter_catalog(
        cat,
        mag_min=args.magmin, mag_max=args.magmax,
        start=args.start, end=args.end,
        lat_min=args.latmin, lat_max=args.latmax,
        lon_min=args.lonmin, lon_max=args.lonmax,
    )

    # 2) build series
    s = build_series(cat_f, bin_rule=args.bin, series_kind=args.series, fillna=args.fillna)
    s.to_csv(os.path.join(args.out, "binned_series.csv"), index=True)
    save_plot_series(s, os.path.join(args.out, "binned_series.png"))

    # 3) preprocess -> numpy vector
    x = preprocess_series(
        s,
        log1p=args.log1p,
        standardize=(not args.nostd),
        detrend=args.detrend,
    )

    # 4) autocorrelation
    lags, acf = autocorr_fft(x, max_lag=args.maxlag)
    pd.DataFrame({"lag_bins": lags, "acf": acf}).to_csv(os.path.join(args.out, "acf.csv"), index=False)
    save_plot_acf(lags, acf, os.path.join(args.out, "acf.png"), bin_rule=args.bin)

    # optional: plot precomputed correlation file
    if args.precomputed_corr is not None:
        try:
            corr_df = pd.read_csv(args.precomputed_corr)
            # heuristics for columns
            lag_col = "lag_days" if "lag_days" in corr_df.columns else (
                      "lag_hours" if "lag_hours" in corr_df.columns else None)
            corr_col = "corr_sign" if "corr_sign" in corr_df.columns else (
                       "corr" if "corr" in corr_df.columns else None)
            if lag_col and corr_col:
                plt.figure()
                plt.plot(corr_df[lag_col].values, corr_df[corr_col].values)
                plt.xlabel(lag_col)
                plt.ylabel(corr_col)
                plt.title("Precomputed correlation vs lag")
                plt.tight_layout()
                plt.savefig(os.path.join(args.out, "precomputed_corr.png"), dpi=180)
                plt.close()
        except Exception as e:
            print(f"[WARN] Could not read/plot precomputed_corr: {e}")

    # 5) MFDFA
    qs = np.arange(args.qmin, args.qmax + 1e-12, args.qstep, dtype=float)

    if args.scales is None or len(args.scales) == 0:
        # log-spaced defaults (in bins), trimmed to N/4
        n = len(x)
        smin = 16
        smax = max(smin + 1, n // 4)
        # about ~14 scales
        scales = np.unique(np.round(np.logspace(np.log10(smin), np.log10(smax), 14)).astype(int))
        scales = [int(v) for v in scales if v >= 4]
    else:
        scales = args.scales

    scales_arr, Fq = mfdfa(x, scales=scales, qs=qs, poly_order=args.poly)

    # save raw Fq
    # columns: scale, then one column per q
    out_Fq = pd.DataFrame({"scale": scales_arr})
    for i, q in enumerate(qs):
        out_Fq[f"Fq_q{q:g}"] = Fq[i, :]
    out_Fq.to_csv(os.path.join(args.out, "mfdfa_Fq.csv"), index=False)
    save_plot_fqs(scales_arr, qs, Fq, os.path.join(args.out, "mfdfa_Fq.png"))

    fit_range = None
    if args.fit_smin is not None or args.fit_smax is not None:
        smin = args.fit_smin if args.fit_smin is not None else int(np.min(scales_arr))
        smax = args.fit_smax if args.fit_smax is not None else int(np.max(scales_arr))
        fit_range = (smin, smax)

    hq, r2 = estimate_hq(scales_arr, Fq, qs, fit_range=fit_range)
    spec = multifractal_spectrum(qs, hq)
    summ = spectrum_summaries(spec)

    # save spectrum table
    spec_df = pd.DataFrame(spec)
    spec_df["r2_fit"] = r2
    spec_df.to_csv(os.path.join(args.out, "multifractal_spectrum.csv"), index=False)
    save_plot_spectrum(spec["alpha"], spec["falpha"], os.path.join(args.out, "spectrum.png"))

    # save summary json
    meta = {
        "catalog": os.path.abspath(args.catalog),
        "out": os.path.abspath(args.out),
        "bin": args.bin,
        "series": args.series,
        "filters": {
            "magmin": args.magmin, "magmax": args.magmax,
            "start": args.start, "end": args.end,
            "latmin": args.latmin, "latmax": args.latmax,
            "lonmin": args.lonmin, "lonmax": args.lonmax,
        },
        "preprocess": {"log1p": args.log1p, "standardize": (not args.nostd), "detrend": args.detrend},
        "acf": {"maxlag_bins": int(args.maxlag)},
        "mfdfa": {
            "poly_order": int(args.poly),
            "qmin": float(args.qmin), "qmax": float(args.qmax), "qstep": float(args.qstep),
            "scales": [int(v) for v in scales_arr.tolist()],
            "fit_range": fit_range,
        },
        "spectrum_summary": summ,
    }
    with open(os.path.join(args.out, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    # quick console summary
    print("\n=== DONE ===")
    print("Output folder:", os.path.abspath(args.out))
    print("Series:", s.name, "| bins:", len(s))
    print("Spectrum summaries:")
    for k, v in summ.items():
        print(f"  - {k}: {v}")


if __name__ == "__main__":
    main()
