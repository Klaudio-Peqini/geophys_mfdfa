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

#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, os
from typing import Dict, Iterable, Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def set_plot_style() -> None:
    plt.rcParams.update({
        "figure.figsize": (10, 6),
        "figure.dpi": 160,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.titlesize": 18,
        "axes.labelsize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 9,
        "lines.markersize": 6,
        "savefig.bbox": "tight",
    })

def load_catalog(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "time" not in df.columns:
        raise ValueError("Catalog CSV must contain a 'time' column.")
    #df["time"] = pd.to_datetime(df["time"], errors="coerce", utc=True)
    #df["time"] = pd.to_datetime(
    #    df["time"],
    #    infer_datetime_format=True,
    #    errors="coerce",
    #    utc=True
    #)
    df["time"] = pd.to_datetime(
        df["time"],
        format="%Y-%m-%d %H:%M:%S",
        errors="coerce",
        utc=True
    )
    if df["time"].isna().any():
        raise ValueError(f"Failed to parse {int(df['time'].isna().sum())} timestamps.")
    for col in ["latitude", "longitude", "mag"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.sort_values("time").reset_index(drop=True)

def filter_catalog(
    df: pd.DataFrame,
    mag_min=None, mag_max=None,
    start=None, end=None,
    lat_min=None, lat_max=None,
    lon_min=None, lon_max=None
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

def build_series(df: pd.DataFrame, bin_rule="1D", series_kind="counts", fillna=0.0) -> pd.Series:
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
    else:
        if "mag" not in g.columns:
            raise ValueError("energy requires a 'mag' column.")
        e = np.power(10.0, 1.5 * g["mag"].astype(float))
        s = e.resample(bin_rule).sum()

    s = s.asfreq(bin_rule).fillna(fillna)
    s.name = f"{series_kind}_{bin_rule}"
    return s

def preprocess_series(s: pd.Series, log1p=False, standardize=True, detrend="none") -> np.ndarray:
    x = s.astype(float).to_numpy()

    if log1p:
        x = np.log1p(np.maximum(x, 0.0))

    if detrend == "linear":
        t = np.arange(len(x), dtype=float)
        A = np.vstack([t, np.ones_like(t)]).T
        a, b = np.linalg.lstsq(A, x, rcond=None)[0]
        x = x - (a * t + b)

    if standardize:
        mu, sig = float(np.mean(x)), float(np.std(x))
        x = x - mu if sig <= 0 else (x - mu) / sig

    return x

def autocorr_fft(x: np.ndarray, max_lag: int):
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 2:
        raise ValueError("Need at least 2 samples.")
    x = x - np.mean(x)
    nfft = 1 << (2 * n - 1).bit_length()
    fx = np.fft.rfft(x, n=nfft)
    acf_full = np.fft.irfft(fx * np.conj(fx), n=nfft)[:n]
    acf_full = (acf_full / np.arange(n, 0, -1, dtype=float)) / acf_full[0]
    max_lag = int(min(max_lag, n - 1))
    lags = np.arange(max_lag + 1, dtype=int)
    return lags, acf_full[:max_lag + 1]

def acf_significance_white_noise(n: int) -> float:
    return 1.96 / np.sqrt(max(n, 1))

def acf_shuffle_band(x: np.ndarray, max_lag: int, n_shuffles=200, seed=42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    sh = []
    for _ in range(n_shuffles):
        xs = np.array(x, copy=True)
        rng.shuffle(xs)
        _, a = autocorr_fft(xs, max_lag=max_lag)
        sh.append(np.abs(a))
    return np.quantile(np.vstack(sh), 0.95, axis=0)

def _poly_detrend(y: np.ndarray, order: int) -> np.ndarray:
    x = np.arange(len(y), dtype=float)
    coeffs = np.polyfit(x, y, deg=order)
    return y - np.polyval(coeffs, x)

def mfdfa(x: np.ndarray, scales: Iterable[int], qs: np.ndarray, poly_order=1, var_floor=1e-30):
    x = np.asarray(x, dtype=float)
    if x.size < 16:
        raise ValueError("Series too short for MFDFA.")

    y = np.cumsum(x - np.mean(x))
    scales_arr = np.array(sorted(set(int(s) for s in scales if s >= 4)), dtype=int)
    if scales_arr.size == 0:
        raise ValueError("No valid scales provided.")

    Fq = np.full((len(qs), len(scales_arr)), np.nan, dtype=float)

    for j, s in enumerate(scales_arr):
        ns = len(x) // s
        if ns < 2:
            continue

        seg_vars = []
        for v in range(ns):
            seg = y[v * s:(v + 1) * s]
            seg_vars.append(np.mean(_poly_detrend(seg, poly_order) ** 2))
        for v in range(ns):
            start = len(x) - (v + 1) * s
            seg = y[start:start + s]
            seg_vars.append(np.mean(_poly_detrend(seg, poly_order) ** 2))

        seg_vars = np.maximum(np.array(seg_vars, dtype=float), var_floor)

        for i, q in enumerate(qs):
            if abs(q) < 1e-14:
                Fq[i, j] = np.exp(0.5 * np.mean(np.log(seg_vars)))
            else:
                val = np.mean(seg_vars ** (q / 2.0))
                if np.isfinite(val) and val > 0:
                    Fq[i, j] = val ** (1.0 / q)

    return scales_arr, Fq

def fit_line_for_q(scales: np.ndarray, Fq_row: np.ndarray, fit_range: Optional[Tuple[int, int]]):
    scales = np.asarray(scales, dtype=float)
    mask = np.ones_like(scales, dtype=bool) if fit_range is None else ((scales >= fit_range[0]) & (scales <= fit_range[1]))
    y = Fq_row[mask]
    x = np.log10(scales[mask])

    ok = np.isfinite(y) & (y > 0)
    if ok.sum() < 3:
        return np.nan, np.nan, np.nan, np.zeros_like(mask, dtype=bool)

    logF, x_ok = np.log10(y[ok]), x[ok]
    A = np.vstack([x_ok, np.ones_like(x_ok)]).T
    m, c = np.linalg.lstsq(A, logF, rcond=None)[0]

    pred = m * x_ok + c
    ss_res = float(np.sum((logF - pred) ** 2))
    ss_tot = float(np.sum((logF - np.mean(logF)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    full_mask = np.zeros_like(mask, dtype=bool)
    idx = np.where(mask)[0]
    full_mask[idx[ok]] = True

    return float(m), float(c), float(r2), full_mask

def estimate_hq(scales: np.ndarray, Fq: np.ndarray, qs: np.ndarray, fit_range=None):
    hq = np.full(qs.shape, np.nan)
    intercept = np.full(qs.shape, np.nan)
    r2 = np.full(qs.shape, np.nan)
    for i in range(qs.size):
        hq[i], intercept[i], r2[i], _ = fit_line_for_q(scales, Fq[i], fit_range)
    return hq, intercept, r2

def multifractal_spectrum(qs: np.ndarray, hq: np.ndarray) -> Dict[str, np.ndarray]:
    tau = qs * hq - 1.0
    alpha = np.gradient(tau, qs)
    falpha = qs * alpha - tau
    return {"q": qs, "h": hq, "tau": tau, "alpha": alpha, "falpha": falpha}

def spectrum_summaries(spec: Dict[str, np.ndarray]) -> Dict[str, float]:
    q, h, alpha, falpha = spec["q"], spec["h"], spec["alpha"], spec["falpha"]
    pick = lambda q0: float(h[int(np.argmin(np.abs(q - q0)))])
    ok = np.isfinite(alpha) & np.isfinite(falpha)

    if ok.any():
        alpha_ok, f_ok = alpha[ok], falpha[ok]
        width = float(np.max(alpha_ok) - np.min(alpha_ok))
        imax = int(np.argmax(f_ok))
        alpha_peak, f_peak = float(alpha_ok[imax]), float(f_ok[imax])
    else:
        width = alpha_peak = f_peak = float("nan")

    return {
        "H_h2": pick(2.0),
        "h0": pick(0.0),
        "h2": pick(2.0),
        "spectrum_width_alpha": width,
        "alpha_at_fmax": alpha_peak,
        "falpha_max": f_peak
    }

def compute_windowed_hq(x, time_index, qs, scales, poly_order, fit_range, window_size_bins, window_step_bins, var_floor):
    rows = []
    n = len(x)

    for start in range(0, n - window_size_bins + 1, window_step_bins):
        end = start + window_size_bins
        xw = x[start:end]
        tw = time_index[start:end]

        scales_arr, Fq = mfdfa(xw, scales=scales, qs=qs, poly_order=poly_order, var_floor=var_floor)
        hq, _, r2 = estimate_hq(scales_arr, Fq, qs, fit_range=fit_range)

        row = {
            "window_start_idx": start,
            "window_end_idx": end - 1,
            "window_start_time": str(tw[0]),
            "window_end_time": str(tw[-1]),
            "window_mid_time": str(tw[len(tw) // 2]),
        }
        for q, h, r in zip(qs, hq, r2):
            row[f"h_q{q:g}"] = h
            row[f"r2_q{q:g}"] = r
        rows.append(row)

    return pd.DataFrame(rows)

def shuffled_multifractal_test(x, qs, scales, poly_order, fit_range, n_shuffles, seed, var_floor):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_shuffles):
        xs = np.array(x, copy=True)
        rng.shuffle(xs)
        scales_arr, Fq = mfdfa(xs, scales=scales, qs=qs, poly_order=poly_order, var_floor=var_floor)
        hq, _, _ = estimate_hq(scales_arr, Fq, qs, fit_range=fit_range)
        row = {"shuffle_id": i}
        row.update(spectrum_summaries(multifractal_spectrum(qs, hq)))
        rows.append(row)
    return pd.DataFrame(rows)

def phase_randomized_surrogate(x: np.ndarray, seed: int | None = None) -> np.ndarray:
    """
    Create a phase-randomized surrogate preserving the linear power spectrum
    approximately, while destroying nonlinear temporal structure.
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    n = len(x)

    # Fourier transform
    Xf = np.fft.rfft(x)
    amp = np.abs(Xf)
    phase = np.angle(Xf)

    # Random phases, keeping DC and Nyquist (if present) unchanged
    random_phases = rng.uniform(0, 2 * np.pi, size=len(Xf))
    random_phases[0] = phase[0]
    if n % 2 == 0 and len(Xf) > 1:
        random_phases[-1] = phase[-1]

    Xf_new = amp * np.exp(1j * random_phases)
    xs = np.fft.irfft(Xf_new, n=n)

    return xs


def phase_randomized_multifractal_test(
    x: np.ndarray,
    qs: np.ndarray,
    scales: np.ndarray,
    poly_order: int,
    fit_range: Optional[Tuple[int, int]],
    n_surrogates: int,
    seed: int,
    var_floor: float,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []

    for i in range(n_surrogates):
        xs = phase_randomized_surrogate(x, seed=int(rng.integers(0, 10**9)))
        scales_arr, Fq = mfdfa(
            xs, scales=scales, qs=qs,
            poly_order=poly_order, var_floor=var_floor
        )
        hq, _, r2 = estimate_hq(scales_arr, Fq, qs, fit_range=fit_range)
        spec = multifractal_spectrum(qs, hq)
        summ = spectrum_summaries(spec)

        row = {"phase_id": i}
        row.update(summ)
        rows.append(row)

    return pd.DataFrame(rows)

def save_plot_series(s, outpng, title):
    plt.figure()
    plt.plot(s.index, s.values, lw=1.0)
    plt.xlabel("Time (UTC)")
    plt.ylabel(s.name or "Series value")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()

def save_plot_acf(lags, acf, outpng, bin_rule, analytic_band=None, shuffle_band=None):
    plt.figure()
    markerline, stemlines, baseline = plt.stem(lags, acf)
    plt.setp(markerline, markersize=4)
    plt.setp(stemlines, linewidth=1.0)
    plt.setp(baseline, linewidth=0.8)

    if analytic_band is not None:
        plt.axhline(analytic_band, ls="--", lw=1.2, label=f"95% white-noise band (+/-{analytic_band:.3f})")
        plt.axhline(-analytic_band, ls="--", lw=1.2)
    if shuffle_band is not None:
        plt.plot(lags, shuffle_band, lw=1.2, label="95% shuffled envelope")
        plt.plot(lags, -shuffle_band, lw=1.2)

    plt.xlabel(f"Lag (bins of {bin_rule})")
    plt.ylabel("Autocorrelation")
    plt.title("Autocorrelation function with significance thresholds")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()

def save_plot_fqs_with_fits(scales, qs, Fq, outpng, fit_range, title="MFDFA fluctuation functions"):
    plt.figure(figsize=(12, 8))
    markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', '<', '>', 'h', '8']
    cmap = plt.cm.turbo
    n_q = len(qs)

    for i, q in enumerate(qs):
        y = Fq[i]
        ok = np.isfinite(y) & (y > 0)
        if ok.sum() < 3:
            continue

        color = cmap(i / max(n_q - 1, 1))
        marker = markers[i % len(markers)]
        x_plot, y_plot = np.log10(scales[ok]), np.log10(y[ok])

        plt.scatter(
            x_plot, y_plot,
            marker=marker, s=28,
            facecolors="none", edgecolors=color,
            linewidths=1.0,
            label=f"q={q:g}"
        )

        m, c, _, fit_mask = fit_line_for_q(scales, y, fit_range)
        if np.isfinite(m):
            x_fit = np.sort(np.log10(scales[fit_mask]))
            plt.plot(x_fit, m * x_fit + c, color=color, lw=1.0, alpha=0.8)

    if fit_range is not None:
        plt.axvline(np.log10(fit_range[0]), color="gray", ls="--", lw=1.0)
        plt.axvline(np.log10(fit_range[1]), color="gray", ls="--", lw=1.0)

    plt.xlabel(r"$\log_{10}(s)$")
    plt.ylabel(r"$\log_{10}(F_q(s))$")
    plt.title(title)
    if n_q <= 25:
        plt.legend(ncol=2)
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()

def save_plot_hq(qs, hq, r2, outpng, title):
    ok = np.isfinite(hq)
    plt.figure()
    plt.scatter(qs[ok], hq[ok], s=28)
    plt.plot(qs[ok], hq[ok], lw=1.0)
    plt.xlabel(r"$q$")
    plt.ylabel(r"$h(q)$")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()

def save_plot_spectrum(alpha, falpha, outpng, title):
    ok = np.isfinite(alpha) & np.isfinite(falpha)
    plt.figure()
    plt.scatter(alpha[ok], falpha[ok], s=28)
    plt.plot(alpha[ok], falpha[ok], lw=1.0)
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$f(\alpha)$")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()

def save_plot_windowed_hq(dfw, q_values, outpng):
    if dfw.empty:
        return
    plt.figure(figsize=(12, 7))
    times = pd.to_datetime(dfw["window_mid_time"])
    for q in q_values:
        col = f"h_q{q:g}"
        if col in dfw.columns:
            plt.plot(times, dfw[col], marker="o", lw=1.0, label=fr"$h({q:g})$")
    plt.xlabel("Window mid-time")
    plt.ylabel(r"$h(q)$")
    plt.title(r"Temporal evolution of selected $h(q)$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()

def save_plot_shuffle_test(real_width, shuffle_df, outpng):
    if shuffle_df.empty or "spectrum_width_alpha" not in shuffle_df.columns:
        return
    vals = shuffle_df["spectrum_width_alpha"].dropna().to_numpy()
    if vals.size == 0:
        return
    plt.figure()
    plt.hist(vals, bins=min(20, max(5, vals.size // 2)), alpha=0.8)
    plt.axvline(real_width, color="red", lw=2.0, label=f"Real width = {real_width:.4g}")
    plt.xlabel(r"Shuffled $\Delta \alpha$")
    plt.ylabel("Count")
    plt.title("Shuffle test for spectrum width")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()

def save_plot_multifractality_strength(dfw: pd.DataFrame, q_left: float, q_right: float, outpng: str) -> None:
    col_left = f"h_q{q_left:g}"
    col_right = f"h_q{q_right:g}"

    if dfw.empty or col_left not in dfw.columns or col_right not in dfw.columns:
        return

    times = pd.to_datetime(dfw["window_mid_time"])
    delta_h = dfw[col_left] - dfw[col_right]

    out_df = pd.DataFrame({
        "window_mid_time": times,
        "delta_h": delta_h
    })
    out_df.to_csv(outpng.replace(".png", ".csv"), index=False)

    plt.figure(figsize=(12, 6))
    plt.plot(times, delta_h, marker="o", lw=1.2)
    plt.xlabel("Window mid-time")
    plt.ylabel(r"$\Delta h$")
    plt.title(fr"Temporal evolution of multifractality strength: $h({q_left:g})-h({q_right:g})$")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()

def save_plot_delta_alpha_windowed(
    x: np.ndarray,
    time_index: pd.Index,
    qs: np.ndarray,
    scales: np.ndarray,
    poly_order: int,
    fit_range: Optional[Tuple[int, int]],
    window_size_bins: int,
    window_step_bins: int,
    var_floor: float,
    outpng: str,
) -> None:
    rows = []
    n = len(x)

    for start in range(0, n - window_size_bins + 1, window_step_bins):
        end = start + window_size_bins
        xw = x[start:end]
        tw = time_index[start:end]

        scales_arr, Fq = mfdfa(xw, scales=scales, qs=qs, poly_order=poly_order, var_floor=var_floor)
        hq, _, _ = estimate_hq(scales_arr, Fq, qs, fit_range=fit_range)
        spec = multifractal_spectrum(qs, hq)

        alpha = spec["alpha"]
        ok = np.isfinite(alpha)
        if np.any(ok):
            delta_alpha = float(np.max(alpha[ok]) - np.min(alpha[ok]))
        else:
            delta_alpha = np.nan

        rows.append({
            "window_mid_time": str(tw[len(tw)//2]),
            "delta_alpha": delta_alpha
        })

    dfa = pd.DataFrame(rows)
    dfa.to_csv(outpng.replace(".png", ".csv"), index=False)

    times = pd.to_datetime(dfa["window_mid_time"])
    plt.figure(figsize=(12, 6))
    plt.plot(times, dfa["delta_alpha"], marker="s", lw=1.2)
    plt.xlabel("Window mid-time")
    plt.ylabel(r"$\Delta \alpha$")
    plt.title(r"Temporal evolution of spectrum width $\Delta \alpha$")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outpng)
    plt.close()

def parse_args():
    p = argparse.ArgumentParser(description="Improved MFDFA + ACF analysis for earthquake catalog CSVs.")
    p.add_argument("--catalog", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--bin", default="1D")
    p.add_argument("--series", default="counts", choices=["counts", "meanmag", "energy"])
    p.add_argument("--fillna", type=float, default=0.0)

    p.add_argument("--magmin", type=float, default=None)
    p.add_argument("--magmax", type=float, default=None)
    p.add_argument("--start", type=str, default=None)
    p.add_argument("--end", type=str, default=None)
    p.add_argument("--latmin", type=float, default=None)
    p.add_argument("--latmax", type=float, default=None)
    p.add_argument("--lonmin", type=float, default=None)
    p.add_argument("--lonmax", type=float, default=None)

    p.add_argument("--log1p", action="store_true")
    p.add_argument("--nostd", action="store_true")
    p.add_argument("--detrend", choices=["none", "linear"], default="none")

    p.add_argument("--maxlag", type=int, default=365)
    p.add_argument("--acf_shuffle_band", action="store_true")
    p.add_argument("--acf_shuffle_count", type=int, default=200)

    p.add_argument("--poly", type=int, default=1)
    p.add_argument("--qmin", type=float, default=-30.0)
    p.add_argument("--qmax", type=float, default=30.0)
    p.add_argument("--qstep", type=float, default=1.0)
    p.add_argument("--scales", nargs="*", type=int, default=None)
    p.add_argument("--fit_smin", type=int, default=None)
    p.add_argument("--fit_smax", type=int, default=None)
    p.add_argument("--var_floor", type=float, default=1e-30)
    p.add_argument("--precomputed_corr", type=str, default=None)

    p.add_argument("--window_size_bins", type=int, default=None)
    p.add_argument("--window_step_bins", type=int, default=None)
    p.add_argument("--window_q_values", nargs="*", type=float, default=[0.0, 2.0, 5.0])

    p.add_argument("--part_start", type=str, default=None)
    p.add_argument("--part_end", type=str, default=None)

    p.add_argument("--n_shuffles", type=int, default=0)
    p.add_argument("--shuffle_seed", type=int, default=42)

    p.add_argument("--n_phase_surrogates", type=int, default=0, help="Number of phase-randomized surrogates for multifractal comparison.")
    return p.parse_args()

def main():
    set_plot_style()
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    cat = load_catalog(args.catalog)
    cat_f = filter_catalog(
        cat,
        mag_min=args.magmin, mag_max=args.magmax,
        start=args.start, end=args.end,
        lat_min=args.latmin, lat_max=args.latmax,
        lon_min=args.lonmin, lon_max=args.lonmax
    )

    s = build_series(cat_f, bin_rule=args.bin, series_kind=args.series, fillna=args.fillna)

    if args.part_start is not None:
        s = s[s.index >= pd.to_datetime(args.part_start, utc=True)]
    if args.part_end is not None:
        s = s[s.index <= pd.to_datetime(args.part_end, utc=True)]
    if s.empty:
        raise ValueError("No data remain after filtering.")

    s.to_csv(os.path.join(args.out, "binned_series.csv"), index=True)
    save_plot_series(s, os.path.join(args.out, "binned_series.png"), title=f"Binned seismic series ({args.series}, bin={args.bin})")

    x = preprocess_series(s, log1p=args.log1p, standardize=(not args.nostd), detrend=args.detrend)

    lags, acf = autocorr_fft(x, max_lag=args.maxlag)
    analytic_band = acf_significance_white_noise(len(x))
    shuffle_band = acf_shuffle_band(x, args.maxlag, args.acf_shuffle_count, args.shuffle_seed) if args.acf_shuffle_band else None

    acf_df = pd.DataFrame({
        "lag_bins": lags,
        "acf": acf,
        "white_noise_95_band": analytic_band,
        "acf_physically_relevant_analytic": np.abs(acf) > analytic_band
    })
    if shuffle_band is not None:
        acf_df["shuffle_95_band"] = shuffle_band
        acf_df["acf_physically_relevant_shuffle"] = np.abs(acf) > shuffle_band

    acf_df.to_csv(os.path.join(args.out, "acf.csv"), index=False)
    save_plot_acf(lags, acf, os.path.join(args.out, "acf.png"), args.bin, analytic_band, shuffle_band)

    if args.precomputed_corr is not None:
        try:
            corr_df = pd.read_csv(args.precomputed_corr)
            lag_col = "lag_days" if "lag_days" in corr_df.columns else ("lag_hours" if "lag_hours" in corr_df.columns else None)
            corr_col = "corr_sign" if "corr_sign" in corr_df.columns else ("corr" if "corr" in corr_df.columns else None)
            if lag_col and corr_col:
                plt.figure()
                plt.scatter(corr_df[lag_col].values, corr_df[corr_col].values, s=18)
                plt.xlabel(lag_col)
                plt.ylabel(corr_col)
                plt.title("Precomputed correlation vs lag")
                plt.tight_layout()
                plt.savefig(os.path.join(args.out, "precomputed_corr.png"))
                plt.close()
        except Exception as e:
            print(f"[WARN] Could not read/plot precomputed_corr: {e}")

    qs = np.arange(args.qmin, args.qmax + 1e-12, args.qstep, dtype=float)

    if args.scales is None or len(args.scales) == 0:
        n = len(x)
        smin = 16
        smax = max(smin + 1, n // 4)
        scales = [int(v) for v in np.unique(np.round(np.logspace(np.log10(smin), np.log10(smax), 14)).astype(int)) if int(v) >= 4]
    else:
        scales = args.scales

    fit_range = None if (args.fit_smin is None and args.fit_smax is None) else (
        (args.fit_smin if args.fit_smin is not None else int(np.min(scales))),
        (args.fit_smax if args.fit_smax is not None else int(np.max(scales)))
    )

    scales_arr, Fq = mfdfa(x, scales=scales, qs=qs, poly_order=args.poly, var_floor=args.var_floor)

    out_Fq = pd.DataFrame({"scale": scales_arr})
    for i, q in enumerate(qs):
        out_Fq[f"Fq_q{q:g}"] = Fq[i, :]
    out_Fq.to_csv(os.path.join(args.out, "mfdfa_Fq.csv"), index=False)

    save_plot_fqs_with_fits(
        scales_arr, qs, Fq,
        os.path.join(args.out, "mfdfa_Fq.png"),
        fit_range,
        title=f"MFDFA fluctuation functions ({args.series}, q in [{args.qmin:g}, {args.qmax:g}])"
    )

    hq, intercept, r2 = estimate_hq(scales_arr, Fq, qs, fit_range=fit_range)
    spec = multifractal_spectrum(qs, hq)
    summ = spectrum_summaries(spec)

    spec_df = pd.DataFrame(spec)
    spec_df["intercept_fit"] = intercept
    spec_df["r2_fit"] = r2
    spec_df.to_csv(os.path.join(args.out, "multifractal_spectrum.csv"), index=False)

    save_plot_hq(qs, hq, r2, os.path.join(args.out, "hq_vs_q.png"), title=r"Generalized Hurst exponent $h(q)$ vs $q$")
    save_plot_spectrum(spec["alpha"], spec["falpha"], os.path.join(args.out, "spectrum.png"), title="Multifractal spectrum")

    if args.window_size_bins is not None:
        step = args.window_step_bins if args.window_step_bins is not None else max(1, args.window_size_bins // 4)
        windowed_df = compute_windowed_hq(x, s.index, qs, np.asarray(scales_arr, dtype=int), args.poly, fit_range, args.window_size_bins, step, args.var_floor)
        windowed_df.to_csv(os.path.join(args.out, "windowed_hq.csv"), index=False)
        save_plot_windowed_hq(windowed_df, args.window_q_values, os.path.join(args.out, "windowed_hq_selected.png"))
        save_plot_multifractality_strength(windowed_df, q_left=0, q_right=5, outpng=os.path.join(args.out, "multifractality_strength_h0_h5.png"),)
        save_plot_delta_alpha_windowed(
            x=x,
            time_index=s.index,
            qs=qs,
            scales=np.asarray(scales_arr, dtype=int),
            poly_order=args.poly,
            fit_range=fit_range,
            window_size_bins=args.window_size_bins,
            window_step_bins=step,
            var_floor=args.var_floor,
            outpng=os.path.join(args.out, "delta_alpha_windowed.png"),)

    shuffle_summary = {}
    if args.n_shuffles > 0:
        shuffle_df = shuffled_multifractal_test(x, qs, np.asarray(scales_arr, dtype=int), args.poly, fit_range, args.n_shuffles, args.shuffle_seed, args.var_floor)
        shuffle_df.to_csv(os.path.join(args.out, "shuffle_test_summary.csv"), index=False)
        save_plot_shuffle_test(summ["spectrum_width_alpha"], shuffle_df, os.path.join(args.out, "shuffle_test_width.png"))
        vals = shuffle_df["spectrum_width_alpha"].dropna().to_numpy()
        if vals.size:
            shuffle_summary = {
                "mean_shuffled_width": float(np.mean(vals)),
                "std_shuffled_width": float(np.std(vals)),
                "p_value_width_ge_real": float(np.mean(vals >= summ["spectrum_width_alpha"]))
            }
    
        phase_df = pd.DataFrame()
    phase_summary = {}

    if args.n_phase_surrogates > 0:
        phase_df = phase_randomized_multifractal_test(
            x=x,
            qs=qs,
            scales=np.asarray(scales_arr, dtype=int),
            poly_order=args.poly,
            fit_range=fit_range,
            n_surrogates=args.n_phase_surrogates,
            seed=args.shuffle_seed + 1000,
            var_floor=args.var_floor,
        )

        phase_df.to_csv(os.path.join(args.out, "phase_surrogate_test_summary.csv"), index=False)

        if "spectrum_width_alpha" in phase_df.columns and not phase_df["spectrum_width_alpha"].dropna().empty:
            phase_widths = phase_df["spectrum_width_alpha"].dropna().to_numpy()
            p_ge_real_phase = float(np.mean(phase_widths >= summ["spectrum_width_alpha"]))

            phase_summary = {
                "mean_phase_width": float(np.mean(phase_widths)),
                "std_phase_width": float(np.std(phase_widths)),
                "p_value_phase_width_ge_real": p_ge_real_phase,
            }

    meta = {
        "catalog": os.path.abspath(args.catalog),
        "out": os.path.abspath(args.out),
        "bin": args.bin,
        "series": args.series,
        "filters": {
            "magmin": args.magmin, "magmax": args.magmax,
            "start": args.start, "end": args.end,
            "latmin": args.latmin, "latmax": args.latmax,
            "lonmin": args.lonmin, "lonmax": args.lonmax
        },
        "subperiod": {
            "part_start": args.part_start,
            "part_end": args.part_end
        },
        "preprocess": {
            "log1p": args.log1p,
            "standardize": (not args.nostd),
            "detrend": args.detrend
        },
        "acf": {
            "maxlag_bins": int(args.maxlag),
            "white_noise_95_band": analytic_band,
            "used_shuffle_band": bool(args.acf_shuffle_band),
            "acf_shuffle_count": int(args.acf_shuffle_count) if args.acf_shuffle_band else 0
        },
        "mfdfa": {
            "poly_order": int(args.poly),
            "qmin": float(args.qmin),
            "qmax": float(args.qmax),
            "qstep": float(args.qstep),
            "var_floor": float(args.var_floor),
            "scales": [int(v) for v in scales_arr.tolist()],
            "fit_range": fit_range
        },
        "windowed_hq": {
            "enabled": args.window_size_bins is not None,
            "window_size_bins": args.window_size_bins,
            "window_step_bins": args.window_step_bins,
            "window_q_values": args.window_q_values
        },
        "shuffle_test": {
            "n_shuffles": int(args.n_shuffles),
            "seed": int(args.shuffle_seed),
            **shuffle_summary
        },
        "phase_surrogate_test": {
            "n_phase_surrogates": int(args.n_phase_surrogates),
            **phase_summary,
        },
        "spectrum_summary": summ
    }

    with open(os.path.join(args.out, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    print("DONE:", os.path.abspath(args.out))

if __name__ == "__main__":
    main()
