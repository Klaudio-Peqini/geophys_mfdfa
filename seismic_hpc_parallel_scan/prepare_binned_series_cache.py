#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from analyze_seismic_multifractal_fast import load_catalog, filter_catalog, build_series

def save_series_npz(series: pd.Series, out_npz: Path) -> None:
    np.savez_compressed(
        out_npz,
        time_ns=np.asarray(series.index.view("int64")),
        value=np.asarray(series.to_numpy(dtype=float)),
        name=np.asarray([str(series.name) if series.name is not None else "series"]),
    )

def main() -> None:
    p = argparse.ArgumentParser(description="Build cached binned series for scan jobs.")
    p.add_argument("--catalog", required=True)
    p.add_argument("--cache-dir", required=True)
    p.add_argument("--bins", nargs="+", default=["1D", "3D", "7D"])
    p.add_argument("--series", nargs="+", default=["energy", "counts"])
    p.add_argument("--magmins", nargs="+", default=["none", "5.0"])
    p.add_argument("--fillna", type=float, default=0.0)
    args = p.parse_args()

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    cat = load_catalog(args.catalog)

    for bin_rule in args.bins:
        for series_kind in args.series:
            for magmin in args.magmins:
                mag_value = None if str(magmin).lower() == "none" else float(magmin)
                cat_f = filter_catalog(cat, mag_min=mag_value)
                s = build_series(cat_f, bin_rule=bin_rule, series_kind=series_kind, fillna=args.fillna)

                out_csv = cache_dir / f"series_{series_kind}_{bin_rule}_m{magmin}.csv"
                out_npz = cache_dir / f"series_{series_kind}_{bin_rule}_m{magmin}.npz"

                pd.DataFrame({
                    "time": s.index.astype(str),
                    "value": s.values,
                }).to_csv(out_csv, index=False)

                save_series_npz(s, out_npz)

                print(f"Wrote {out_csv}")
                print(f"Wrote {out_npz}")

if __name__ == "__main__":
    main()