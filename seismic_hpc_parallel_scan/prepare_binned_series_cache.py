#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from analyze_seismic_multifractal_fast import load_catalog, filter_catalog, build_series


def main() -> None:
    p = argparse.ArgumentParser(description='Precompute reusable binned series for parallel scans.')
    p.add_argument('--catalog', required=True)
    p.add_argument('--cache-dir', required=True)
    p.add_argument('--bins', nargs='+', default=['1D', '3D', '7D'])
    p.add_argument('--series', nargs='+', default=['energy', 'counts'])
    p.add_argument('--magmins', nargs='+', default=['none', '5.0'])
    p.add_argument('--fillna', type=float, default=0.0)
    args = p.parse_args()

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    cat = load_catalog(args.catalog)
    for bin_rule in args.bins:
        for series_kind in args.series:
            for magmin in args.magmins:
                mag_value = None if str(magmin).lower() == 'none' else float(magmin)
                cat_f = filter_catalog(cat, mag_min=mag_value)
                s = build_series(cat_f, bin_rule=bin_rule, series_kind=series_kind, fillna=args.fillna)
                out = cache_dir / f'series_{series_kind}_{bin_rule}_m{magmin}.csv'
                df = pd.DataFrame({'time': s.index.astype(str), 'value': s.values})
                df.to_csv(out, index=False)
                print(f'Wrote {out}')


if __name__ == '__main__':
    main()
