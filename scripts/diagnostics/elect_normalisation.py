"""Elect a temporal normalisation on evidence, not on habit.

`etna-local-msnoise-run.md` section 8 predicted this would matter at Etna:
persistent, strongly localised volcanic tremor at 0.5-5 Hz radiating from the
summit craters violates the diffuse-field assumption, and aggressive temporal
normalisation is the standard defence. So run both and measure.

Usage:  python scripts/diagnostics/elect_normalisation.py A.npz B.npz ...
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from etnalib import config as C, dvv

paths = sys.argv[1:] or [str(C.PROC / "daily_acf.npz"), "/tmp/acf_1bit.npz"]
LAGWIN = (3.0, 15.0)

for p in paths:
    Z = np.load(p, allow_pickle=True)
    lags = Z["lags"]
    m = (lags >= LAGWIN[0]) & (lags <= LAGWIN[1])
    print(f"\n=== {os.path.basename(p)}   [{', '.join(map(str, Z['meta']))}]")
    print(f"{'sta':6s} {'adjacent':>9s} {'vs grand':>9s} {'median cc':>10s} "
          f"{'cc>0.7':>7s} {'dv/v std':>9s} {'saturated':>10s}")
    tot_days = tot_good = 0
    for sta in C.STATIONS:
        k = f"{sta}_acf"
        if k not in Z:
            continue
        A, days = Z[k], list(Z[f"{sta}_days"])
        adj = np.median([np.corrcoef(A[i][m], A[i + 1][m])[0, 1]
                         for i in range(len(A) - 1)])
        gm = A.mean(axis=0)
        vs = np.median([np.corrcoef(a[m], gm[m])[0, 1] for a in A])
        out, _, _ = dvv.dvv_series(A, days, lags, LAGWIN, ref_days=None)
        v = np.array([o[1] for o in out])
        cc = np.array([o[2] for o in out])
        good = cc > 0.7
        sat = int(np.sum(np.abs(v) > 4.9))
        tot_days += len(days)
        tot_good += int(good.sum())
        print(f"{sta:6s} {adj:9.3f} {vs:9.3f} {np.nanmedian(cc):10.3f} "
              f"{int(good.sum()):7d} "
              f"{np.nanstd(v[good]) if good.any() else np.nan:9.3f} {sat:10d}")
    print(f"{'TOTAL':6s} {'':9s} {'':9s} {'':10s} {tot_good:4d}/{tot_days:<4d}")
