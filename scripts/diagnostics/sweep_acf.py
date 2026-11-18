"""Find a configuration whose daily autocorrelations are actually stable.

Measures ADJACENT-DAY correlation, not correlation with a reference built
from the same days -- a day inside the reference window correlates with it
by construction, which flatters the result and hides a pipeline that is not
converging.
"""
import sys, os, itertools, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from etnalib import io, acorr, config as C

STA = sys.argv[1] if len(sys.argv) > 1 else "EMSG"
DAYS = io.available_days(STA)[:21]

traces = {}
for d in DAYS:
    tr = io.load_day(STA, d, fs_out=C.ACORR_FS, min_coverage=0.5,
                     remove_response=False)
    if tr is not None:
        traces[d] = tr
print(f"{STA}: {len(traces)} days loaded\n")

CONFIGS = [
    dict(band=(1.0, 3.0), lagwin=(3, 15), norm="ram", whiten=True),
    dict(band=(1.0, 3.0), lagwin=(1, 6), norm="ram", whiten=True),
    dict(band=(0.5, 2.0), lagwin=(3, 15), norm="ram", whiten=True),
    dict(band=(0.5, 2.0), lagwin=(5, 20), norm="ram", whiten=True),
    dict(band=(0.2, 1.0), lagwin=(5, 20), norm="ram", whiten=True),
    dict(band=(0.1, 0.5), lagwin=(5, 20), norm="ram", whiten=True),
    dict(band=(1.0, 3.0), lagwin=(3, 15), norm="onebit", whiten=True),
    dict(band=(1.0, 3.0), lagwin=(3, 15), norm="ram", whiten=False),
    dict(band=(0.5, 2.0), lagwin=(3, 15), norm="ram", whiten=False),
    dict(band=(0.5, 2.0), lagwin=(3, 15), norm="onebit", whiten=True),
]

print(f"{'band':12s} {'lag win':10s} {'norm':7s} {'whit':5s} "
      f"{'adj-day cc':>11s} {'cc vs grand mean':>17s} {'t/day':>7s}")
for cfg in CONFIGS:
    t0 = time.time()
    acfs, lags = [], None
    for d, tr in traces.items():
        lg, c, nw = acorr.daily_autocorr(
            tr, cfg["band"], C.ACORR_MAXLAG, window_s=C.ACORR_WINDOW,
            norm=cfg["norm"], whiten=cfg["whiten"], hours=C.PARTIAL_HOURS)
        if c is None or nw < 4:
            continue
        lags = lg
        acfs.append(c)
    if len(acfs) < 5:
        print(f"{str(cfg['band']):12s} -- too few days")
        continue
    A = np.vstack(acfs)
    m = (lags >= cfg["lagwin"][0]) & (lags <= cfg["lagwin"][1])
    adj = [np.corrcoef(A[i][m], A[i + 1][m])[0, 1] for i in range(len(A) - 1)]
    gm = A.mean(axis=0)
    vs = [np.corrcoef(a[m], gm[m])[0, 1] for a in A]
    dt = (time.time() - t0) / len(acfs)
    print(f"{str(cfg['band']):12s} {str(cfg['lagwin']):10s} {cfg['norm']:7s} "
          f"{str(cfg['whiten']):5s} {np.median(adj):11.3f} {np.median(vs):17.3f} "
          f"{dt:7.2f}s")
