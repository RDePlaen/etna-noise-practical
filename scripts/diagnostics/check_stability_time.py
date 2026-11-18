"""Where in the window do the autocorrelations stop resembling each other?"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from etnalib import config as C

Z = np.load(C.PROC / "daily_acf.npz", allow_pickle=True)
lags = Z["lags"]
m = (lags >= 3) & (lags <= 15)

for sta in ("EMSG", "EFIU"):
    A = Z[f"{sta}_acf"]
    days = list(Z[f"{sta}_days"])
    nw = Z[f"{sta}_nwin"]
    print(f"\n=== {sta}: adjacent-day correlation, and windows stacked per day")
    print(f"{'day':12s} {'nwin':>5s} {'adj cc':>7s} {'vs Dec mean':>12s} "
          f"{'rms(coda)':>10s}")
    decmask = np.array([d <= "2025-12-23" for d in days])
    decmean = A[decmask].mean(axis=0)
    for i, d in enumerate(days):
        adj = (np.corrcoef(A[i][m], A[i + 1][m])[0, 1]
               if i + 1 < len(A) else np.nan)
        vs = np.corrcoef(A[i][m], decmean[m])[0, 1]
        rms = np.sqrt(np.mean(A[i][m] ** 2))
        flag = ""
        if np.isfinite(adj) and adj < 0.6:
            flag = "  <-- breaks"
        print(f"{d:12s} {nw[i]:5d} {adj:7.3f} {vs:12.3f} {rms:10.4f}{flag}")
