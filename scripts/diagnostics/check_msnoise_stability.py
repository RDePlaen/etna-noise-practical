"""Why are MSNoise's daily autocorrelations less stable than etnalib's?

Measures adjacent-day correlation, which separates "the waveform drifts"
from "each day is noisy".
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from etnalib import config as C

Z = np.load(C.PROC / "daily_acf.npz", allow_pickle=True)
lags = Z["lags"]
print("produced by:", list(Z["meta"]))

for win in [(1, 6), (3, 15), (5, 20)]:
    m = (lags >= win[0]) & (lags <= win[1])
    print(f"\n--- coda {win[0]}-{win[1]} s")
    for sta in C.STATIONS:
        k = f"{sta}_acf"
        if k not in Z:
            continue
        A = Z[k]
        adj = [np.corrcoef(A[i][m], A[i + 1][m])[0, 1]
               for i in range(len(A) - 1)]
        gm = A.mean(axis=0)
        vs = [np.corrcoef(a[m], gm[m])[0, 1] for a in A]
        # how much of each day's energy sits at very short lag? a daily ACF
        # dominated by the zero-lag spike has little usable coda
        e_coda = np.sqrt(np.mean(A[:, m] ** 2, axis=1))
        e_near = np.sqrt(np.mean(A[:, lags < 1] ** 2, axis=1))
        print(f"  {sta:6s} adjacent {np.median(adj):6.3f}   "
              f"vs grand mean {np.median(vs):6.3f}   "
              f"coda/near-lag rms {np.median(e_coda/e_near):7.4f}")
