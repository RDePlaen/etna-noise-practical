"""Which coda window keeps enough coherence with a December reference?

Late coda samples more medium and is more sensitive to a velocity change,
but it also decorrelates first. Find the window that is still usable across
the whole three months.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from etnalib import dvv, config as C

Z = np.load(C.PROC / "daily_acf.npz", allow_pickle=True)
lags = Z["lags"]
REF = ("2025-12-10", "2025-12-23")
WINS = [(1, 4), (1, 6), (2, 8), (3, 10), (3, 15), (5, 20), (8, 20)]

for sta in ("EMSG", "EFIU", "ECPN"):
    A, days = Z[f"{sta}_acf"], list(Z[f"{sta}_days"])
    print(f"\n=== {sta} ({len(days)} days)")
    print(f"{'lag window':12s} {'median cc':>10s} {'days cc>0.7':>12s} "
          f"{'days cc>0.5':>12s} {'dv/v std':>9s} {'dv/v range':>18s}")
    for w in WINS:
        out, ref, sel = dvv.dvv_series(A, days, lags, w, ref_days=REF)
        v = np.array([o[1] for o in out])
        cc = np.array([o[2] for o in out])
        g = cc > 0.5
        print(f"{str(w):12s} {np.nanmedian(cc):10.3f} {int((cc>0.7).sum()):12d} "
              f"{int(g.sum()):12d} {np.nanstd(v[g]) if g.any() else np.nan:9.3f} "
              f"{np.nanmin(v[g]) if g.any() else np.nan:+8.3f} ..{np.nanmax(v[g]) if g.any() else np.nan:+8.3f}")
