import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from etnalib import dvv, config as C

Z = np.load(C.PROC / "daily_acf.npz", allow_pickle=True)
lags = Z["lags"]
REF = None   # grand mean of all days

for sta in C.STATIONS:
    k = f"{sta}_acf"
    if k not in Z:
        continue
    acfs, days = Z[k], list(Z[f"{sta}_days"])
    out, ref, sel = dvv.dvv_series(acfs, days, lags, (3.0, 15.0),
                                   ref_days=REF)
    v = np.array([o[1] for o in out])
    cc = np.array([o[2] for o in out])
    g = cc > 0.7
    print(f"\n=== {sta}: {len(days)} days, reference from {sel.sum()} days, "
          f"{g.sum()} days with cc>0.7")
    print(f"    median cc {np.nanmedian(cc):.3f}   "
          f"dv/v range {np.nanmin(v[g]) if g.any() else np.nan:+.3f} .. "
          f"{np.nanmax(v[g]) if g.any() else np.nan:+.3f} %   "
          f"std {np.nanstd(v[g]) if g.any() else np.nan:.3f} %")
    # weekly means, to see the shape without the noise
    print("    weekly means (cc>0.7 only):")
    for i in range(0, len(days), 7):
        dd = days[i:i + 7]
        vv = v[i:i + 7][g[i:i + 7]]
        if len(vv) == 0:
            print(f"      {dd[0]} .. {dd[-1]}:   (no usable days)")
            continue
        print(f"      {dd[0]} .. {dd[-1]}: {np.nanmean(vv):+6.3f} %  "
              f"(n={len(vv)}, cc={np.nanmean(cc[i:i+7][g[i:i+7]]):.2f})")
