"""Night-only (00-06 UTC) against whole days, same everything else.

The archive originally held only night hours on most days, so the
autocorrelations were built from 00:00-06:00. Now that every day is
complete the same measurement can be made from all 24 hours, and the two
can be compared directly: same band (1-2 Hz), same 2-day stacks, same
whole-period reference, same MWCS window (5-35 s).

Four times more data per stack should reduce the scatter. Whether it does
is the question -- the extra hours are daytime, which is noisier and whose
noise field points in different directions.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from etnalib import config as C, mwcs

SETS = [("night 00-06", C.PROC / "pcc_1.0-2.0.npz"),
        ("whole day  ", C.PROC / "pcc_1.0-2.0_24h.npz")]
BAND, TMIN, TMAX = (1.0, 2.0), 5.0, 35.0


def run(path):
    Z = np.load(path, allow_pickle=True)
    lags = Z["lags"]
    fs = 1.0 / (lags[1] - lags[0])
    out = {}
    for sta in C.STATIONS:
        k = f"{sta}_acf"
        if k not in Z:
            continue
        A, days = Z[k], np.array(Z[f"{sta}_days"])
        ref = A.mean(axis=0)
        v, e, c = [], [], []
        for a in A:
            dv, se, n, coh = mwcs.dvv_mwcs(a, ref, fs, BAND[0], BAND[1],
                                           TMIN, TMAX, window_length=12.0,
                                           step=4.0, max_err=0.1,
                                           min_coh=0.6)
            v.append(dv); e.append(se); c.append(coh)
        # day-to-day stability of the autocorrelations themselves
        m = (lags >= TMIN) & (lags <= TMAX)
        adj = np.median([np.corrcoef(A[i][m], A[i + 1][m])[0, 1]
                         for i in range(len(A) - 1)])
        out[sta] = dict(days=days, v=np.array(v), e=np.array(e),
                        coh=np.array(c), adj=adj, nwin=Z[f"{sta}_nwin"])
    return out


res = {lab: run(p) for lab, p in SETS}

print(f"{'station':8s} {'coverage':12s} {'win/stack':>9s} {'adj-day cc':>11s} "
      f"{'ok':>5s} {'median err':>11s} {'dv/v std':>9s} {'med coh':>8s}")
for sta in C.STATIONS:
    for lab, _ in SETS:
        d = res[lab].get(sta)
        if d is None:
            continue
        ok = np.isfinite(d["v"])
        print(f"{sta:8s} {lab:12s} {int(np.median(d['nwin'])):9d} "
              f"{d['adj']:11.3f} {int(ok.sum()):3d}/{len(d['v']):<3d} "
              f"{np.nanmedian(d['e'][ok]):11.3f} "
              f"{np.nanstd(d['v'][ok]):9.3f} "
              f"{np.nanmedian(d['coh'][ok]):8.3f}")
    print()

print("Do the two agree, stack by stack?  (1-2 Hz, common dates)")
print(f"{'station':8s} {'n common':>9s} {'correlation':>12s} "
      f"{'median |diff|':>14s} {'rms diff':>9s}")
for sta in C.STATIONS:
    a, b = res[SETS[0][0]].get(sta), res[SETS[1][0]].get(sta)
    if a is None or b is None:
        continue
    common = np.intersect1d(a["days"], b["days"])
    ia = np.array([np.where(a["days"] == d)[0][0] for d in common])
    ib = np.array([np.where(b["days"] == d)[0][0] for d in common])
    va, vb = a["v"][ia], b["v"][ib]
    k = np.isfinite(va) & np.isfinite(vb)
    if k.sum() < 3:
        continue
    print(f"{sta:8s} {int(k.sum()):9d} "
          f"{np.corrcoef(va[k], vb[k])[0,1]:12.3f} "
          f"{np.median(np.abs(va[k]-vb[k])):14.3f} "
          f"{np.sqrt(np.mean((va[k]-vb[k])**2)):9.3f}")
