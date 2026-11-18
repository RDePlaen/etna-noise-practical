"""PCC2 against one-bit CC, focused on the 24 December 2025 eruption onset.

The question is specific: the autocorrelations were already known to drift
away from their December shape starting 24-27 December, which is the onset
of the continuous eruptive phase, and with classical CC that drift destroys
any dv/v measured against a pre-eruption reference. PCC is
amplitude-unbiased, and the thing that changes at the onset is above all the
AMPLITUDE of the summit tremor. So PCC ought to survive where CC did not.

That is a prediction, so measure it.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from etnalib import config as C, dvv

ONSET = "2025-12-24"
PRE = ("2025-12-10", "2025-12-23")      # the 14 days before the eruption
LAGWIN = (3.0, 15.0)

SETS = [("CC  one-bit", C.PROC / "daily_acf.npz"),
        ("PCC2  none ", C.PROC / "daily_acf_pcc.npz")]


def series(path, sta, ref):
    Z = np.load(path, allow_pickle=True)
    k = f"{sta}_acf"
    if k not in Z:
        return None
    A, days = Z[k], list(Z[f"{sta}_days"])
    out, _, sel = dvv.dvv_series(A, days, Z["lags"], LAGWIN, ref_days=ref)
    d = np.array([o[0] for o in out])
    v = np.array([o[1] for o in out])
    cc = np.array([o[2] for o in out])
    return d, v, cc, A, Z["lags"], sel


print("=" * 78)
print("1. Day-to-day stability of the autocorrelations themselves")
print("   (adjacent-day correlation, 3-15 s coda; no reference involved)")
print("=" * 78)
print(f"{'station':9s} {'estimator':12s} {'all days':>9s} "
      f"{'before onset':>13s} {'after onset':>12s}")
for sta in C.STATIONS:
    for lab, p in SETS:
        r = series(p, sta, None)
        if r is None:
            continue
        d, v, cc, A, lags, _ = r
        m = (lags >= LAGWIN[0]) & (lags <= LAGWIN[1])
        adj = np.array([np.corrcoef(A[i][m], A[i + 1][m])[0, 1]
                        for i in range(len(A) - 1)])
        dd = np.array(d[:-1])
        pre = adj[dd < ONSET]
        post = adj[dd >= ONSET]
        print(f"{sta:9s} {lab:12s} {np.median(adj):9.3f} "
              f"{(np.median(pre) if pre.size else np.nan):13.3f} "
              f"{(np.median(post) if post.size else np.nan):12.3f}")

print()
print("=" * 78)
print("2. dv/v against a PRE-ERUPTION reference (10-23 Dec)")
print("   This is the measurement that matters for the onset, and the one")
print("   that classical CC could not deliver.")
print("=" * 78)
print(f"{'station':9s} {'estimator':12s} {'median cc':>10s} {'cc>0.7':>8s} "
      f"{'usable after onset':>19s} {'dv/v std':>9s}")
for sta in C.STATIONS:
    for lab, p in SETS:
        r = series(p, sta, PRE)
        if r is None:
            continue
        d, v, cc, A, lags, sel = r
        post = np.array(d) >= ONSET
        good = cc > 0.7
        print(f"{sta:9s} {lab:12s} {np.nanmedian(cc):10.3f} "
              f"{int(good.sum()):4d}/{len(cc):<3d} "
              f"{int((good & post).sum()):10d}/{int(post.sum()):<8d} "
              f"{np.nanstd(v[good]) if good.any() else np.nan:9.3f}")

print()
print("=" * 78)
print(f"3. What happens across the onset itself  (reference = {PRE[0]}..{PRE[1]})")
print("=" * 78)
for sta in ("ECPN", "EMSG"):
    print(f"\n--- {sta}")
    print(f"{'day':12s} {'CC dv/v':>9s} {'cc':>6s}   {'PCC dv/v':>9s} {'cc':>6s}")
    rc = series(SETS[0][1], sta, PRE)
    rp = series(SETS[1][1], sta, PRE)
    if rc is None or rp is None:
        continue
    dc, vc, ccc, *_ = rc
    dp, vp, ccp, *_ = rp
    days = [x for x in dp if "2025-12-18" <= x <= "2026-01-06"]
    for day in days:
        ic = np.where(dc == day)[0]
        ip = np.where(dp == day)[0]
        a = f"{vc[ic[0]]:+9.3f} {ccc[ic[0]]:6.3f}" if ic.size else f"{'-':>9s} {'-':>6s}"
        b = f"{vp[ip[0]]:+9.3f} {ccp[ip[0]]:6.3f}" if ip.size else f"{'-':>9s} {'-':>6s}"
        mark = "  <-- onset" if day == ONSET else ""
        print(f"{day:12s} {a}   {b}{mark}")
