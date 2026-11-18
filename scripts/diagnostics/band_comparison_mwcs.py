"""dv/v by MWCS in three bands, following De Plaen et al. (2019).

Method matched to the paper where it is stated:
  PCC autocorrelation, no normalisation      (Front. Earth Sci. 6:251)
  100 Hz -> 20 Hz
  daily autocorrelations stacked over 2 days
  reference = stack of the entire period
  MWCS, lag window 5-35 s
  reject a window if dt error > 0.1 s or coherence < 0.6

The paper's band is 1.0-2.0 Hz; 0.5-1.0 Hz and 1.0-3.0 Hz are run alongside
for comparison. Note the paper uses causal AND acausal lags: for an
AUTOcorrelation those are identical by symmetry, so the non-negative lags
used here are the same measurement.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from etnalib import config as C, mwcs

BANDS = [("1.0-2.0", (1.0, 2.0), "paper"),
         ("0.5-1.0", (0.5, 1.0), ""),
         ("1.0-3.0", (1.0, 3.0), "")]
TMIN, TMAX = 5.0, 35.0
ONSET = "2025-12-24"


def run(tag, band):
    Z = np.load(C.PROC / f"pcc_{tag}.npz", allow_pickle=True)
    lags = Z["lags"]
    fs = 1.0 / (lags[1] - lags[0])
    res = {}
    for sta in C.STATIONS:
        k = f"{sta}_acf"
        if k not in Z:
            continue
        A, days = Z[k], list(Z[f"{sta}_days"])
        ref = A.mean(axis=0)                     # whole-period stack
        v, e, n, c = [], [], [], []
        for a in A:
            dv, se, npts, coh = mwcs.dvv_mwcs(
                a, ref, fs, band[0], band[1], TMIN, TMAX,
                window_length=12.0, step=4.0, max_err=0.1, min_coh=0.6)
            v.append(dv); e.append(se); n.append(npts); c.append(coh)
        res[sta] = dict(days=np.array(days), v=np.array(v, float),
                        e=np.array(e, float), n=np.array(n),
                        coh=np.array(c, float))
    return res


print("=" * 86)
print("dv/v by MWCS on PCC2 autocorrelations, lag 5-35 s, whole-period "
      "reference, 2-day stacks")
print("=" * 86)
print(f"{'band [Hz]':11s} {'station':8s} {'n meas':>7s} {'ok':>5s} "
      f"{'median err':>11s} {'dv/v range':>18s} {'dv/v std':>9s} {'med coh':>8s}")
allres = {}
for tag, band, note in BANDS:
    r = run(tag, band)
    allres[tag] = r
    for sta, d in r.items():
        ok = np.isfinite(d["v"])
        lab = f"{tag}{' *' if note else ''}"
        if ok.sum() == 0:
            print(f"{lab:11s} {sta:8s} {len(d['v']):7d} {0:5d}  no usable measurements")
            continue
        print(f"{lab:11s} {sta:8s} {len(d['v']):7d} {int(ok.sum()):5d} "
              f"{np.nanmedian(d['e'][ok]):11.3f} "
              f"{np.nanmin(d['v'][ok]):+8.3f} ..{np.nanmax(d['v'][ok]):+8.3f} "
              f"{np.nanstd(d['v'][ok]):9.3f} {np.nanmedian(d['coh'][ok]):8.3f}")
    print()

print("=" * 86)
print(f"Across the onset ({ONSET}) -- ECPN, the station in the paper")
print("=" * 86)
hdr = f"{'stack ending':14s}"
for tag, _b, _n in BANDS:
    hdr += f"{tag + ' Hz':>22s}"
print(hdr)
sub = f"{'':14s}" + "".join(f"{'dv/v':>10s}{'+/-':>7s}{'coh':>5s}"
                            for _ in BANDS)
print(sub)
ref_days = None
for tag, _b, _n in BANDS:
    if "ECPN" in allres[tag]:
        ref_days = allres[tag]["ECPN"]["days"]
        break
for day in [d for d in ref_days if "2025-12-14" <= d <= "2026-01-08"]:
    line = f"{day:14s}"
    for tag, _b, _n in BANDS:
        d = allres[tag].get("ECPN")
        i = np.where(d["days"] == day)[0] if d is not None else []
        if len(i) == 0 or not np.isfinite(d["v"][i[0]]):
            line += f"{'-':>10s}{'-':>7s}{'-':>5s}"
        else:
            j = i[0]
            line += f"{d['v'][j]:+10.3f}{d['e'][j]:7.3f}{d['coh'][j]:5.2f}"
    if day >= ONSET and (day <= "2025-12-26"):
        line += "   <-- onset"
    print(line)
