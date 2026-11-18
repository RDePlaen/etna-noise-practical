#!/usr/bin/env python3
"""Daily PCC2 autocorrelations for the whole window.

No temporal normalisation anywhere in here, on purpose: PCC is
amplitude-unbiased by construction, so one-bit or running-absolute-mean
would be redundant. That is the comparison this script exists to make --
PCC with nothing against CC with one-bit.

Writes processed/daily_acf_pcc.npz in the same layout as daily_acf.npz, so
every downstream tool works on it unchanged.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from etnalib import config as C, io, pcc2  # noqa: E402


def main(stations=None, band=None, maxlag=None, out=None, stack_days=1):
    band = band or C.ACORR_BAND
    maxlag = maxlag or C.ACORR_MAXLAG
    store = {}
    lags = None
    t0 = time.time()
    for sta in (stations or C.STATIONS):
        days = io.available_days(sta)
        acfs, used, nwins = [], [], []
        for d in days:
            tr = io.load_day(sta, d, fs_out=C.ACORR_FS, min_coverage=0.5,
                             remove_response=False)
            if tr is None:
                continue
            lg, c, nw = pcc2.daily_pcc_autocorr(
                tr, band, maxlag, window_s=C.ACORR_WINDOW,
                hours=C.ACORR_HOURS)
            if c is None or nw < 4:
                continue
            lags = lg
            acfs.append(c.astype(np.float32))
            used.append(d)
            nwins.append(nw)
        if not acfs:
            print(f"  {sta}: nothing", flush=True)
            continue
        A = np.vstack(acfs)
        used = np.array(used)
        nwins = np.array(nwins)
        if stack_days > 1:
            # De Plaen et al. (2019) stack the daily autocorrelations over
            # 2 days. Only consecutive calendar days are combined, so a gap
            # never silently merges days that are weeks apart.
            dts = used.astype("datetime64[D]")
            keepA, keepd, keepn = [], [], []
            i = 0
            while i < len(used):
                j = i
                while (j + 1 < len(used) and j - i + 1 < stack_days
                       and (dts[j + 1] - dts[j]).astype(int) == 1):
                    j += 1
                keepA.append(A[i:j + 1].mean(axis=0))
                keepd.append(str(dts[j]))       # label with the last day
                keepn.append(int(nwins[i:j + 1].sum()))
                i = j + 1
            A = np.vstack(keepA)
            used = np.array(keepd)
            nwins = np.array(keepn)
        store[f"{sta}_acf"] = A
        store[f"{sta}_days"] = used
        store[f"{sta}_nwin"] = nwins
        print(f"  {sta}: {len(A)} PCC2 ACFs (from {len(acfs)} days), "
              f"median {int(np.median(nwins))} windows/day "
              f"[{time.time()-t0:.0f} s]", flush=True)
    store["lags"] = lags
    store["band"] = np.array(band)
    store["meta"] = np.array(["pcc2-ventosa2019", "none (amplitude-unbiased)",
                              f"linear x{stack_days}d", str(C.ACORR_WINDOW)])
    dest = Path(out) if out else (C.PROC / "daily_acf_pcc.npz")
    np.savez_compressed(dest, **store)
    print(f"wrote {dest} ({dest.stat().st_size/1e6:.1f} MB) "
          f"in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default=None, help="e.g. 1.0,2.0")
    ap.add_argument("--maxlag", type=float, default=None)
    ap.add_argument("--stack-days", type=int, default=1)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    b = tuple(float(x) for x in a.band.split(",")) if a.band else None
    main(band=b, maxlag=a.maxlag, out=a.out, stack_days=a.stack_days)
