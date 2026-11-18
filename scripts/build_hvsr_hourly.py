#!/usr/bin/env python3
"""Hourly H/V across the whole window, so the diurnal cycle is visible.

The earlier continuous H/V used one curve per day built from 00:00-06:00,
because that was all the horizontal data there was. With complete days we
can compute H/V per HOUR, which turns the same measurement into something
that resolves the daily cycle: 24 curves a day, ~1700 curves in the window.

Two products come out of it:

  hvsr_hourly.npz       curve per station per hour
  the stacked-by-hour-of-day view, built in the figure script

Caching is OFF here. The decimated day traces would be float32 at 50 Hz for
three channels on 72 days, which is tens of gigabytes of cache for data that
is read exactly once.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import obspy

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from etnalib import config as C, hvsr, io  # noqa: E402

FMIN, FMAX, NF = 0.2, 25.0, 140
FS = 50.0                      # Nyquist 25 Hz, matches FMAX


def load_3c_day(sta, day, fs_out=FS):
    out = []
    for cha in C.CHANNELS_3C:
        tr = io.load_day(sta, day, cha=cha, fs_out=fs_out,
                         remove_response=False, min_coverage=0.2,
                         cache=False)
        if tr is None:
            return None
        out.append(tr)
    n = min(t.stats.npts for t in out)
    for t in out:
        t.data = t.data[:n]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", default=",".join(C.STATIONS))
    ap.add_argument("--win", type=float, default=60.0)
    ap.add_argument("--stride", type=int, default=1,
                    help="use every Nth day. The hour-of-day stack averages "
                         "over days, so every 2nd day still gives ~35 "
                         "samples per hour and halves the run time.")
    ap.add_argument("--out", default=str(C.PROC / "hvsr_hourly.npz"))
    a = ap.parse_args()

    fcs = np.logspace(np.log10(FMIN), np.log10(FMAX), NF)
    store = {"f": fcs}
    t0 = time.time()
    for sta in a.stations.split(","):
        days = [d for d in io.available_days(sta, "HHN")
                if d in io.available_days(sta, "HHZ")][::a.stride]
        curves, stamps = [], []
        for d in days:
            tr = load_3c_day(sta, d)
            if tr is None:
                continue
            z, n, e = tr
            day0 = obspy.UTCDateTime(f"{d}T00:00:00")
            for h in range(24):
                sl = [t.slice(day0 + h * 3600.0, day0 + (h + 1) * 3600.0)
                      for t in (z, n, e)]
                if min(s.stats.npts for s in sl) < int(0.6 * 3600 * FS):
                    continue
                r = hvsr.compute_hvsr(*sl, win_s=a.win, overlap=0.5, b=40.0,
                                      fmin=FMIN, fmax=FMAX, nf=NF)
                if r is None:
                    continue
                curves.append(r["hv"].astype(np.float32))
                stamps.append(f"{d}T{h:02d}")
        if not curves:
            print(f"  {sta}: nothing", flush=True)
            continue
        store[f"{sta}_hv"] = np.vstack(curves)
        store[f"{sta}_t"] = np.array(stamps)
        print(f"  {sta}: {len(curves)} hourly curves from {len(days)} days "
              f"[{(time.time()-t0)/60:.1f} min]", flush=True)
    np.savez_compressed(a.out, **store)
    print(f"wrote {a.out} in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
