#!/usr/bin/env python3
"""Export MSNoise daily stacks into the .npz the notebooks read.

MSNoise 1.6.5 writes one miniSEED per station-day-component:

    STACKS/<filter>/001_DAYS/<COMP>/<NET>_<STA>_<NET>_<STA>/<YYYY-MM-DD>.MSEED

with the correlation on a symmetric lag axis from -maxlag to +maxlag.

What this writes
----------------
`processed/daily_acf.npz`:

    lags        non-negative lag axis, seconds           (nlag,)
    lags_full   symmetric lag axis, seconds              (2*nlag-1,)
    <STA>_acf   ZZ autocorrelation, folded               (ndays, nlag)
    <STA>_days  the days, as YYYY-MM-DD                  (ndays,)
    <STA>_nwin  windows stacked per day, or -1           (ndays,)
    <STA>_<CC>  EN / EZ / NZ, two-sided                  (ndays, 2*nlag-1)

**ZZ is folded, the cross-components are not.** An autocorrelation is an
even function of lag, so its two branches are the same measurement twice and
folding them is free signal-to-noise. A single-station cross-component is
NOT even: C_EN(t) = C_NE(-t), so the negative and positive branches describe
different things and folding one would destroy exactly the asymmetry that
makes cross-components interesting.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import obspy

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from etnalib import config as C  # noqa: E402

DATE = re.compile(r"(\d{4}-\d{2}-\d{2})\.MSEED$")


def read_one(path):
    tr = obspy.read(str(path))[0]
    d = tr.data.astype(np.float64)
    fs = float(tr.stats.sampling_rate)
    n = d.size
    if n % 2 == 0:                      # should be odd; guard anyway
        d = d[:-1]
        n -= 1
    lags_full = (np.arange(n) - (n - 1) // 2) / fs
    return lags_full, d


def fold(d):
    """Average the two branches of an even function."""
    half = (d.size - 1) // 2
    pos = d[half:]
    neg = d[: half + 1][::-1]
    return 0.5 * (pos + neg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--filter", default="01")
    ap.add_argument("--out", default=None)
    ap.add_argument("--components", default="ZZ,EN,EZ,NZ")
    a = ap.parse_args()

    root = Path(a.project) / "STACKS" / a.filter / "001_DAYS"
    if not root.exists():
        sys.exit(f"no stacks at {root} -- did compute_cc run?")
    comps = a.components.split(",")
    out = {}
    lags_full = None

    for sta in C.STATIONS:
        per_comp = {}
        for comp in comps:
            d = root / comp / f"{C.NETWORK}_{sta}_{C.NETWORK}_{sta}"
            if not d.exists():
                continue
            rows, days = [], []
            for f in sorted(d.glob("*.MSEED")):
                m = DATE.search(f.name)
                if not m:
                    continue
                try:
                    lf, arr = read_one(f)
                except Exception:
                    continue
                if not np.all(np.isfinite(arr)) or not np.any(arr):
                    continue
                lags_full = lf
                rows.append(arr)
                days.append(m.group(1))
            if rows:
                per_comp[comp] = (np.array(days), np.vstack(rows))

        if "ZZ" not in per_comp:
            print(f"  {sta}: no ZZ, skipped")
            continue

        days, zz = per_comp["ZZ"]
        out[f"{sta}_acf"] = np.vstack([fold(r) for r in zz]).astype(np.float32)
        out[f"{sta}_days"] = days
        out[f"{sta}_nwin"] = np.full(len(days), -1)
        msg = [f"ZZ {len(days)}d"]
        for comp, (dd, arr) in per_comp.items():
            if comp == "ZZ":
                continue
            # keep only days that also have ZZ, so every array lines up
            keep = np.isin(dd, days)
            if keep.sum() == 0:
                continue
            out[f"{sta}_{comp}"] = arr[keep].astype(np.float32)
            out[f"{sta}_{comp}_days"] = dd[keep]
            msg.append(f"{comp} {int(keep.sum())}d")
        print(f"  {sta}: " + ", ".join(msg))

    if lags_full is None:
        sys.exit("nothing exported")
    half = (lags_full.size - 1) // 2
    out["lags"] = lags_full[half:]
    out["lags_full"] = lags_full
    out["band"] = np.array(C.ACORR_BAND)
    # read the normalisation actually used, rather than assuming it: a
    # mislabelled provenance string is worse than none.
    wind = "windsorizing-?"
    try:
        import sqlite3
        db = Path(a.project) / "msnoise.sqlite"
        if db.exists():
            con = sqlite3.connect(db)
            r = con.execute("select value from config where name=?",
                            ("windsorizing",)).fetchone()
            con.close()
            if r:
                wind = {"-1": "one-bit", "0": "none"}.get(
                    str(r[0]).strip(), f"windsorizing-{r[0]}")
    except Exception:
        pass
    out["meta"] = np.array(["msnoise-1.6.5", wind, "linear",
                            str(C.ACORR_WINDOW)])

    dest = Path(a.out) if a.out else (C.PROC / "daily_acf.npz")
    np.savez_compressed(dest, **out)
    print(f"\nwrote {dest} ({dest.stat().st_size/1e6:.1f} MB)")
    print(f"lag axis: 0 .. {out['lags'][-1]:.0f} s "
          f"({out['lags'].size} samples, {1/(out['lags'][1]-out['lags'][0]):.0f} Hz)")


if __name__ == "__main__":
    main()
