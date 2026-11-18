#!/usr/bin/env python3
"""Build everything the practical loads from disk, from the raw archive.

Run once, offline. Students never run this: they get the output in the data
archive. Products:

  processed/ppsd_summary.pkl   PPSD histograms + Peterson models, per station
  processed/daily_acf.npz      daily autocorrelations, per station
  processed/band_series.npz    band amplitudes through the whole window
  processed/hvsr_ready.npz     the three-component HVSR window, decimated
"""
import pickle
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from etnalib import acorr, config as C, io, noise   # noqa: E402


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


# ------------------------------------------------------------------ PPSD
def build_ppsd(stations=None, max_days=18):
    from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
    out = {}
    for sta in (stations or C.STATIONS):
        days = io.available_days(sta)
        if not days:
            continue
        sel = days[:: max(1, len(days) // max_days)][:max_days]
        p = None
        for d in sel:
            try:
                p = noise.ppsd_for_day(sta, d, ppsd=p)
            except Exception as e:
                log(f"  {sta} {d}: {e}")
        if p is None or len(p.times_processed) == 0:
            log(f"  {sta}: no PPSD")
            continue
        hist, db = p.current_histogram, p.db_bin_centers
        T = p.period_bin_centers
        nl_T, nl = get_nlnm()
        nh_T, nh = get_nhnm()
        out[sta] = dict(periods=T, db=db,
                        hist=hist / max(p.current_histogram_count, 1),
                        n=len(p.times_processed), days=sel,
                        nlnm_T=nl_T, nlnm=nl, nhnm_T=nh_T, nhnm=nh)
        log(f"  {sta}: PPSD from {len(p.times_processed)} segments, {len(sel)} days")
    with open(C.PROC / "ppsd_summary.pkl", "wb") as fh:
        pickle.dump(out, fh)
    return out


# ------------------------------------------------------------- daily ACF
def build_acf(stations=None, band=None, maxlag=None, estimator="ccgn",
              norm="ram", stack="linear"):
    band = band or C.ACORR_BAND
    maxlag = maxlag or C.ACORR_MAXLAG
    store = {}
    for sta in (stations or C.STATIONS):
        days = io.available_days(sta)
        acfs, used, nwins = [], [], []
        lags = None
        for d in days:
            # remove_response=False on purpose: for a SINGLE-station
            # autocorrelation the instrument phase cancels exactly
            # (C(f) = |H(f)|^2 S(f)), so deconvolution buys nothing and
            # costs minutes per day. See io.load_day's docstring.
            tr = io.load_day(sta, d, fs_out=C.ACORR_FS, min_coverage=0.5,
                             remove_response=False)
            if tr is None:
                continue
            lg, c, nw = acorr.daily_autocorr(
                tr, band, maxlag, window_s=C.ACORR_WINDOW,
                estimator=estimator, norm=norm, stack=stack,
                hours=C.ACORR_HOURS)     # None = whole day; see config
            if c is None or nw < 4:   # a day stacked from <4 windows is noise
                continue
            lags = lg
            acfs.append(c.astype(np.float32))
            used.append(d)
            nwins.append(nw)
        if not acfs:
            log(f"  {sta}: no ACFs")
            continue
        store[f"{sta}_acf"] = np.vstack(acfs)
        store[f"{sta}_days"] = np.array(used)
        store[f"{sta}_nwin"] = np.array(nwins)
        log(f"  {sta}: {len(acfs)} daily ACFs, median {int(np.median(nwins))} windows/day")
    store["lags"] = lags
    store["band"] = np.array(band)
    store["meta"] = np.array([estimator, norm, stack, str(C.ACORR_WINDOW)])
    np.savez_compressed(C.PROC / "daily_acf.npz", **store)
    return store


# ----------------------------------------------------------- band series
# Band amplitudes need the instrument response removed (counts/sensitivity
# is within 1 % for the Trillium stations but 10 % out at EFIU, measured in
# scripts/check_sensitivity_shortcut.py), and deconvolution costs ~21 s per
# station-day. So build only what the notebooks actually plot:
#   notebook 01: all four stations, but only across the Cyclone Harry window
#   notebook 04: ECPN and EFIU, across the whole window
BAND_SCOPE = {
    "ECPN": None,                              # None = every available day
    "EFIU": None,
    "ESLN": ("2026-01-13", "2026-01-28"),
    "EMSG": ("2026-01-13", "2026-01-28"),
}


def build_bands(stations=None, scope=None):
    scope = BAND_SCOPE if scope is None else scope
    store = {}
    for sta in (stations or C.STATIONS):
        days = io.available_days(sta)
        win = scope.get(sta) if scope else None
        if win:
            days = [d for d in days if win[0] <= d <= win[1]]
        if not days:
            continue
        s = noise.band_series(sta, days)
        for k, (t, v) in s.items():
            if len(t) == 0:
                continue
            store[f"{sta}_{k}_t"] = t.astype("datetime64[s]").astype(np.int64)
            store[f"{sta}_{k}_v"] = v.astype(np.float32)
        log(f"  {sta}: band series over {len(days)} days")
    np.savez_compressed(C.PROC / "band_series.npz", **store)
    return store


# ------------------------------------------------------------------ HVSR
def build_hvsr_ready(stations=None):
    store = {}
    for sta in (stations or C.STATIONS):
        tr = io.load_hvsr_window(sta)
        if tr is None:
            log(f"  {sta}: no HVSR window")
            continue
        z, n, e = tr
        for lab, t in zip("ZNE", (z, n, e)):
            store[f"{sta}_{lab}"] = t.data.astype(np.float32)
        store[f"{sta}_fs"] = np.array([z.stats.sampling_rate])
        log(f"  {sta}: HVSR window {z.stats.npts} pts @ {z.stats.sampling_rate} Hz")
    np.savez_compressed(C.PROC / "hvsr_ready.npz", **store)
    return store


if __name__ == "__main__":
    what = sys.argv[1:] or ["hvsr", "bands", "acf", "ppsd"]
    if "hvsr" in what:
        log("HVSR window"); build_hvsr_ready()
    if "bands" in what:
        log("band series"); build_bands()
    if "acf" in what:
        log("daily autocorrelations"); build_acf()
    if "ppsd" in what:
        log("PPSD"); build_ppsd()
    log("done")
