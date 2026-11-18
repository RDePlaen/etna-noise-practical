#!/usr/bin/env python3
"""Spectrogram, PPSD and continuous HVSR for the whole dataset.

The continuous HVSR is new: with HHN/HHE now downloaded for 00:00-06:00 UTC
on every day, H/V can be computed per day across the whole window instead of
once on a single quiet window. That turns a site-characterisation number
into a time series, and the question becomes whether f0 moves -- with the
eruption, with the storm, or at all.

Night hours only, so there is no day/night cycle in here by construction.
"""
from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from etnalib import config as C, io, noise, hvsr, plots  # noqa: E402

plots.style()
FMIN, FMAX, NF = 0.2, 25.0, 140

# Tracking band per station, centred on the resonance found in the single
# quiet window (notebook 02). Picking f0 day after day over the full 0.2-25 Hz
# makes the pick hop between the fundamental and an overtone whenever their
# amplitudes cross, which looks like tens of per cent of site change and is
# nothing of the sort. EFIU has no resonance in band at all, so it gets none
# and is reported as "no f0" rather than given a spurious one.
TRACK = {"ECPN": (0.5, 1.8), "ESLN": (0.35, 1.0), "EMSG": (0.9, 2.5),
         "EFIU": None}
CACHE = C.PROC / "hvsr_continuous.npz"


# ----------------------------------------------------------- continuous HVSR
def load_3c_day(sta, day, fs_out=50.0):
    out = []
    for cha in C.CHANNELS_3C:
        # cache=False: a 24 h, 3-channel, 50 Hz float32 cache over 72 days
        # is tens of GB for data read once.
        tr = io.load_day(sta, day, cha=cha, fs_out=fs_out,
                         remove_response=False, min_coverage=0.3,
                         cache=False)
        if tr is None:
            return None
        out.append(tr)
    n = min(t.stats.npts for t in out)
    for t in out:
        t.data = t.data[:n]
    return out


def build_continuous(stations=None, force=False):
    if CACHE.exists() and not force:
        return dict(np.load(CACHE, allow_pickle=True))
    store = {}
    fcs = np.logspace(np.log10(FMIN), np.log10(FMAX), NF)
    store["f"] = fcs
    for sta in (stations or C.STATIONS):
        days = [d for d in io.available_days(sta, "HHN")
                if d in io.available_days(sta, "HHZ")]
        curves, f0s, used = [], [], []
        for d in days:
            tr = load_3c_day(sta, d)
            if tr is None:
                continue
            r = hvsr.compute_hvsr(*tr, win_s=60.0, overlap=0.5, b=40.0,
                                  fmin=FMIN, fmax=FMAX, nf=NF)
            if r is None:
                continue
            curves.append(r["hv"].astype(np.float32))
            used.append(d)
            tb = TRACK.get(sta)
            if tb is None:
                f0s.append(np.nan)
            else:
                f, _a, _i = hvsr.pick_f0(fcs, r["hv"], fmin=tb[0], fmax=tb[1],
                                         min_amp=1.5, min_prominence=1.1)
                f0s.append(f)
        if not curves:
            print(f"  {sta}: nothing")
            continue
        store[f"{sta}_hv"] = np.vstack(curves)
        store[f"{sta}_f0"] = np.array(f0s)
        store[f"{sta}_days"] = np.array(used)
        f0a = np.array(f0s, float)
        ok = np.isfinite(f0a)
        if ok.sum() == 0:
            print(f"  {sta}: {len(used)} days, no resolvable f0 in band",
                  flush=True)
        else:
            med = np.nanmedian(f0a)
            print(f"  {sta}: {len(used)} days, f0 on {ok.sum()}, "
                  f"median {med:.3f} Hz, "
                  f"scatter {np.nanstd(f0a)/med*100:.1f} % "
                  f"(IQR {np.nanpercentile(f0a,75)-np.nanpercentile(f0a,25):.3f} Hz)",
                  flush=True)
    np.savez_compressed(CACHE, **store)
    return store


def fig_continuous_hvsr(store):
    stas = [s for s in C.STATIONS if f"{s}_hv" in store]
    fig, axes = plt.subplots(len(stas), 1, figsize=(11, 2.5 * len(stas)),
                             sharex=True)
    axes = np.atleast_1d(axes)
    f = store["f"]
    for ax, sta in zip(axes, stas):
        hv = store[f"{sta}_hv"]
        days = np.array([np.datetime64(d) for d in store[f"{sta}_days"]])
        m = ax.pcolormesh(days, f, np.log10(hv.T), shading="auto",
                          cmap="magma", vmin=-0.3, vmax=1.0)
        ax.set_yscale("log")
        ax.set_ylabel(f"{sta}\nfreq [Hz]")
        f0 = store[f"{sta}_f0"]
        ax.plot(days, f0, ".", ms=3, color="#7fe3ff")
        for date, kind, _s, _l in C.CHRONOLOGY:
            d = np.datetime64(date)
            if days.min() <= d <= days.max():
                ax.axvline(d, color=C.COLORS[kind], lw=1.0, ls="--", alpha=.8)
    plt.colorbar(m, ax=list(axes), label="log10( H/V )", pad=0.01)
    plots.date_axis(axes[-1])
    axes[0].set_title("Continuous H/V — one curve per day, 00:00–06:00 UTC "
                      "(cyan dots: picked f0)")
    fig.savefig(C.FIG / "overview_hvsr_continuous.png", bbox_inches="tight")
    plt.close(fig)


def fig_f0_stability(store):
    stas = [s for s in C.STATIONS if f"{s}_f0" in store]
    fig, ax = plt.subplots(figsize=(11, 3.6))
    for sta in stas:
        f0 = store[f"{sta}_f0"]
        if not np.any(np.isfinite(f0)):
            continue
        if np.nanstd(f0) == 0:
            pass
        days = np.array([np.datetime64(d) for d in store[f"{sta}_days"]])
        med = np.nanmedian(f0)
        ax.plot(days, 100 * (f0 - med) / med, ".-", ms=4, lw=0.8,
                color=plots.station_color(sta),
                label=f"{sta}  median f0 = {med:.2f} Hz")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_ylabel("f0 departure from its\nown median  [%]")
    ax.set_ylim(-12, 12)
    plots.mark_events(ax); plots.date_axis(ax)
    ax.legend(fontsize=8, ncol=2, loc="lower left")
    ax.set_title("Is the site resonance stable through an eruption and a cyclone?")
    fig.savefig(C.FIG / "overview_f0_stability.png", bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------- PPSD
def fig_ppsd():
    import pickle
    with open(C.PROC / "ppsd_summary.pkl", "rb") as fh:
        P = pickle.load(fh)
    stas = [s for s in C.STATIONS if s in P]
    fig, axes = plt.subplots(1, len(stas), figsize=(3.3 * len(stas), 4),
                             sharey=True, sharex=True)
    axes = np.atleast_1d(axes)
    vmax = plots.ppsd_vmax([P[s]["hist"] for s in stas])
    for ax, sta in zip(axes, stas):
        d = P[sta]
        m = plots.ppsd_panel(ax, d, vmax)
        ax.set_title(f"{sta}  ({C.STATIONS[sta]['elev']:.0f} m, "
                     f"{plots.ppsd_hours(d):.0f} h)", fontsize=10)
    axes[0].set_ylabel("power [dB rel. 1 (m/s)²/Hz]")
    plt.colorbar(m, ax=list(axes), label="probability", pad=0.01)
    fig.suptitle("PPSD — blue: secondary microseism,  orange: volcanic tremor "
                 "band.  Dashed: Peterson NLNM / NHNM", y=1.02, fontsize=10)
    fig.savefig(C.FIG / "overview_ppsd.png", bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------ spectrogram
def fig_spectrogram(sta="EFIU"):
    days = [d for d in (f"2026-01-{x:02d}" for x in range(14, 27))
            if d in io.available_days(sta)]
    t, f, S = noise.long_spectrogram(sta, days, nperseg=4096, noverlap=2048,
                                     fmax=9.0)
    fig, ax = plt.subplots(figsize=(12, 4.4))
    k = f > 0.03
    m = ax.pcolormesh(t, f[k], 10 * np.log10(S[k]), shading="auto",
                      cmap="magma", vmin=-175, vmax=-105)
    ax.set_yscale("log"); ax.set_ylabel("frequency [Hz]")
    for name, (f1, f2) in C.BANDS.items():
        if f2 > 9:
            continue
        ax.axhline(f1, color="w", lw=0.4, alpha=0.35)
        ax.text(t[int(0.02 * len(t))], np.sqrt(f1 * f2), name, color="w",
                fontsize=8, va="center")
    plots.date_axis(ax); plots.mark_events(ax)
    plt.colorbar(m, ax=ax, label="power [dB]", pad=0.01)
    ax.set_title(f"{sta} ({C.STATIONS[sta]['coast_km']:.0f} km from the coast)"
                 " — Cyclone Harry")
    fig.savefig(C.FIG / "overview_spectrogram.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    what = sys.argv[1:] or ["spec", "ppsd", "hvsr"]
    if "spec" in what:
        print("spectrogram"); fig_spectrogram()
    if "ppsd" in what:
        print("ppsd"); fig_ppsd()
    if "hvsr" in what:
        print("continuous hvsr")
        st = build_continuous(force="--force" in sys.argv)
        fig_continuous_hvsr(st)
        fig_f0_stability(st)
    for p in sorted(C.FIG.glob("overview_*.png")):
        print(f"  {p.name}  {p.stat().st_size/1e3:.0f} kB")
