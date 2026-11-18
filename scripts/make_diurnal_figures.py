#!/usr/bin/env python3
"""The diurnal cycle: in the noise, in H/V, and in dv/v.

Three things that only became visible once every day was complete:

  1. band amplitudes folded onto the hour of day -- which bands have a
     daily cycle and which do not
  2. hourly H/V, and whether f0 moves with the clock
  3. dv/v measured from night hours only against the whole day
"""
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from etnalib import config as C, hvsr, mwcs, plots  # noqa: E402

plots.style()


# ------------------------------------------------- 1. noise vs hour of day
def fig_bands_diurnal():
    B = np.load(C.PROC / "band_series.npz")
    stas = [s for s in C.STATIONS if f"{s}_secondary_t" in B]
    bands = [b for b in C.BANDS if C.BANDS[b][1] < 9]
    fig, axes = plt.subplots(1, len(bands), figsize=(3.3 * len(bands), 3.6),
                             sharex=True)
    for ax, bname in zip(np.atleast_1d(axes), bands):
        for sta in stas:
            tk, vk = f"{sta}_{bname}_t", f"{sta}_{bname}_v"
            if tk not in B:
                continue
            t = B[tk].astype("datetime64[s]")
            v = B[vk].astype(float)
            hr = ((t - t.astype("datetime64[D]")) /
                  np.timedelta64(1, "h")).astype(float)
            ok = np.isfinite(v) & (v > 0)
            if ok.sum() < 50:
                continue
            med = np.array([np.nanmedian(v[ok][(hr[ok] >= h) & (hr[ok] < h + 1)])
                            for h in range(24)])
            ref = np.nanmedian(med)
            ax.plot(np.arange(24) + 0.5, 100 * (med / ref - 1), "o-", ms=3,
                    lw=1.2, color=plots.station_color(sta), label=sta)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(f"{bname}  {C.BANDS[bname][0]}–{C.BANDS[bname][1]} Hz",
                     fontsize=10)
        ax.set_xlabel("hour of day [UTC]")
        ax.set_xlim(0, 24); ax.set_xticks([0, 6, 12, 18, 24])
        ax.axvspan(0, 6, color="#1f6f8b", alpha=0.10, zorder=0)
    np.atleast_1d(axes)[0].set_ylabel("median amplitude\ndeparture from daily median [%]")
    np.atleast_1d(axes)[-1].legend(fontsize=7, ncol=2)
    fig.suptitle("Diurnal cycle by band — shaded: the 00:00–06:00 window the "
                 "first version of this dataset used", y=1.03, fontsize=10)
    fig.savefig(C.FIG / "diurnal_bands.png", bbox_inches="tight")
    plt.close(fig)
    print("  diurnal_bands.png")


# ---------------------------------------------------------- 2. hourly H/V
def fig_hvsr_diurnal(station=None):
    p = C.PROC / "hvsr_hourly.npz"
    if not p.exists():
        print("  (no hvsr_hourly.npz yet)")
        return
    Z = np.load(p, allow_pickle=True)
    f = Z["f"]
    stas = [s for s in C.STATIONS if f"{s}_hv" in Z]
    if station:
        stas = [station]
    TRACK = {"ECPN": (0.5, 1.8), "ESLN": (0.35, 1.0), "EMSG": (0.9, 2.5),
             "EFIU": None}

    fig, axes = plt.subplots(1, len(stas), figsize=(3.6 * len(stas), 3.8),
                             sharey=True)
    for ax, sta in zip(np.atleast_1d(axes), stas):
        hv = Z[f"{sta}_hv"]
        ts = Z[f"{sta}_t"]
        hr = np.array([int(s[-2:]) for s in ts])
        stack = np.vstack([np.nanmedian(hv[hr == h], axis=0)
                           for h in range(24)])
        # both axes as CENTRES: mixing edges on one axis with centres on
        # the other makes pcolormesh fall back to shading="flat" and raise.
        m = ax.pcolormesh(np.arange(24) + 0.5, f, np.log10(stack.T),
                          shading="auto", cmap="magma", vmin=-0.3, vmax=1.0)
        ax.set_yscale("log")
        ax.set_xlabel("hour of day [UTC]")
        ax.set_xticks([0, 6, 12, 18, 24])
        tb = TRACK.get(sta)
        if tb:
            f0 = []
            for h in range(24):
                v, _a, _i = hvsr.pick_f0(f, stack[h], fmin=tb[0], fmax=tb[1],
                                         min_amp=1.5, min_prominence=1.1)
                f0.append(v)
            ax.plot(np.arange(24) + 0.5, f0, ".-", ms=5, color="#7fe3ff")
            good = np.isfinite(f0)
            if good.sum() > 4:
                spread = 100 * (np.nanmax(f0) - np.nanmin(f0)) / np.nanmedian(f0)
                ax.set_title(f"{sta}   f0 varies {spread:.1f} % over the day",
                             fontsize=9)
        else:
            ax.set_title(f"{sta}   (no resonance in band)", fontsize=9)
    np.atleast_1d(axes)[0].set_ylabel("frequency [Hz]")
    plt.colorbar(m, ax=list(np.atleast_1d(axes)), label="log10( H/V )", pad=0.01)
    fig.suptitle("H/V stacked by hour of day — median of ~70 days per hour",
                 y=1.04, fontsize=10)
    fig.savefig(C.FIG / "diurnal_hvsr.png", bbox_inches="tight")
    plt.close(fig)
    print("  diurnal_hvsr.png")


# ------------------------------------------------------- 3. dv/v, two ways
def fig_dvv_night_vs_day():
    sets = [("night 00–06 UTC", C.PROC / "pcc_1.0-2.0.npz", "#1f6f8b"),
            ("whole day", C.PROC / "pcc_1.0-2.0_24h.npz", "#c1440e")]
    stas = ["ECPN", "EMSG"]
    fig, axes = plt.subplots(len(stas), 1, figsize=(11, 6), sharex=True)
    for ax, sta in zip(np.atleast_1d(axes), stas):
        for lab, path, col in sets:
            if not path.exists():
                continue
            Z = np.load(path, allow_pickle=True)
            if f"{sta}_acf" not in Z:
                continue
            lags = Z["lags"]
            fs = 1.0 / (lags[1] - lags[0])
            A, days = Z[f"{sta}_acf"], Z[f"{sta}_days"]
            ref = A.mean(axis=0)
            v, e = [], []
            for a in A:
                dv, se, n, coh = mwcs.dvv_mwcs(a, ref, fs, 1.0, 2.0, 5.0, 35.0,
                                               window_length=12.0, step=4.0,
                                               max_err=0.1, min_coh=0.6)
                v.append(dv); e.append(se)
            d = np.array([np.datetime64(x) for x in days])
            v = np.array(v); e = np.array(e)
            ok = np.isfinite(v)
            ax.errorbar(d[ok], v[ok], yerr=e[ok], fmt="o-", ms=3.5, lw=1.0,
                        elinewidth=0.8, color=col, label=lab, alpha=0.9)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_ylim(-0.8, 0.8)
        ax.set_ylabel(f"{sta}\ndv/v [%]")
        plots.mark_events(ax, label=(sta == stas[0]), fontsize=7)
    np.atleast_1d(axes)[0].legend(fontsize=8, loc="upper left", ncol=2)
    plots.date_axis(np.atleast_1d(axes)[-1])
    np.atleast_1d(axes)[0].set_title(
        "Same measurement from six night hours and from whole days — "
        "PCC2, 1–2 Hz, MWCS 5–35 s, 2-day stacks", fontsize=10)
    fig.savefig(C.FIG / "diurnal_dvv.png", bbox_inches="tight")
    plt.close(fig)
    print("  diurnal_dvv.png")


if __name__ == "__main__":
    what = sys.argv[1:] or ["bands", "hvsr", "dvv"]
    if "bands" in what:
        fig_bands_diurnal()
    if "hvsr" in what:
        fig_hvsr_diurnal()
    if "dvv" in what:
        fig_dvv_night_vs_day()
