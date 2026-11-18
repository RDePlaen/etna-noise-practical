"""Render the key figures of each notebook to figures/ for inspection."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import warnings
warnings.filterwarnings("ignore")
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from etnalib import config as C, io, noise, hvsr, dvv, plots

plots.style()

# ---------------------------------------------------------- nb01 spectrogram
days = [d for d in (f"2026-01-{x:02d}" for x in range(14, 27))
        if d in io.available_days("EFIU")]
t, f, S = noise.long_spectrogram("EFIU", days, nperseg=4096, noverlap=2048,
                                 fmax=5.0)
fig, ax = plt.subplots(figsize=(12, 4.2))
k = f > 0.03
m = ax.pcolormesh(t, f[k], 10 * np.log10(S[k]), shading="auto", cmap="magma",
                  vmin=-175, vmax=-105)
ax.set_yscale("log"); ax.set_ylabel("frequency [Hz]")
plots.date_axis(ax); plots.mark_events(ax)
plt.colorbar(m, ax=ax, label="power [dB]", pad=0.01)
ax.set_title("EFIU — the arrival of Cyclone Harry")
fig.savefig(C.FIG / "nb01_harry_spectrogram.png", bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------- nb01 band amplitude
B = np.load(C.PROC / "band_series.npz")
fig, ax = plt.subplots(figsize=(12, 3.8))
for sta in C.STATIONS:
    tk, vk = f"{sta}_secondary_t", f"{sta}_secondary_v"
    if tk not in B:
        continue
    t_, a_ = B[tk].astype("datetime64[s]"), B[vk]
    w = (t_ >= np.datetime64("2026-01-13")) & (t_ <= np.datetime64("2026-01-28"))
    if w.sum() == 0:
        continue
    ax.semilogy(t_[w], a_[w] * 1e6, lw=1.2, color=plots.station_color(sta),
                label=f"{sta} ({C.STATIONS[sta]['coast_km']:.0f} km from coast)")
ax.set_ylabel("0.1–0.35 Hz RMS [µm/s]")
plots.date_axis(ax); plots.mark_events(ax, kinds=("storm",))
ax.legend(fontsize=8); ax.set_title("Secondary microseism, four stations")
fig.savefig(C.FIG / "nb01_microseism_stations.png", bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------ nb02 HVSR
res = {}
for sta in C.STATIONS:
    tr = io.load_hvsr_window(sta)
    if tr is None:
        continue
    r = hvsr.compute_hvsr(*tr, win_s=60.0, fmin=0.2, fmax=25.0, nf=160)
    if r:
        res[sta] = r
fig, axes = plt.subplots(1, len(res), figsize=(3.4 * len(res), 3.8),
                         sharey=True, sharex=True)
for ax, (sta, r) in zip(np.atleast_1d(axes), res.items()):
    ax.fill_between(r["f"], r["hv"] / r["sd"], r["hv"] * r["sd"],
                    color=plots.station_color(sta), alpha=0.25)
    ax.loglog(r["f"], r["hv"], color=plots.station_color(sta), lw=1.8)
    ax.axhline(1.0, color="k", lw=0.6, ls=":")
    if np.isfinite(r["f0"]):
        ax.plot(r["f0"], r["A0"], "v", color="#c1440e", ms=8)
        ax.annotate(f"{r['f0']:.2f} Hz", (r["f0"], r["A0"]),
                    textcoords="offset points", xytext=(4, 6),
                    color="#c1440e", fontsize=9)
    ax.set_title(f"{sta} ({C.STATIONS[sta]['elev']:.0f} m)", fontsize=10)
    ax.set_xlabel("frequency [Hz]")
np.atleast_1d(axes)[0].set_ylabel("H / V")
np.atleast_1d(axes)[0].set_ylim(0.2, 20)
fig.savefig(C.FIG / "nb02_hvsr.png", bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------ nb03 dv/v
Z = np.load(C.PROC / "daily_acf.npz", allow_pickle=True)
lags = Z["lags"]
STA = os.environ.get("DVV_STA", "ECPN")
A, days_ = Z[f"{STA}_acf"], list(Z[f"{STA}_days"])
out, ref, sel = dvv.dvv_series(A, days_, lags, (3.0, 15.0), ref_days=None)
d = np.array([np.datetime64(o[0]) for o in out])
v = np.array([o[1] for o in out]); cc = np.array([o[2] for o in out])
fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 5.6), sharex=True,
                             height_ratios=[2, 1])
g = cc > 0.7
a1.plot(d[g], v[g], "o-", ms=3.5, lw=1, color="#333", label="cc > 0.7")
a1.plot(d[~g], v[~g], "x", ms=5, color="#bbb", label="cc < 0.7")
a1.axhline(0, color="k", lw=0.6); a1.set_ylabel("dv/v [%]")
plots.robust_ylim(a1, v, g, symmetric=True)
a1.legend(fontsize=8, loc="lower left"); plots.mark_events(a1)
a2.plot(d, cc, ".-", ms=3, lw=0.8, color="#1f6f8b")
a2.axhline(0.7, color="#c1440e", ls="--", lw=1)
a2.set_ylabel("cc"); a2.set_ylim(0, 1.02); plots.date_axis(a2)
a1.set_title(f"{STA} — dv/v, 1–3 Hz, coda 3–15 s, grand-mean reference")
fig.savefig(C.FIG / f"nb03_dvv_{STA}.png", bbox_inches="tight")
plt.close(fig)

# ----------------------------------------------------------- nb03 ACF stack
fig, ax = plt.subplots(figsize=(11, 4))
tt = np.array([np.datetime64(x) for x in days_])
mm = np.nanpercentile(np.abs(A[:, lags > 1]), 98)
im = ax.pcolormesh(lags, tt, A, cmap="RdBu_r", vmin=-mm, vmax=mm,
                   shading="auto")
ax.set_xlim(0.5, 20); ax.set_xlabel("lag time [s]")
for date, kind, _s, _l in C.CHRONOLOGY:
    ax.axhline(np.datetime64(date), color=C.COLORS[kind], lw=1.1, ls="--")
plt.colorbar(im, ax=ax, label="autocorrelation", pad=0.01)
ax.set_title(f"{STA} — daily autocorrelations")
fig.savefig(C.FIG / f"nb03_acf_stack_{STA}.png", bbox_inches="tight")
plt.close(fig)

print("figures written to", C.FIG)
for p in sorted(C.FIG.glob("*.png")):
    print(f"  {p.name}  {p.stat().st_size/1e3:.0f} kB")
