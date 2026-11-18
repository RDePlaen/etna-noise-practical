#!/usr/bin/env python3
"""dv/v across the 24 December 2025 eruption onset: PCC2 against one-bit CC."""
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
from etnalib import config as C, dvv, plots  # noqa: E402

plots.style()
ONSET = np.datetime64("2025-12-24")
PRE = ("2025-12-10", "2025-12-23")
LAGWIN = (3.0, 15.0)
SETS = [("CC, one-bit", C.PROC / "daily_acf.npz", "#333333"),
        ("PCC2, no normalisation", C.PROC / "daily_acf_pcc.npz", "#c1440e")]
STAS = ["ECPN", "EMSG"]


def get(path, sta, ref):
    Z = np.load(path, allow_pickle=True)
    k = f"{sta}_acf"
    if k not in Z:
        return None
    out, _, _ = dvv.dvv_series(Z[k], list(Z[f"{sta}_days"]), Z["lags"],
                               LAGWIN, ref_days=ref)
    d = np.array([np.datetime64(o[0]) for o in out])
    return d, np.array([o[1] for o in out]), np.array([o[2] for o in out])


fig, axes = plt.subplots(4, 1, figsize=(11, 9), sharex=True,
                         height_ratios=[2, 1, 2, 1])
for row, sta in enumerate(STAS):
    a_v, a_c = axes[2 * row], axes[2 * row + 1]
    for lab, path, col in SETS:
        r = get(path, sta, PRE)
        if r is None:
            continue
        d, v, cc = r
        good = cc > 0.7
        a_v.plot(d[good], v[good], "o-", ms=4, lw=1.1, color=col, label=lab)
        a_v.plot(d[~good], v[~good], "x", ms=4, color=col, alpha=0.3)
        a_c.plot(d, cc, ".-", ms=4, lw=1.0, color=col, label=lab)
    a_v.axhline(0, color="k", lw=0.6)
    a_v.set_ylim(-0.6, 0.6)
    a_v.set_ylabel(f"{sta}\ndv/v  [%]")
    a_c.axhline(0.7, color="#1f6f8b", ls="--", lw=1.0)
    a_c.set_ylim(0, 1.02)
    a_c.set_ylabel("cc")
    for ax in (a_v, a_c):
        ax.axvspan(np.datetime64(PRE[0]), np.datetime64(PRE[1]),
                   color="#1f6f8b", alpha=0.07, zorder=0)
        ax.axvline(ONSET, color=C.COLORS["volcanic"], lw=1.6, ls="--")
    if row == 0:
        a_v.legend(fontsize=8, loc="lower left")
        a_v.annotate("eruptive phase begins\n24 Dec",
                     xy=(ONSET, 0.52), xycoords=("data", "data"),
                     fontsize=8, color=C.COLORS["volcanic"], ha="right")
        a_v.annotate("reference window", xy=(np.datetime64("2025-12-16"), -0.5),
                     fontsize=8, color="#1f6f8b", ha="center")
axes[-1].set_xlim(np.datetime64("2025-12-10"), np.datetime64("2026-01-10"))
plots.date_axis(axes[-1])
axes[0].set_title("dv/v across the 24 December eruption onset, measured against "
                  "a pre-eruption reference\n"
                  "crosses = rejected (cc < 0.7);  1–3 Hz, coda 3–15 s",
                  fontsize=10)
fig.savefig(C.FIG / "pcc_vs_cc_december.png", bbox_inches="tight")
plt.close(fig)

# ---- the autocorrelations themselves, before and after
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
for ax, (lab, path, col) in zip(axes, SETS):
    Z = np.load(path, allow_pickle=True)
    lags = Z["lags"]
    A, days = Z["EMSG_acf"], np.array(Z["EMSG_days"])
    pre = A[days <= "2025-12-23"].mean(axis=0)
    post = A[(days >= "2025-12-25") & (days <= "2026-01-07")].mean(axis=0)
    ax.plot(lags, pre, lw=1.4, color="#1f6f8b", label="10–23 Dec (before)")
    ax.plot(lags, post, lw=1.4, color=C.COLORS["volcanic"],
            label="25 Dec – 7 Jan (after)")
    m = (lags >= LAGWIN[0]) & (lags <= LAGWIN[1])
    ax.set_xlim(0, 20)
    ax.set_xlabel("lag [s]")
    ax.set_title(f"{lab}   (coda cc = {np.corrcoef(pre[m], post[m])[0,1]:.3f})",
                 fontsize=10)
    ax.axvspan(LAGWIN[0], LAGWIN[1], color="k", alpha=0.05, zorder=0)
axes[0].set_ylabel("autocorrelation")
axes[0].legend(fontsize=8)
fig.suptitle("EMSG — stacked autocorrelation before and after the onset "
             "(shaded: the measurement window)", y=1.03, fontsize=10)
fig.savefig(C.FIG / "pcc_vs_cc_waveforms.png", bbox_inches="tight")
plt.close(fig)

for p in sorted(C.FIG.glob("pcc_*.png")):
    print(f"  {p.name}  {p.stat().st_size/1e3:.0f} kB")
