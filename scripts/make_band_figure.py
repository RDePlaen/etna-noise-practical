#!/usr/bin/env python3
"""dv/v by MWCS in three bands, De Plaen et al. (2019) method."""
import sys, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from etnalib import config as C, mwcs, plots  # noqa: E402

plots.style()
BANDS = [("1.0-2.0", (1.0, 2.0), "#c1440e", "1–2 Hz  (De Plaen et al. 2019)"),
         ("0.5-1.0", (0.5, 1.0), "#1f6f8b", "0.5–1 Hz"),
         ("1.0-3.0", (1.0, 3.0), "#555555", "1–3 Hz")]
TMIN, TMAX = 5.0, 35.0
STAS = ["ECPN", "EMSG"]

fig, axes = plt.subplots(len(STAS), 1, figsize=(11, 6.4), sharex=True)
for ax, sta in zip(np.atleast_1d(axes), STAS):
    for tag, band, col, lab in BANDS:
        Z = np.load(C.PROC / f"pcc_{tag}.npz", allow_pickle=True)
        if f"{sta}_acf" not in Z:
            continue
        lags = Z["lags"]
        fs = 1.0 / (lags[1] - lags[0])
        A, days = Z[f"{sta}_acf"], list(Z[f"{sta}_days"])
        ref = A.mean(axis=0)
        v, e = [], []
        for a in A:
            dv, se, n, coh = mwcs.dvv_mwcs(a, ref, fs, band[0], band[1],
                                           TMIN, TMAX, window_length=12.0,
                                           step=4.0, max_err=0.1, min_coh=0.6)
            v.append(dv); e.append(se)
        d = np.array([np.datetime64(x) for x in days])
        v = np.array(v); e = np.array(e)
        ok = np.isfinite(v)
        ax.errorbar(d[ok], v[ok], yerr=e[ok], fmt="o-", ms=3.5, lw=1.0,
                    elinewidth=0.8, capsize=0, color=col, label=lab, alpha=0.9)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_ylim(-1.0, 1.0)
    ax.set_ylabel(f"{sta}\ndv/v  [%]")
    plots.mark_events(ax, label=(sta == STAS[0]), fontsize=7)
np.atleast_1d(axes)[0].legend(fontsize=8, ncol=3, loc="upper left")
plots.date_axis(np.atleast_1d(axes)[-1])
np.atleast_1d(axes)[0].set_title(
    "dv/v by MWCS on PCC2 autocorrelations — 2-day stacks, whole-period "
    "reference, coda 5–35 s\nerror bars are the MWCS regression error",
    fontsize=10)
fig.savefig(C.FIG / "dvv_bands_mwcs.png", bbox_inches="tight")
plt.close(fig)
print("  dvv_bands_mwcs.png")
