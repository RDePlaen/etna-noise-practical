"""Shared plotting helpers, so every figure in the practical reads the same way.

The one rule worth stating: volcanic events are drawn in one colour and
meteorological events in another, on every single time-series figure. The
whole point of this dataset is that both kinds of event are present in the
same window, and a student who sees them in the same colour will not ask the
question the practical is built around.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import dates as mdates

from . import config as C


def style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 160,
        "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "figure.autolayout": False,
    })


def mark_events(ax, kinds=("volcanic", "storm"), label=True, y=0.985,
                fontsize=7.0, clip=True):
    """Draw the chronology on a time axis.

    `clip=True` keeps the existing x limits: an event line outside the data
    would otherwise stretch the axis and squash the thing being plotted into
    a corner. Call this AFTER plotting, not before.
    """
    import matplotlib.dates as _md
    x0, x1 = ax.get_xlim()
    lo, hi = _md.num2date(x0).replace(tzinfo=None), \
        _md.num2date(x1).replace(tzinfo=None)
    for date, kind, short, _long in C.CHRONOLOGY:
        if kind not in kinds:
            continue
        d = np.datetime64(date)
        if clip and not (np.datetime64(lo) <= d <= np.datetime64(hi)):
            continue
        ax.axvline(d, color=C.COLORS[kind], lw=1.1, ls="--", alpha=0.85,
                   zorder=0)
        if label:
            # outline the text: these labels often sit over a bright
            # spectrogram, where plain coloured text is unreadable
            import matplotlib.patheffects as pe
            ax.annotate(short, xy=(d, y), xycoords=("data", "axes fraction"),
                        rotation=90, va="top", ha="right", fontsize=fontsize,
                        color=C.COLORS[kind], alpha=1.0, zorder=5,
                        path_effects=[pe.withStroke(linewidth=2.2,
                                                    foreground="white")])
    if clip:
        ax.set_xlim(x0, x1)


def robust_ylim(ax, values, mask=None, pad=0.35, symmetric=False):
    """Set y limits from the GOOD points, not from the outliers.

    A single rejected day at -4.5 % will otherwise compress a +/-0.3 %
    signal into one pixel row, which is how a real measurement gets mistaken
    for a flat line.
    """
    v = np.asarray(values, float)
    if mask is not None:
        v = v[np.asarray(mask, bool)]
    v = v[np.isfinite(v)]
    if v.size == 0:
        return
    lo, hi = np.nanmin(v), np.nanmax(v)
    span = max(hi - lo, 1e-9)
    lo, hi = lo - pad * span, hi + pad * span
    if symmetric:
        m = max(abs(lo), abs(hi))
        lo, hi = -m, m
    ax.set_ylim(lo, hi)


def event_legend(ax, **kw):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=C.COLORS["volcanic"], ls="--", label="volcanic"),
         Line2D([], [], color=C.COLORS["storm"], ls="--", label="storm")]
    ax.legend(handles=h, loc="upper left", fontsize=8, ncol=2, **kw)


def date_axis(ax, fmt="%d %b"):
    ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=4, maxticks=10))
    ax.xaxis.set_major_formatter(mdates.DateFormatter(fmt))


def break_gaps(t, v, max_gap_h=2.0):
    """Insert NaN where the series jumps by more than `max_gap_h` hours.

    matplotlib joins consecutive points with a straight line regardless of
    how far apart they are in time, so a multi-day outage is drawn as a
    smooth ramp between the last sample before it and the first after. On a
    tremor amplitude plot that ramp looks exactly like a slow eruption
    onset. Breaking the line at real gaps is not cosmetic.
    """
    t = np.asarray(t)
    v = np.asarray(v, dtype=float)
    if t.size < 2:
        return t, v
    dt = np.diff(t).astype("timedelta64[s]").astype(float) / 3600.0
    bad = np.where(dt > max_gap_h)[0]
    if bad.size == 0:
        return t, v
    # insert a NaN sample just after each gap start
    t_out = np.insert(t, bad + 1, t[bad] + np.timedelta64(1, "s"))
    v_out = np.insert(v, bad + 1, np.nan)
    return t_out, v_out


def ppsd_vmax(hists, pct=97.0):
    """Colour-scale maximum for a PPSD.

    ObsPy's histogram is normalised so each PERIOD COLUMN sums to 1, spread
    over a few hundred dB bins, so occupied cells hold ~0.01-0.1 and the
    maximum is 1.0 only where a column is in a single bin. Leaving vmax at
    its default renders the whole thing as blank paper. Take a high
    percentile of the occupied cells instead.
    """
    a = np.concatenate([np.asarray(h).ravel() for h in hists])
    a = a[a > 0]
    return float(np.percentile(a, pct)) if a.size else 1.0


def ppsd_panel(ax, d, vmax, bands=True):
    """One PPSD panel, drawn the same way everywhere.

    `d` is an entry of processed/ppsd_summary.pkl. Note `d['n']` counts
    1800 s segments at 50 % overlap, so the time span is n*900/3600 hours --
    labelling n as hours overstates the data about fourfold.
    """
    m = ax.pcolormesh(d["periods"], d["db"], d["hist"].T, cmap="YlOrRd",
                      shading="auto", vmin=0, vmax=vmax)
    ax.plot(d["nlnm_T"], d["nlnm"], "k--", lw=1)
    ax.plot(d["nhnm_T"], d["nhnm"], "k--", lw=1)
    ax.set_xscale("log")
    ax.set_xlim(0.05, 50)
    ax.set_ylim(-190, -90)
    if bands:
        # vertical lines, not axvspan: a span drawn under an opaque
        # pcolormesh is invisible, and drawn after it hides the data
        for (f1, f2), col, lab in (((0.10, 0.35), "#1f6f8b", "microseism"),
                                   ((0.5, 5.0), "#c1440e", "tremor")):
            for fb in (f1, f2):
                ax.axvline(1 / fb, color=col, lw=1.0, ls=":", alpha=0.9)
            ax.text(1 / np.sqrt(f1 * f2), -96, lab, color=col, fontsize=7.5,
                    ha="center", va="top")
    ax.set_xlabel("period [s]")
    return m


def ppsd_hours(d):
    """Hours of data behind a PPSD entry (n is segments, not hours)."""
    return d["n"] * 900.0 / 3600.0


def station_color(sta):
    order = list(C.STATIONS)
    cmap = plt.get_cmap("viridis")
    return cmap(0.1 + 0.75 * order.index(sta) / max(len(order) - 1, 1))
