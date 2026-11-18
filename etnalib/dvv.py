"""Measuring a relative velocity change from autocorrelation functions.

The stretching method (Sens-Schoenfelder & Wegler 2006, Lobkis & Weaver 2003)
assumes the medium's velocity changes by the same relative amount everywhere
that the coda samples. A uniform velocity change dv/v stretches or squeezes
the whole coda in time:

    if  v -> v (1 + dv/v)   then   t -> t / (1 + dv/v)

so the current waveform is a time-stretched copy of the reference. We search
over trial stretches, correlate, and take the one that fits best. The peak
correlation coefficient is reported alongside: a dv/v with a CC of 0.3 is not
a measurement, it is a number.

Sign convention: dv/v > 0 means the medium got FASTER, so arrivals move to
EARLIER times. Checked in test_dvv() below against a synthetic.
"""
from __future__ import annotations

import numpy as np


def _corr(a, b):
    a = a - a.mean()
    b = b - b.mean()
    d = np.sqrt(np.sum(a * a) * np.sum(b * b))
    return float(np.sum(a * b) / d) if d > 0 else 0.0


def stretching(cur, ref, lags, lag_win, eps=0.05, neps=201):
    """Best-fitting uniform stretch between `cur` and `ref`.

    Parameters
    ----------
    cur, ref : the current and reference autocorrelations, same time base
    lags     : lag times in seconds
    lag_win  : (t1, t2), the coda window used for the measurement. Never
               include the zero-lag peak: it is identical in every trace by
               construction and will force the answer towards zero.
    eps      : half-width of the search, as a fraction (0.05 = +/- 5 %)
    neps     : number of trial values

    Returns
    -------
    dict with dvv (in per cent), cc, and the full (eps_grid, cc_curve) so the
    student can look at the shape of the misfit rather than trusting a number.
    """
    cur = np.asarray(cur, float)
    ref = np.asarray(ref, float)
    lags = np.asarray(lags, float)
    m = (lags >= lag_win[0]) & (lags <= lag_win[1])
    if m.sum() < 16:
        return dict(dvv=np.nan, cc=np.nan, eps=np.array([]), ccs=np.array([]))

    grid = np.linspace(-eps, eps, neps)
    ccs = np.empty(neps)
    for i, e in enumerate(grid):
        # a medium faster by e has its coda at times t/(1+e); to compare with
        # the reference we resample the reference onto the stretched axis
        ref_s = np.interp(lags[m] * (1.0 + e), lags, ref,
                          left=np.nan, right=np.nan)
        ok = np.isfinite(ref_s)
        ccs[i] = _corr(cur[m][ok], ref_s[ok]) if ok.sum() > 16 else np.nan

    if not np.any(np.isfinite(ccs)):
        return dict(dvv=np.nan, cc=np.nan, eps=grid, ccs=ccs)
    k = int(np.nanargmax(ccs))
    # parabolic refinement on the three points around the maximum
    if 0 < k < neps - 1 and np.all(np.isfinite(ccs[k - 1:k + 2])):
        y0, y1, y2 = ccs[k - 1], ccs[k], ccs[k + 1]
        den = y0 - 2 * y1 + y2
        shift = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        best = grid[k] + shift * (grid[1] - grid[0])
    else:
        best = grid[k]
    return dict(dvv=100.0 * best, cc=float(ccs[k]), eps=grid, ccs=ccs)


def reference_stack(acfs, days, ref_days=None):
    """The reference against which every day is measured.

    `ref_days` lets you choose the reference period explicitly, which is the
    single most consequential choice in the whole measurement.

    **The default -- the mean of every day -- is the default for a reason.**
    A reference built from one end of the window guarantees that the other
    end is least similar to it, so the correlation coefficient decays with
    distance from the reference and the measurement quietly dies. Measured
    on this dataset (EMSG, 1-3 Hz, 3-15 s coda, 72 days): a reference taken
    from the first fortnight gives a median cc of 0.43 and 15 usable days,
    while the grand mean gives 0.79 and 44. It also removes a circularity --
    days inside the reference window correlate with it partly because they
    built it.

    The cost is that dv/v is then relative to an average state that has no
    physical meaning, so only CHANGES are interpretable, never the absolute
    level. That is true of almost every published dv/v curve.
    """
    acfs = np.asarray(acfs, float)
    if ref_days is None:
        sel = np.ones(len(days), bool)
    else:
        d0, d1 = ref_days
        sel = np.array([d0 <= d <= d1 for d in days])
    if sel.sum() == 0:
        raise ValueError("reference window contains no days")
    return np.nanmean(acfs[sel], axis=0), sel


def dvv_series(acfs, days, lags, lag_win, ref_days=None, eps=0.05, neps=201):
    """dv/v for every day against a common reference."""
    ref, sel = reference_stack(acfs, days, ref_days)
    out = []
    for i, d in enumerate(days):
        r = stretching(acfs[i], ref, lags, lag_win, eps=eps, neps=neps)
        out.append((d, r["dvv"], r["cc"]))
    return out, ref, sel


def test_dvv(verbose=True):
    """Plant a known velocity change in a synthetic coda and recover it.

    Nothing in this module is trusted because it ran.
    """
    fs, maxlag = 20.0, 20.0
    t = np.arange(0, maxlag, 1 / fs)
    rng = np.random.default_rng(0)
    env = np.exp(-t / 6.0)
    ref = env * np.convolve(rng.standard_normal(t.size),
                            np.hanning(9), mode="same")
    errs = []
    for true in (-2.0, -0.5, 0.0, 0.5, 2.0):
        e = true / 100.0
        cur = np.interp(t, t / (1.0 + e), ref, left=0, right=0)
        r = stretching(cur, ref, t, (3.0, 15.0), eps=0.05, neps=201)
        errs.append(abs(r["dvv"] - true))
        if verbose:
            print(f"  planted {true:+5.2f} %   recovered {r['dvv']:+6.3f} %"
                  f"   cc={r['cc']:.4f}")
    assert max(errs) < 0.05, f"stretching is biased: worst {max(errs):.3f} %"
    if verbose:
        print(f"  worst error {max(errs):.4f} % -- OK")
    return max(errs)


if __name__ == "__main__":
    test_dvv()
