"""Moving-window cross-spectral analysis (MWCS) for dv/v.

This is the method used in De Plaen et al. (2019), *Temporal changes of
seismic velocity caused by volcanic activity at Mt. Etna revealed by the
autocorrelation of ambient seismic noise*, Front. Earth Sci. 6:251,
doi:10.3389/feart.2018.00251 -- so it is implemented here to match, rather
than substituting the stretching in `etnalib.dvv`.

How it differs from stretching, and why it matters
--------------------------------------------------
Stretching assumes a single uniform relative velocity change and searches
for the best global time warp. MWCS makes no such assumption up front: it
measures a time shift `dt` independently in each short window of the coda,
by fitting the slope of the cross-spectral PHASE against frequency, and only
then regresses dt against lag time t. For a homogeneous velocity change the
regression is a straight line through the origin with slope dt/t = -dv/v.

The advantages are that you see whether that line is actually straight --
if it is not, the change is not homogeneous -- and that each point carries
its own error bar and coherence, which is where the quality control comes
from.

Conventions follow Clarke et al. (2011) and the MSNoise implementation:
dv/v = -dt/t, so a POSITIVE dv/v (faster medium) means arrivals move to
EARLIER times (negative dt).
"""
from __future__ import annotations

import numpy as np

try:
    from scipy.signal.windows import hann as _hann
except Exception:                      # pragma: no cover
    def _hann(n, sym=True):
        return np.hanning(n)


def _smooth(x, half):
    """Running mean over 2*half+1 samples, edges handled by reflection."""
    if half <= 0:
        return x
    k = np.ones(2 * half + 1) / (2 * half + 1)
    return np.convolve(np.pad(x, half, mode="reflect"), k, mode="valid")


def mwcs(cur, ref, fs, fmin, fmax, tmin, tmax,
         window_length=12.0, step=4.0, smoothing_half=5):
    """Per-window time shifts between `cur` and `ref`.

    Parameters
    ----------
    cur, ref : 1-D arrays on a common lag axis starting at lag 0
    fs       : sampling rate, Hz
    fmin,fmax: the band over which the phase slope is fitted
    tmin,tmax: the coda window, seconds of lag
    window_length, step : moving window geometry, seconds

    Returns a dict of arrays: t (window centre), dt, err, coh.
    """
    cur = np.asarray(cur, float)
    ref = np.asarray(ref, float)
    n = min(cur.size, ref.size)
    lags = np.arange(n) / fs

    wl = int(round(window_length * fs))
    st = int(round(step * fs))
    if wl < 16:
        raise ValueError("window_length too short for this sampling rate")
    nfft = 1 << int(np.ceil(np.log2(wl * 2)))
    freqs = np.fft.rfftfreq(nfft, d=1.0 / fs)
    band = (freqs >= fmin) & (freqs <= fmax)
    if band.sum() < 4:
        raise ValueError("fewer than 4 frequency samples in the band")

    taper = _hann(wl)
    out_t, out_dt, out_err, out_coh = [], [], [], []

    starts = range(0, n - wl + 1, st)
    for s in starts:
        tc = lags[s] + window_length / 2.0
        if tc < tmin or tc > tmax:
            continue
        a = cur[s:s + wl] * taper
        b = ref[s:s + wl] * taper
        if not np.any(a) or not np.any(b):
            continue
        A = np.fft.rfft(a, nfft)
        B = np.fft.rfft(b, nfft)

        # coherence from smoothed spectra
        saa = _smooth(np.abs(A) ** 2, smoothing_half)
        sbb = _smooth(np.abs(B) ** 2, smoothing_half)
        sab = _smooth(A * np.conj(B), smoothing_half)
        den = np.sqrt(saa * sbb)
        with np.errstate(divide="ignore", invalid="ignore"):
            coh = np.abs(sab) / den
        coh[~np.isfinite(coh)] = 0.0
        coh = np.clip(coh, 0.0, 1.0)

        phi = np.angle(sab)
        phi = np.unwrap(phi)

        f = freqs[band]
        p = phi[band]
        w = coh[band] ** 2
        x = 2 * np.pi * f
        sw = np.sum(w * x * x)
        if sw <= 0:
            continue
        # SIGN. With S_ab = A conj(B) and A = cur, B = ref, a current trace
        # DELAYED by tau relative to the reference gives
        #     A(f) = B(f) exp(-2 i pi f tau)   ->   phase(S_ab) = -2 pi f tau
        # so the least-squares slope of phi against 2*pi*f is MINUS the
        # delay. Returning the slope directly would report every velocity
        # change with the wrong sign -- and since the magnitude is right,
        # nothing downstream would look broken. Hence the explicit minus,
        # and the planted-change test that pins it down.
        dt = -float(np.sum(w * x * p) / sw)
        resid = p + dt * x
        dof = max(band.sum() - 1, 1)
        sigma2 = float(np.sum(w * resid ** 2) / (np.sum(w) * dof))
        err = float(np.sqrt(sigma2 * band.sum() / sw)) if sw > 0 else np.nan

        out_t.append(tc)
        out_dt.append(dt)
        out_err.append(err)
        out_coh.append(float(np.mean(coh[band])))

    return dict(t=np.array(out_t), dt=np.array(out_dt),
                err=np.array(out_err), coh=np.array(out_coh))


def dtt(m, max_err=0.1, min_coh=0.6, min_points=5):
    """Regress dt on t and return dv/v in per cent.

    Quality control follows De Plaen et al. (2019): a window is dropped if
    its dt error exceeds `max_err` seconds or its coherence is below
    `min_coh`.

    Returns (dvv_percent, error_percent, n_used, coh_median).
    """
    t, dt_, err, coh = m["t"], m["dt"], m["err"], m["coh"]
    if t.size == 0:
        return np.nan, np.nan, 0, np.nan
    keep = (err <= max_err) & (coh >= min_coh) & np.isfinite(dt_)
    if keep.sum() < min_points:
        return np.nan, np.nan, int(keep.sum()), float(np.nanmedian(coh))
    tt, dd = t[keep], dt_[keep]
    w = 1.0 / np.maximum(err[keep], 1e-6) ** 2
    # slope through the origin: dt = s * t, dv/v = -s
    s = np.sum(w * tt * dd) / np.sum(w * tt * tt)
    resid = dd - s * tt
    dof = max(keep.sum() - 1, 1)
    se = np.sqrt(np.sum(w * resid ** 2) / (np.sum(w) * dof) /
                 np.sum(w * tt * tt) * np.sum(w))
    return -100.0 * s, 100.0 * se, int(keep.sum()), float(np.median(coh[keep]))


def dvv_mwcs(cur, ref, fs, fmin, fmax, tmin, tmax, **kw):
    """Convenience: MWCS then regression. Returns the dtt() tuple."""
    qc = dict(max_err=kw.pop("max_err", 0.1), min_coh=kw.pop("min_coh", 0.6),
              min_points=kw.pop("min_points", 5))
    return dtt(mwcs(cur, ref, fs, fmin, fmax, tmin, tmax, **kw), **qc)


# ------------------------------------------------------------------ test
def _test(verbose=True):
    """Plant a known homogeneous velocity change and recover it."""
    fs, maxlag = 20.0, 40.0
    t = np.arange(0, maxlag, 1 / fs)
    rng = np.random.default_rng(3)
    from scipy.signal import butter, sosfiltfilt
    sos = butter(4, [1.0 / (fs / 2), 2.0 / (fs / 2)], btype="band",
                 output="sos")
    ref = sosfiltfilt(sos, rng.standard_normal(t.size)) * np.exp(-t / 18.0)

    errs = []
    for true in (-0.4, -0.1, 0.0, 0.1, 0.4):
        e = true / 100.0
        cur = np.interp(t, t / (1.0 + e), ref, left=0, right=0)
        v, se, n, coh = dvv_mwcs(cur, ref, fs, 1.0, 2.0, 5.0, 35.0,
                                 window_length=12.0, step=4.0)
        errs.append(abs(v - true))
        if verbose:
            print(f"  planted {true:+5.2f} %   recovered {v:+6.3f} "
                  f"+/- {se:.3f} %   n={n:2d}  coh={coh:.3f}")
    worst = max(errs)
    assert worst < 0.05, f"MWCS is biased: worst {worst:.3f} %"
    if verbose:
        print(f"  worst error {worst:.4f} % -- OK")
    return worst


if __name__ == "__main__":
    _test()
