"""Single-station ambient-noise autocorrelation.

Why autocorrelation rather than cross-correlation, for a one-hour practical:
it needs ONE station. No pairs, no inter-station geometry, no arguing about
why one branch of the correlation is bigger than the other. The physics is
the same -- Claerbout (1968) for the 1-D case, Wapenaar (2004) in general --
and so is the monitoring workflow.

Sign convention, which trips everyone up once
---------------------------------------------
Claerbout (1968): the autocorrelation of the transmission response recorded
at the free surface equals MINUS the reflection response at non-zero lag. So
an impedance INCREASE at depth appears as a NEGATIVE peak (a trough) in the
autocorrelation, at the two-way travel time of the reflector.

Estimators
----------
CCGN : classical cross-correlation, geometrically normalised, optionally
       after spectral whitening. Fast (FFT), and what most monitoring
       studies use.
PCC2 : phase cross-correlation of Schimmel (1999) with nu = 2, which reduces
       to the correlation of unit phasors and is therefore also an FFT.
       Amplitude-unbiased: a single loud transient cannot dominate it.

The nu = 1 PCC of Romero & Schimmel (2018) is a direct O(N*M) sum and needs
numba to be tolerable, so it is deliberately left out of the teaching stack.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import butter, hilbert, sosfiltfilt


# -------------------------------------------------------- preprocessing
def bandpass(x, fs, f1, f2, order=4):
    """Zero-phase Butterworth band-pass."""
    nyq = 0.5 * fs
    sos = butter(order, [f1 / nyq, min(f2, 0.98 * nyq) / nyq],
                 btype="band", output="sos")
    return sosfiltfilt(sos, np.asarray(x, dtype=np.float64))


def onebit(x):
    """Keep only the sign. The bluntest possible defence against transients."""
    return np.sign(np.asarray(x, dtype=np.float64))


def ram_norm(x, nwin):
    """Running-absolute-mean normalisation (Bensen et al. 2007).

    Gentler than one-bit: divides by a smoothed envelope, so it suppresses
    earthquakes and tremor bursts without destroying relative amplitude
    everywhere else.
    """
    x = np.asarray(x, dtype=np.float64)
    ker = np.ones(2 * nwin + 1) / (2 * nwin + 1)
    w = np.convolve(np.abs(x), ker, mode="same")
    w[w == 0] = 1.0
    return x / w


def spectral_whiten(x, fs, band, water_level=1e-3, taper_frac=0.25,
                    smooth_hz=None):
    """Whiten the amplitude spectrum inside `band`.

    Two modes, and the difference matters enormously for a SINGLE-station
    autocorrelation:

    smooth_hz=None  -- full whitening: |X| is flattened to 1. The
        autocorrelation then tends to a band-limited delta, because the
        reverberations of the structure ARE spectral peaks and troughs (a
        comb) and flattening the spectrum erases them. Good for monitoring
        the coda shape, destructive if you want to image a reflector.
        tests/test_acorr.py measures this: it moves a 4.0 s trough to 2.05 s.

    smooth_hz=w     -- running-mean whitening: |X| is divided by its own
        average over a window of width w Hz. Choose w WIDER than the
        spectral features you want to keep (the comb spacing is 1/TWT) and
        the structure survives while the overall colour of the noise is
        removed. This is what ambient-noise imaging actually does.
    """
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    nfft = 1 << int(np.ceil(np.log2(n)))
    X = np.fft.rfft(x, nfft)
    f = np.fft.rfftfreq(nfft, d=1.0 / fs)
    f1, f2 = band
    amp = np.abs(X)
    if smooth_hz is None:
        den = amp
    else:
        k = max(3, int(round(smooth_hz / (f[1] - f[0]))) | 1)
        den = np.convolve(amp, np.ones(k) / k, mode="same")
    wl = water_level * amp.max() if amp.max() > 0 else 1.0
    Xw = X / np.maximum(den, wl)
    w = np.zeros_like(f)
    dt1 = taper_frac * (f2 - f1)
    w[(f >= f1) & (f <= f2)] = 1.0
    lo = (f > max(f1 - dt1, 0.0)) & (f < f1)
    if dt1 > 0:
        w[lo] = 0.5 * (1 - np.cos(np.pi * (f[lo] - (f1 - dt1)) / dt1))
        hi = (f > f2) & (f < f2 + dt1)
        w[hi] = 0.5 * (1 + np.cos(np.pi * (f[hi] - f2) / dt1))
    return np.fft.irfft(Xw * w, nfft)[:n]


# ----------------------------------------------------------- estimators
def ccgn_autocorr(x, maxlag, whiten=False, fs=None, band=None,
                  water_level=1e-3, smooth_hz=0.5):
    """Classical autocorrelation for lags 0..maxlag samples, c[0] = 1.

    `smooth_hz` is passed to spectral_whiten; the default of 0.5 Hz is a
    running-mean whitening that preserves structure with a two-way time
    longer than ~2 s. Pass smooth_hz=None for full whitening and watch what
    happens to the reflector.
    """
    x = np.asarray(x, dtype=np.float64)
    if maxlag >= x.size:
        raise ValueError("maxlag must be < len(x)")
    if whiten:
        if fs is None or band is None:
            raise ValueError("whiten=True needs fs and band")
        x = spectral_whiten(x, fs, band, water_level=water_level,
                            smooth_hz=smooth_hz)
    nfft = 1 << int(np.ceil(np.log2(2 * x.size)))
    X = np.fft.rfft(x, nfft)
    acf = np.fft.irfft(np.abs(X) ** 2, nfft)[:maxlag + 1]
    return acf / acf[0] if acf[0] else np.zeros(maxlag + 1)


def pcc2_autocorr(x, maxlag, theta=None):
    """Phase autocorrelation with nu = 2, by FFT. Amplitude-unbiased."""
    x = np.asarray(x, dtype=np.float64)
    if theta is None:
        theta = np.angle(hilbert(x))
    z = np.exp(1j * theta)
    n = z.size
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    Z = np.fft.fft(z, nfft)
    acf = np.real(np.fft.ifft(Z * np.conj(Z))[:maxlag + 1])
    acf = acf / (n - np.arange(maxlag + 1))      # unbiased overlap count
    return acf / acf[0] if acf[0] else np.zeros(maxlag + 1)


# -------------------------------------------------------------- stacking
def pws(traces, nu=2.0):
    """Time-domain phase-weighted stack (Schimmel & Paulssen 1997).

    weight(t) = | mean_j exp(i phi_j(t)) |^nu, where phi_j is the
    instantaneous phase of trace j. Incoherent arrivals get a small weight.
    Cheaper than the time-frequency version and enough to make the point.
    """
    traces = np.atleast_2d(np.asarray(traces, dtype=np.float64))
    lin = traces.mean(axis=0)
    ph = np.exp(1j * np.angle(hilbert(traces, axis=-1)))
    w = np.abs(ph.mean(axis=0)) ** nu
    return lin * w, lin


def daily_autocorr(tr, band, maxlag_s, window_s=1800.0, overlap=0.5,
                   estimator="ccgn", norm="ram", whiten=True, stack="linear",
                   hours=(0, 6)):
    """One day of data -> one autocorrelation function.

    The day is cut into `window_s` windows, each window is normalised and
    autocorrelated, and the windows are stacked. Stacking short windows
    rather than correlating the whole day is what makes the result robust:
    one noisy window can be down-weighted instead of contaminating
    everything.

    **`hours` is not a detail, it is the measurement.** Every day must be
    built from the SAME hours of the clock, or you are not comparing like
    with like: the ambient noise field has a strong daily cycle, so an
    autocorrelation stacked over 24 hours and one stacked over six night
    hours sample different source distributions and are simply different
    waveforms. Mixing the two drops the day-to-day correlation coefficient
    from ~0.9 to ~0.4, which looks exactly like a loss of data quality and
    is in fact a self-inflicted wound. Pass hours=None only if every day in
    your archive has identical coverage.

    Returns (lags_in_seconds, acf, n_windows_used).
    """
    if hours is not None:
        t0 = tr.stats.starttime
        # A day's trace can begin a few seconds BEFORE midnight, because the
        # first record of the archive straddles the boundary. Truncating
        # starttime to its own date would then select the previous day and
        # return almost nothing, so take the intended day from stats.day
        # when io.load_day has set it, and otherwise round to the nearest
        # midnight rather than flooring.
        day = getattr(tr.stats, "day", None)
        if day:
            day0 = type(t0)(f"{day}T00:00:00")
        else:
            day0 = type(t0)((t0 + 3600.0).strftime("%Y-%m-%dT00:00:00"))
        tr = tr.slice(day0 + hours[0] * 3600.0, day0 + hours[1] * 3600.0)
        if tr.stats.npts == 0:
            return None, None, 0
    fs = tr.stats.sampling_rate
    n = int(window_s * fs)
    maxlag = int(maxlag_s * fs)
    if n <= 2 * maxlag or tr.stats.npts < n:
        return None, None, 0
    step = max(int(n * (1 - overlap)), 1)
    acc = []
    for s in range(0, tr.stats.npts - n + 1, step):
        seg = tr.data[s:s + n].astype(np.float64)
        if not np.any(seg) or not np.all(np.isfinite(seg)):
            continue
        # a window that is mostly zero-fill is not data
        if np.count_nonzero(seg) < 0.8 * seg.size:
            continue
        seg = seg - seg.mean()
        seg = bandpass(seg, fs, band[0], band[1])
        if norm == "onebit":
            seg = onebit(seg)
        elif norm == "ram":
            seg = ram_norm(seg, int(fs / band[0]))
        if estimator == "pcc2":
            c = pcc2_autocorr(seg, maxlag)
        else:
            c = ccgn_autocorr(seg, maxlag, whiten=whiten, fs=fs, band=band)
        if np.all(np.isfinite(c)):
            acc.append(c)
    if not acc:
        return None, None, 0
    acc = np.vstack(acc)
    out = pws(acc)[0] if stack == "pws" else acc.mean(axis=0)
    return np.arange(maxlag + 1) / fs, out, len(acc)


def twt_to_depth(twt, v):
    """Zero-offset two-way time to depth for a constant velocity."""
    return 0.5 * v * twt
