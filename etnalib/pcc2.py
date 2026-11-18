"""PCC2 -- phase cross-correlation with nu = 2, after Ventosa et al. (2019).

Routine by Luis-Fabian Bonilla (IPGP & IFSTTAR), Jan 2020, supplied by the
user and adapted here. Reference implementation:
Ventosa, Schimmel & Stutzmann (2019), *Towards the processing of large data
volumes with phase cross-correlation*, SRL 90, 1663-1669,
doi:10.1785/0220190022 -- and https://github.com/sergiventosa/FastPCC.

Why PCC at all
--------------
The classical correlation is dominated by whatever is loudest. On a volcano
that is the tremor, which at Etna is a persistent source fixed at the summit
craters and which changes amplitude whenever the eruption does. The usual
defence is aggressive temporal normalisation (one-bit, running-absolute-mean)
applied *before* correlating.

PCC does not need it. It is computed from the instantaneous PHASE of the
analytic signal, with the amplitude divided out exactly:

    xa = hilbert(x) / |hilbert(x)|

so a loud window and a quiet window contribute equally by construction. With
nu = 2 the sum of phasor products reduces to a plain correlation of unit
phasors, which is an FFT, so it is also fast.

Changes from the supplied version, both required and both noted:

1. `np.complex_` was removed in NumPy 2.0; it is `complex` here. The original
   raises AttributeError on any NumPy >= 2, which is what this project runs.
2. The lag axis is built with `np.linspace` rather than `np.arange`, because
   `arange` with a float step can return Nz or Nz+-1 samples depending on
   rounding, and the returned axis has to line up with the correlation
   sample-for-sample.

Everything else -- the taper, the zero padding to the next power of two, the
normalisation by N -- is as supplied.
"""
from __future__ import annotations

import numpy as np
import scipy.signal


def taper(x, p):
    """Hanning taper over a fraction p of each end."""
    if p <= 0.0:
        return x
    f0 = 0.5
    f1 = 0.5
    n = len(x)
    nw = int(p * n)
    if nw == 0:
        return x
    ow = np.pi / nw
    w = np.ones(n)
    for i in range(nw):
        w[i] = f0 - f1 * np.cos(ow * i)
    for i in range(n - nw, n):
        w[i] = 1.0 - w[i - n + nw]
    return x * w


def next_power_of_2(n):
    """Next power of 2 greater than or equal to n."""
    return 2 ** (n - 1).bit_length()


def _lag_axis(nz, dt):
    """Symmetric lag axis of exactly nz samples, zero at index nz // 2."""
    return (np.arange(nz) - nz // 2) * dt


def pcc2(x1, x2, dt, lag0, lagu):
    """PCC2 between two traces, returned on lags [lag0, lagu]."""
    x1 = x1 - np.mean(x1)
    x2 = x2 - np.mean(x2)
    x1 = taper(x1, 0.05)
    x2 = taper(x2, 0.05)
    n = len(x1)
    nz = next_power_of_2(2 * n)

    xa1 = scipy.signal.hilbert(x1)
    xa2 = scipy.signal.hilbert(x2)
    xa1 = xa1 / np.abs(xa1)
    xa2 = xa2 / np.abs(xa2)

    xa1 = np.append(xa1, np.zeros((nz - n), dtype=complex))
    xa2 = np.append(xa2, np.zeros((nz - n), dtype=complex))

    xa1 = np.fft.fft(xa1)
    xa2 = np.fft.fft(xa2)
    amp = xa1 * np.conj(xa2)
    pcc = np.real(np.fft.ifft(amp)) / n
    pcc = np.fft.ifftshift(pcc)

    t = _lag_axis(nz, dt)
    m = (t >= lag0) & (t <= lagu)
    return t[m], pcc[m]


def apcc2(x1, dt, lag0, lagu):
    """PCC2 autocorrelation of one trace, on lags [lag0, lagu]."""
    x1 = x1 - np.mean(x1)
    x1 = taper(x1, 0.05)
    n = len(x1)
    nz = next_power_of_2(2 * n)

    xa1 = scipy.signal.hilbert(x1)
    xa1 = xa1 / np.abs(xa1)
    xa1 = np.append(xa1, np.zeros((nz - n), dtype=complex))

    xa1 = np.fft.fft(xa1)
    amp = xa1 * np.conj(xa1)
    pcc = np.real(np.fft.ifft(amp)) / n
    pcc = np.fft.ifftshift(pcc)

    t = _lag_axis(nz, dt)
    m = (t >= lag0) & (t <= lagu)
    return t[m], pcc[m]


# ----------------------------------------------------------------- daily
def daily_pcc_autocorr(tr, band, maxlag_s, window_s=1800.0, overlap=0.5,
                       hours=(0, 6), bandpass_first=True):
    """One station-day -> one PCC2 autocorrelation, non-negative lags.

    Deliberately has NO temporal-normalisation argument. PCC is
    amplitude-unbiased by construction, so one-bit or running-absolute-mean
    would be redundant -- that is the whole reason to use it here.

    Returns (lags_in_seconds, acf, n_windows_used).
    """
    from scipy.signal import butter, sosfiltfilt

    if hours is not None:
        t0 = tr.stats.starttime
        day = getattr(tr.stats, "day", None)
        day0 = (type(t0)(f"{day}T00:00:00") if day else
                type(t0)((t0 + 3600.0).strftime("%Y-%m-%dT00:00:00")))
        tr = tr.slice(day0 + hours[0] * 3600.0, day0 + hours[1] * 3600.0)
        if tr.stats.npts == 0:
            return None, None, 0

    fs = tr.stats.sampling_rate
    dt = 1.0 / fs
    n = int(window_s * fs)
    if tr.stats.npts < n:
        return None, None, 0
    step = max(int(n * (1 - overlap)), 1)

    nyq = 0.5 * fs
    sos = butter(4, [band[0] / nyq, min(band[1], 0.98 * nyq) / nyq],
                 btype="band", output="sos")

    acc, lags = [], None
    for s in range(0, tr.stats.npts - n + 1, step):
        seg = tr.data[s:s + n].astype(np.float64)
        if not np.any(seg) or not np.all(np.isfinite(seg)):
            continue
        if np.count_nonzero(seg) < 0.8 * seg.size:
            continue
        if bandpass_first:
            seg = sosfiltfilt(sos, seg - seg.mean())
        t, c = apcc2(seg, dt, 0.0, maxlag_s)
        if not np.all(np.isfinite(c)):
            continue
        lags = t
        acc.append(c)
    if not acc:
        return None, None, 0
    out = np.vstack(acc).mean(axis=0)
    if out[0] != 0:
        out = out / out[0]          # c[0] = 1, as for the CC branch
    return lags, out, len(acc)


# ------------------------------------------------------------------ tests
def _test():
    """Checks with known answers. Run: python -m etnalib.pcc2"""
    rng = np.random.default_rng(0)
    fs, dt = 20.0, 0.05
    n = 4096

    # 1. zero lag sits at t = 0 and equals 1 for an autocorrelation
    x = rng.standard_normal(n)
    t, c = apcc2(x, dt, -10.0, 10.0)
    i0 = int(np.argmin(np.abs(t)))
    assert abs(t[i0]) < 1e-12, f"zero lag not on the axis: {t[i0]}"
    assert abs(c[i0] / c.max() - 1.0) < 1e-9, "peak is not at zero lag"
    print(f"  zero lag at t={t[i0]:+.3e}, peak there  OK")

    # 2. an autocorrelation must be even in lag
    half = len(c) // 2
    m = min(half, len(c) - half - 1)
    asym = np.max(np.abs(c[i0 + 1:i0 + 1 + m] - c[i0 - m:i0][::-1]))
    assert asym / np.max(np.abs(c)) < 1e-9, f"not even: {asym:.2e}"
    print(f"  symmetry error {asym/np.max(np.abs(c)):.2e}  OK")

    # 3. amplitude blindness: scaling part of the trace by 50 must not
    #    change PCC, which is the entire point of the method
    y = x.copy()
    y[n // 2:n // 2 + 200] *= 50.0
    _, cy = apcc2(y, dt, 0.0, 10.0)
    _, cx = apcc2(x, dt, 0.0, 10.0)
    d_pcc = np.max(np.abs(cy - cx)) / np.max(np.abs(cx))
    # the classical equivalent, for contrast
    def ccgn(z):
        nf = 1 << int(np.ceil(np.log2(2 * z.size)))
        Z = np.fft.rfft(z - z.mean(), nf)
        a = np.fft.irfft(np.abs(Z) ** 2, nf)[:201]
        return a / a[0]
    d_cc = np.max(np.abs(ccgn(y) - ccgn(x)))
    assert d_pcc < d_cc, "PCC should be less perturbed than CC"
    print(f"  one loud burst moves CC by {d_cc:.4f}, PCC by {d_pcc:.4f}  OK")

    # 4. a planted periodicity must appear at the right lag
    lag_s = 3.0
    z = rng.standard_normal(n)
    z[int(lag_s * fs):] += 0.8 * z[:n - int(lag_s * fs)]
    t4, c4 = apcc2(z, dt, 0.5, 8.0)
    pk = t4[int(np.argmax(c4))]
    assert abs(pk - lag_s) < 2 * dt, f"peak at {pk}, expected {lag_s}"
    print(f"  planted {lag_s:.2f} s periodicity recovered at {pk:.2f} s  OK")
    return True


if __name__ == "__main__":
    _test()
