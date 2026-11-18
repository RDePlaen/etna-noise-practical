"""Checks for the autocorrelation stack. Run: python -m pytest tests -q"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from etnalib import acorr


def _reverberation_series(r_base, twt_samp, n, n_mult=60):
    """Free-surface multiples for one layer over a half-space.

    Impulse from below, free surface reflection coefficient -1, base of
    layer reflection coefficient r_base (positive = impedance increase
    downwards). The surface record is a series of arrivals spaced by the
    two-way time, with alternating-sign amplitudes (-r)^n.
    """
    g = np.zeros(n)
    for k in range(n_mult):
        i = k * twt_samp
        if i >= n:
            break
        g[i] += (-r_base) ** k
    return g


def test_claerbout_sign():
    """An impedance INCREASE must appear as a NEGATIVE peak at the two-way time."""
    fs, n = 20.0, 2 ** 16
    twt_s, r = 4.0, 0.4                    # 4 s two-way time, r > 0
    twt_samp = int(twt_s * fs)
    rng = np.random.default_rng(1)
    src = rng.standard_normal(n)
    g = _reverberation_series(r, twt_samp, n)
    rec = np.convolve(src, g, mode="same")
    c = acorr.ccgn_autocorr(rec, maxlag=int(12 * fs), whiten=False)
    lags = np.arange(c.size) / fs
    i = int(np.argmin(np.abs(lags - twt_s)))
    win = c[i - 3:i + 4]
    assert c[0] == 1.0
    assert win.min() < -0.1, f"expected a trough at {twt_s} s, got {win}"
    assert abs(win.min()) > abs(c[i + 10]), "trough is not above the background"
    return float(win.min())


def test_impedance_decrease_flips_sign():
    """r < 0 (impedance decrease) must give a POSITIVE peak."""
    fs, n = 20.0, 2 ** 16
    twt_s, r = 4.0, -0.4
    twt_samp = int(twt_s * fs)
    rng = np.random.default_rng(2)
    rec = np.convolve(rng.standard_normal(n),
                      _reverberation_series(r, twt_samp, n), mode="same")
    c = acorr.ccgn_autocorr(rec, maxlag=int(12 * fs), whiten=False)
    lags = np.arange(c.size) / fs
    i = int(np.argmin(np.abs(lags - twt_s)))
    assert c[i - 3:i + 4].max() > 0.1


def test_pcc2_is_amplitude_unbiased():
    """One loud transient must move PCC2 much less than it moves CCGN."""
    fs, n = 20.0, 2 ** 15
    rng = np.random.default_rng(3)
    x = acorr.bandpass(rng.standard_normal(n), fs, 1.0, 3.0)
    y = x.copy()
    y[n // 2:n // 2 + 200] += 50 * np.hanning(200) * x[n // 2:n // 2 + 200]
    ml = int(10 * fs)
    dc = np.abs(acorr.ccgn_autocorr(y, ml) - acorr.ccgn_autocorr(x, ml)).max()
    dp = np.abs(acorr.pcc2_autocorr(y, ml) - acorr.pcc2_autocorr(x, ml)).max()
    assert dp < dc, f"PCC2 ({dp:.4f}) should be less perturbed than CCGN ({dc:.4f})"
    return dc, dp


def test_full_whitening_destroys_the_reflector():
    """Full whitening flattens away the comb that encodes the multiples.

    This is not a bug, it is what whitening means, and it is the reason the
    default is a running-mean whitening instead.
    """
    fs, n = 20.0, 2 ** 15
    twt_samp = int(4.0 * fs)
    rng = np.random.default_rng(4)
    rec = np.convolve(rng.standard_normal(n),
                      _reverberation_series(0.4, twt_samp, n), mode="same")
    ml = int(12 * fs)
    a = acorr.ccgn_autocorr(rec, ml, whiten=False)
    full = acorr.ccgn_autocorr(rec, ml, whiten=True, fs=fs, band=(0.5, 5.0),
                               smooth_hz=None)
    la = np.argmin(a[int(2 * fs):]) + int(2 * fs)
    lf = np.argmin(full[int(2 * fs):]) + int(2 * fs)
    assert abs(la - lf) > 5, "full whitening was expected to move the trough"
    return la / fs, lf / fs


def test_smoothed_whitening_preserves_lag():
    """Running-mean whitening must NOT move an arrival in time."""
    fs, n = 20.0, 2 ** 15
    twt_samp = int(4.0 * fs)
    rng = np.random.default_rng(4)
    rec = np.convolve(rng.standard_normal(n),
                      _reverberation_series(0.4, twt_samp, n), mode="same")
    ml = int(12 * fs)
    a = acorr.ccgn_autocorr(rec, ml, whiten=False)
    b = acorr.ccgn_autocorr(rec, ml, whiten=True, fs=fs, band=(0.5, 5.0),
                            smooth_hz=0.5)
    la = np.argmin(a[int(2 * fs):]) + int(2 * fs)
    lb = np.argmin(b[int(2 * fs):]) + int(2 * fs)
    assert abs(la - lb) <= 1, f"smoothed whitening moved the trough: {la} vs {lb}"
    return la / fs, lb / fs


if __name__ == "__main__":
    v = test_claerbout_sign()
    print(f"test_claerbout_sign          OK  (trough amplitude {v:+.3f})")
    test_impedance_decrease_flips_sign()
    print("test_impedance_decrease      OK")
    dc, dp = test_pcc2_is_amplitude_unbiased()
    print(f"test_pcc2_amplitude_unbiased OK  (CCGN moved {dc:.4f}, PCC2 {dp:.4f})")
    ra, rf = test_full_whitening_destroys_the_reflector()
    print(f"test_full_whitening_destroys OK  (trough {ra:.2f} s -> {rf:.2f} s)")
    la, lb = test_smoothed_whitening_preserves_lag()
    print(f"test_smoothed_whitening      OK  (trough {la:.2f} s -> {lb:.2f} s)")
