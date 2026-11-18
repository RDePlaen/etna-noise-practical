"""Characterising the noise field: spectrograms, PPSDs, and band amplitudes.

Three ways of looking at the same continuous record, in increasing order of
compression:

  * spectrogram          -- power as a function of time AND frequency
  * PPSD (McNamara &     -- the probability density of power at each frequency,
    Buland 2004)            with time thrown away
  * band amplitude       -- one number per band per hour, with frequency
                            thrown away

The practical uses all three, because each one makes a different thing
obvious: the spectrogram shows that a storm arrives, the PPSD shows that the
station is otherwise normal, and the band amplitude is what you would actually
monitor operationally.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import spectrogram as _spec
from scipy.signal import welch

from . import config as C
from .io import inventory, load_day


# ------------------------------------------------------------------ PPSD
def ppsd_for_day(sta, day, cha=C.CHANNEL, ppsd=None):
    """Add one station-day to an ObsPy PPSD, creating it if needed.

    The trace must NOT be response-corrected -- PPSD removes the response
    itself, and doing it twice is a classic and silent error. So this reads
    the day with remove_response=False and at full sampling rate.
    """
    from obspy.signal import PPSD
    tr = load_day(sta, day, cha, fs_out=1e9, remove_response=False,
                  min_coverage=0.3)
    if tr is None:
        return ppsd
    inv = inventory()
    if ppsd is None:
        resp = inv.get_response(f"{C.NETWORK}.{sta}..{cha}", tr.stats.starttime)
        ppsd = PPSD(tr.stats, metadata=inv,
                    ppsd_length=1800.0, overlap=0.5,
                    period_smoothing_width_octaves=0.25,
                    period_step_octaves=0.05,
                    db_bins=(-200, -50, 0.5))
    ppsd.add(tr)
    return ppsd


# ------------------------------------------------- spectrogram, long window
def long_spectrogram(sta, days, cha=C.CHANNEL, fs=20.0,
                     nperseg=2048, noverlap=1024, fmax=9.0):
    """A multi-day spectrogram, assembled day by day.

    Missing days become columns of NaN rather than being silently skipped,
    so that the time axis of the figure is real time. A spectrogram with the
    gaps closed up is a lie about when things happened.
    """
    cols, times = [], []
    for d in days:
        tr = load_day(sta, d, cha, fs_out=fs, remove_response=True,
                      min_coverage=0.3)
        t_day = np.datetime64(d)
        if tr is None:
            f = np.fft.rfftfreq(nperseg, 1.0 / fs)
            keep = f <= fmax
            cols.append(np.full((keep.sum(), 1), np.nan))
            times.append(np.array([t_day]))
            continue
        f, t, S = _spec(tr.data, fs=fs, nperseg=nperseg, noverlap=noverlap,
                        scaling="density", mode="psd")
        keep = f <= fmax
        cols.append(S[keep])
        t0 = np.datetime64(str(tr.stats.starttime)[:19])
        times.append(t0 + (t * 1e9).astype("timedelta64[ns]"))
    f = np.fft.rfftfreq(nperseg, 1.0 / fs)
    f = f[f <= fmax]
    S = np.hstack(cols)
    t = np.concatenate(times)
    order = np.argsort(t)
    return t[order], f, S[:, order]


# ------------------------------------------------------- band amplitudes
def band_amplitude(tr, band, seg_s=1800.0):
    """RMS ground velocity in a frequency band, per segment, in m/s.

    This is the quantity behind both RSAM (volcano tremor monitoring, in the
    tremor band) and the microseism amplitude used in storm studies (in the
    secondary band). Same arithmetic, different band, completely different
    literature -- worth saying out loud.
    """
    fs = tr.stats.sampling_rate
    n = int(seg_s * fs)
    if n < 16 or tr.stats.npts < n:
        return np.zeros(0), np.zeros(0)
    nseg = tr.stats.npts // n
    out, tmid = [], []
    f1, f2 = band
    for i in range(nseg):
        seg = tr.data[i * n:(i + 1) * n].astype(float)
        # A segment that straddles a data gap is part zero-fill, and its RMS
        # is biased low in proportion to how much of it is missing. That
        # shows up as a downward spike in an otherwise smooth time series,
        # which looks like a sudden quiet period and is nothing of the kind.
        # Reject rather than report.
        if np.count_nonzero(seg) < 0.9 * seg.size:
            out.append(np.nan)
        else:
            seg = seg - seg.mean()
            nf = 1 << int(np.ceil(np.log2(seg.size)))
            X = np.fft.rfft(seg, nf)
            fr = np.fft.rfftfreq(nf, 1.0 / fs)
            m = (fr >= f1) & (fr <= f2)
            # Parseval: rms in band = sqrt(2 * sum|X|^2 / (N*nf)) ... use the
            # one-sided PSD integral instead, which is unambiguous
            psd = (np.abs(X) ** 2) * 2.0 / (fs * seg.size)
            out.append(float(np.sqrt(np.trapezoid(psd[m], fr[m]))))
        tmid.append(float((i + 0.5) * seg_s))
    t0 = np.datetime64(str(tr.stats.starttime)[:19])
    t = t0 + (np.array(tmid) * 1e9).astype("timedelta64[ns]")
    return t, np.array(out)


def band_series(sta, days, bands=None, cha=C.CHANNEL, fs=20.0, seg_s=1800.0):
    """Band amplitudes for a list of days -> dict of (times, values)."""
    bands = bands or C.BANDS
    acc = {k: ([], []) for k in bands}
    for d in days:
        tr = load_day(sta, d, cha, fs_out=fs, remove_response=True,
                      min_coverage=0.3)
        if tr is None:
            continue
        for k, b in bands.items():
            if b[1] >= 0.5 * fs:
                continue
            t, a = band_amplitude(tr, b, seg_s=seg_s)
            acc[k][0].append(t)
            acc[k][1].append(a)
    out = {}
    for k, (ts, vs) in acc.items():
        if not ts:
            out[k] = (np.array([], dtype="datetime64[ns]"), np.array([]))
            continue
        t = np.concatenate(ts)
        v = np.concatenate(vs)
        o = np.argsort(t)
        out[k] = (t[o], v[o])
    return out


def mra(amplitude_ms, period_s=5.0):
    """Microseism Reduced Amplitude, in micrometres.

    The amplitude measure used in the Mediterranean storm literature
    (Cannata et al. 2022 and after): displacement amplitude scaled so that
    stations at different distances can be compared. Here we only convert
    the band-limited velocity RMS to a displacement scale, which is the part
    that matters for the shape of the time series.
    """
    return amplitude_ms / (2 * np.pi / period_s) * 1e6
