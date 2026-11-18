"""Horizontal-to-vertical spectral ratio (Nakamura 1989), SESAME conventions.

Adapted for teaching from the implementation used in this group's own work.
Conventions, unchanged:
  * horizontals combined as sqrt((E^2 + N^2) / 2)
  * Konno-Ohmachi smoothing with b = 40
  * geometric mean and geometric standard deviation over windows

The quantity of interest is f0, the fundamental resonance of a soft layer,
which for a layer of thickness h over a stiff half-space sits at Vs / (4h).

References
----------
Nakamura (1989), QR RTRI 30, 25-33
SESAME (2004), European guidelines for H/V, deliverable D23.12
Konno & Ohmachi (1998), BSSA 88, 228-241
"""
from __future__ import annotations

import numpy as np


def konno_ohmachi_matrix(freqs, fcs, b=40.0):
    """Smoothing matrix: row i smooths a spectrum onto output frequency fcs[i].

        W(f, fc) = [ sin(b log10(f/fc)) / (b log10(f/fc)) ]^4

    Built once and reused for every window and every station -- the matrix is
    the whole cost of Konno-Ohmachi smoothing and it depends only on the two
    frequency grids.
    """
    freqs = np.asarray(freqs, dtype=float)
    fcs = np.asarray(fcs, dtype=float)
    W = np.zeros((len(fcs), len(freqs)))
    with np.errstate(divide="ignore", invalid="ignore"):
        for i, fc in enumerate(fcs):
            if fc <= 0:
                continue
            x = b * np.log10(np.where(freqs > 0, freqs, np.nan) / fc)
            w = (np.sin(x) / x) ** 4
            w[~np.isfinite(w)] = 0.0
            w[np.isclose(freqs, fc)] = 1.0
            w[np.abs(np.log10(np.where(freqs > 0, freqs, np.nan) / fc)) > 0.5] = 0.0
            s = w.sum()
            if s > 0:
                W[i] = w / s
    return W


def window_spectra(data, dt, win_s, overlap=0.5, taper=0.05):
    """Split into tapered windows and return their amplitude spectra.

    Returns (spectra (nwin, nfreq), freqs, window rms)."""
    n = int(round(win_s / dt))
    if n < 16 or len(data) < n:
        return np.zeros((0, 0)), np.zeros(0), np.zeros(0)
    step = max(int(n * (1 - overlap)), 1)
    starts = range(0, len(data) - n + 1, step)
    nfft = int(2 ** np.ceil(np.log2(n)))
    freqs = np.fft.rfftfreq(nfft, d=dt)
    m = int(n * taper)
    w = np.ones(n)
    if m > 0:
        ramp = 0.5 * (1 - np.cos(np.pi * np.arange(m) / m))
        w[:m] = ramp
        w[-m:] = ramp[::-1]
    out, rms = [], []
    idx = np.arange(n)
    for s in starts:
        seg = np.asarray(data[s:s + n], dtype=float)
        if not np.all(np.isfinite(seg)):
            continue
        seg = seg - seg.mean()
        seg = seg - np.polyval(np.polyfit(idx, seg, 1), idx)
        rms.append(float(np.sqrt(np.mean(seg ** 2))))
        out.append(np.abs(np.fft.rfft(seg * w, nfft)) * dt)
    if not out:
        return np.zeros((0, len(freqs))), freqs, np.zeros(0)
    return np.vstack(out), freqs, np.asarray(rms)


def select_windows(rms_z, lo=0.4, hi=2.5):
    """Anti-triggering: keep windows whose rms is close to the median.

    SESAME uses a short-term/long-term average criterion to throw out
    transients; on an hour of continuous noise, rejecting windows whose rms
    departs from the station's own median does the same job and needs no
    extra parameters. Transients (earthquakes, machinery, footsteps) raise
    the rms; data gaps filled with zeros lower it.
    """
    if len(rms_z) == 0:
        return np.zeros(0, bool)
    pos = rms_z[rms_z > 0]
    med = np.median(pos) if pos.size else 0.0
    if med <= 0:
        return np.zeros(len(rms_z), bool)
    r = rms_z / med
    return (r > lo) & (r < hi)


def hv_from_spectra(sz, sn, se, W):
    """Per-window H/V on the smoothed grid. Horizontals: sqrt((E^2+N^2)/2)."""
    Z = sz @ W.T
    N = sn @ W.T
    E = se @ W.T
    H = np.sqrt((N ** 2 + E ** 2) / 2.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        hv = H / Z
    hv[~np.isfinite(hv)] = np.nan
    return hv, Z, H


def geometric_stats(hv):
    """Geometric mean and geometric standard deviation over windows."""
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.log(np.where(hv > 0, hv, np.nan))
    n = np.sum(np.isfinite(L), axis=0)
    return np.exp(np.nanmean(L, axis=0)), np.exp(np.nanstd(L, axis=0)), n


def find_peaks(fcs, curve, fmin=0.1, fmax=20.0, min_amp=2.0,
               min_prominence=1.3):
    """All local maxima in the band, with their prominence.

    Prominence is the peak amplitude divided by the higher of the two
    adjacent minima. Absolute amplitude alone is not enough: a station whose
    horizontal channels are noisier than its vertical shows an H/V that
    simply climbs towards the end of the band, and the top of that climb is
    a larger number than the real peak without being a peak at all.
    """
    m = (fcs >= fmin) & (fcs <= fmax) & np.isfinite(curve)
    if m.sum() < 5:
        return []
    idx = np.where(m)[0]
    out = []
    for i in idx[1:-1]:
        if not (curve[i] > curve[i - 1] and curve[i] >= curve[i + 1]):
            continue
        if curve[i] < min_amp:
            continue
        lo = i
        while lo > idx[0] and curve[lo - 1] <= curve[lo]:
            lo -= 1
        hi = i
        while hi < idx[-1] and curve[hi + 1] <= curve[hi]:
            hi += 1
        base = max(curve[lo], curve[hi])
        prom = curve[i] / base if base > 0 else np.inf
        if prom >= min_prominence:
            out.append(dict(f=float(fcs[i]), A=float(curve[i]),
                            prominence=float(prom), index=int(i)))
    return out


def refine_peak(fcs, curve, i):
    """Sub-bin peak position by a parabola through the three points, in
    log-frequency.

    Without this, f0 can only take the value of a grid point. On a
    140-point log grid from 0.2 to 25 Hz the spacing is ~3.5 %, so a
    sequence of daily picks looks either perfectly constant or jumps by
    3.5 % -- and a "stability" plot built from it shows quantisation, not
    geology. Returns (f_refined, A_refined).
    """
    if i <= 0 or i >= len(fcs) - 1:
        return float(fcs[i]), float(curve[i])
    y0, y1, y2 = curve[i - 1], curve[i], curve[i + 1]
    den = y0 - 2 * y1 + y2
    if not np.isfinite(den) or den == 0:
        return float(fcs[i]), float(curve[i])
    d = 0.5 * (y0 - y2) / den                 # in bins, within (-1, 1)
    d = float(np.clip(d, -1.0, 1.0))
    lf = np.log10(fcs)
    f = 10 ** (lf[i] + d * (lf[i + 1] - lf[i - 1]) / 2.0)
    A = y1 - 0.25 * (y0 - y2) * d
    return float(f), float(A)


def pick_f0(fcs, curve, fmin=0.1, fmax=20.0, min_amp=2.0,
            min_prominence=1.3, prefer="strongest", refine=True):
    """Fundamental peak, or (nan, nan, -1) if nothing qualifies.

    A curve that merely rises to the edge of the band has its resonance
    outside what the data can see, and reporting the edge as f0 would invent
    a layer.

    `fmin`/`fmax` also serve as a TRACKING band when f0 is picked day after
    day: without one, the pick hops between the fundamental and an overtone
    whenever their amplitudes cross, and the resulting time series looks
    like the site is changing by tens of per cent when it is not.
    """
    pk = find_peaks(fcs, curve, fmin, fmax, min_amp, min_prominence)
    if not pk:
        return np.nan, np.nan, -1
    p = (max(pk, key=lambda d: d["A"]) if prefer == "strongest"
         else min(pk, key=lambda d: d["f"]))
    if refine:
        f, A = refine_peak(fcs, curve, p["index"])
        return f, A, p["index"]
    return p["f"], p["A"], p["index"]


def sesame_checks(fcs, curve, sd, nwin, win_s, f0, A0, i0):
    """The SESAME (2004) clarity criteria that can be evaluated here.

    Reported per station rather than used as a hard filter, so a peak that
    fails one criterion is visible rather than silently dropped.
    """
    out = dict(n_windows=int(nwin), f0=f0, A0=A0)
    if not np.isfinite(f0) or i0 < 0:
        return {**out, "ok_f0_window": False, "ok_amplitude": False,
                "ok_lowside": False, "ok_highside": False, "ok_sigma": False,
                "f_minus": np.nan, "f_plus": np.nan, "sigma_f": np.nan,
                "n_sesame_passed": 0}
    out["ok_f0_window"] = bool(f0 > 10.0 / win_s)
    out["ok_amplitude"] = bool(A0 > 2.0)
    half = A0 / 2.0
    lo = np.where((fcs < f0) & (fcs > f0 / 4.0) & (curve < half))[0]
    hi = np.where((fcs > f0) & (fcs < f0 * 4.0) & (curve < half))[0]
    out["ok_lowside"] = bool(len(lo) > 0)
    out["ok_highside"] = bool(len(hi) > 0)
    out["f_minus"] = float(fcs[lo[-1]]) if len(lo) else np.nan
    out["f_plus"] = float(fcs[hi[0]]) if len(hi) else np.nan
    out["sigma_f"] = float(sd[i0]) if np.isfinite(sd[i0]) else np.nan
    out["ok_sigma"] = bool(np.isfinite(out["sigma_f"]) and out["sigma_f"] < 2.0)
    out["n_sesame_passed"] = int(sum([out["ok_f0_window"], out["ok_amplitude"],
                                      out["ok_lowside"], out["ok_highside"],
                                      out["ok_sigma"]]))
    return out


def compute_hvsr(z, n, e, win_s=60.0, overlap=0.5, b=40.0,
                 fmin=0.1, fmax=25.0, nf=200):
    """End-to-end H/V for three Traces on a common time base.

    Returns a dict with the smoothed frequency grid, the geometric mean
    curve, its geometric standard deviation, the number of windows kept, and
    the picked f0 with its SESAME report.
    """
    dt = 1.0 / z.stats.sampling_rate
    nmin = min(len(z.data), len(n.data), len(e.data))
    sz, freqs, rms = window_spectra(z.data[:nmin], dt, win_s, overlap)
    sn, _, _ = window_spectra(n.data[:nmin], dt, win_s, overlap)
    se, _, _ = window_spectra(e.data[:nmin], dt, win_s, overlap)
    k = min(len(sz), len(sn), len(se))
    if k == 0:
        return None
    keep = select_windows(rms[:k])
    sz, sn, se = sz[:k][keep], sn[:k][keep], se[:k][keep]
    if len(sz) < 5:
        return None
    fcs = np.logspace(np.log10(fmin), np.log10(fmax), nf)
    W = konno_ohmachi_matrix(freqs, fcs, b=b)
    hv, Zs, Hs = hv_from_spectra(sz, sn, se, W)
    mu, sd, nn = geometric_stats(hv)
    f0, A0, i0 = pick_f0(fcs, mu, fmin=max(fmin, 0.2), fmax=min(fmax, 20.0))
    rep = sesame_checks(fcs, mu, sd, len(sz), win_s, f0, A0, i0)
    return dict(f=fcs, hv=mu, sd=sd, n=nn, nwin=len(sz), hv_all=hv,
                f0=f0, A0=A0, i0=i0, sesame=rep,
                Z=np.exp(np.nanmean(np.log(np.where(Zs > 0, Zs, np.nan)), 0)),
                H=np.exp(np.nanmean(np.log(np.where(Hs > 0, Hs, np.nan)), 0)))


def depth_from_f0(f0, vs):
    """Quarter-wavelength depth to the impedance contrast: h = Vs / (4 f0)."""
    return vs / (4.0 * f0)


def layer_transfer(f, h, vs, vs_half, rho=2.0, rho_half=2.5, Q=30.0):
    """Exact SH transfer function of one layer over a half-space.

    Used to check the implementation: the fundamental peak must land at
    Vs/(4h) and the overtones at the odd multiples.
    """
    f = np.asarray(f, dtype=float)
    # NOTE the sign. With +i/(2Q) the imaginary part of cos(kh) subtracts
    # from the impedance term, the denominator passes through zero near
    # sinh(b) = imp, and the transfer function BLOWS UP instead of being
    # attenuated -- giving overtones larger than the fundamental. The
    # physical branch, attenuation with increasing f, is -i/(2Q).
    k = 2 * np.pi * f / vs * (1 - 1j / (2 * Q))
    imp = (rho * vs) / (rho_half * vs_half)
    return np.abs(1.0 / (np.cos(k * h) + 1j * imp * np.sin(k * h)))
