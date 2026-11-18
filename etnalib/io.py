"""Loading and preprocessing the downloaded chunks.

The raw archive is a flat set of hourly-ish chunks written by download.py.
Everything downstream asks for "one station, one day" and gets back a single
ObsPy Trace, already response-corrected and decimated, or None if that day is
unusable. Deciding what "unusable" means is part of the exercise, so the rule
is explicit and in one place: see `load_day`.
"""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import obspy
from obspy import UTCDateTime

from . import config as C

_INV = None


def inventory():
    """The StationXML, read once and cached."""
    global _INV
    if _INV is None:
        p = C.PROC / "inventory.xml"
        if not p.exists():
            raise FileNotFoundError(
                f"{p} is missing. Run notebook 00, or unpack the data archive."
            )
        _INV = obspy.read_inventory(str(p))
    return _INV


def chunk_files(sta: str, day, cha: str = C.CHANNEL) -> list[Path]:
    if isinstance(day, str):
        day = datetime.strptime(day, "%Y-%m-%d").date()
    d = C.RAW / sta / f"{day:%Y}"
    if not d.exists():
        return []
    return sorted(d.glob(f"{sta}.{cha}.{day:%Y%m%d}.*.mseed"))


def load_day(sta: str, day, cha: str = C.CHANNEL,
             fs_out: float = C.ACORR_FS,
             remove_response: bool = True,
             min_coverage: float = 0.5,
             prefilt=(0.02, 0.05, 8.0, 9.5),
             cache: bool = True):
    """One station-day as a single gap-free Trace, or None.

    Steps, in the order that matters:

    1. read every chunk for that day and merge, filling gaps with zeros
       (`fill_value=0`). Zeros are honest here: they contribute no energy and
       they keep the time axis rigid, which every later step assumes.
    2. reject the day if less than `min_coverage` of it is real data. A day
       that is 80 % zeros will still produce a beautiful-looking, meaningless
       autocorrelation, so this test is not optional.
    3. detrend, taper
    4. anti-alias low-pass and decimate to `fs_out`
    5. optionally remove the instrument response, to velocity

    **Decimate before deconvolving.** Removing the response evaluates the
    instrument transfer function at every FFT frequency, so the cost scales
    with the number of samples -- at 100 Hz a single day takes minutes, at
    20 Hz it takes seconds. The anti-alias filter is applied to raw counts,
    which is harmless, and the pre-filter corners sit well inside the new
    Nyquist. Doing it the other way round is the single easiest way to make
    an ambient-noise pipeline unusably slow.

    **When you need the response removed, and when you do not.** For band
    amplitudes in m/s, you need it. For a single-station AUTOCORRELATION you
    do not, and this is exact rather than approximate: the correlation of a
    trace with itself is

        C(f) = H*(f) H(f) S(f) = |H(f)|^2 S(f)

    so the instrument PHASE cancels identically, and a phase error is the
    travel-time error that monitoring would otherwise measure. Only the
    amplitude is coloured, and whitening flattens that anyway. (This is not
    true of cross-correlation between two stations with different sensors,
    where the phases do not cancel and the response must be removed.)
    """
    # The cache is consulted BEFORE the raw archive, because the packaged
    # student dataset ships cached days and no raw miniSEED at all.
    #
    # The key carries the NUMBER OF SOURCE CHUNKS. Without it, a cache
    # written when only 00:00-06:00 was on disk keeps being served after the
    # rest of the day is downloaded, and every product silently rebuilds
    # from six hours while reporting success. Including the count means new
    # data invalidates the entry by itself.
    if isinstance(day, str):
        _dkey = day
    else:
        _dkey = f"{day:%Y-%m-%d}"
    CACHE_VERSION = 3
    files = [f for f in chunk_files(sta, day, cha) if f.stat().st_size > 0]
    _rate = int(fs_out) if fs_out < 1e6 else "raw"
    _resp = "resp" if remove_response else "counts"
    _base = f"{sta}.{cha}.{_dkey}.{_rate}.{_resp}"
    cdir = C.DATA / "cache"
    ckey = cdir / f"{_base}.n{len(files)}.v{CACHE_VERSION}.npz"

    if cache:
        hit = ckey if ckey.exists() else None
        if hit is None and not files:
            # no raw at all: this is a student install, take whatever
            # chunk count the packaged cache was built with
            g = sorted(cdir.glob(f"{_base}.n*.v{CACHE_VERSION}.npz")) \
                if cdir.exists() else []
            hit = g[0] if g else None
        if hit is not None:
            try:
                z = np.load(hit)
                tr = obspy.Trace(z["data"])
                tr.stats.network, tr.stats.station, tr.stats.channel = \
                    C.NETWORK, sta, cha
                tr.stats.sampling_rate = float(z["fs"])
                tr.stats.starttime = UTCDateTime(str(z["t0"]))
                tr.stats.coverage = float(z["coverage"])
                tr.stats.day = _dkey
                return tr
            except Exception:
                pass

    if not files:
        return None

    st = obspy.Stream()
    for f in files:
        try:
            st += obspy.read(str(f), format="MSEED")
        except Exception:
            continue
    if len(st) == 0:
        return None

    # how much real data is there, before any zero filling?
    real = sum(tr.stats.npts / tr.stats.sampling_rate for tr in st)

    st.merge(method=1, fill_value=0)
    if len(st) != 1:
        st = obspy.Stream([st[0]])
    tr = st[0]

    # the expected span is the part of the day we actually asked for
    hours = sorted({int(f.stem.split(".")[-1]) for f in files})
    expected = (len(hours) * 3) * 3600.0        # chunks are 3 h
    if expected <= 0 or real / expected < min_coverage:
        return None

    tr.detrend("demean")
    tr.detrend("linear")
    # Taper a FIXED 60 s, not a percentage. `taper(0.02)` on a 24 h trace
    # tapers 29 minutes at each end -- about one analysis segment -- so the
    # first and last segment of every day come out artificially quiet, and a
    # band-amplitude time series grows a downward spike at every midnight.
    # Sixty seconds is ample for the filters that follow.
    tr.taper(max_percentage=None, max_length=60.0, type="cosine")
    tr.data = tr.data.astype(np.float64)

    # decimate FIRST -- see the docstring
    if tr.stats.sampling_rate > fs_out:
        tr.filter("lowpass", freq=0.45 * fs_out, corners=4, zerophase=True)
        factor = tr.stats.sampling_rate / fs_out
        if abs(factor - round(factor)) < 1e-9:
            tr.decimate(int(round(factor)), no_filter=True)
        else:
            tr.resample(fs_out, no_filter=True)

    if remove_response:
        nyq = 0.5 * tr.stats.sampling_rate
        pf = tuple(min(p, 0.95 * nyq) for p in prefilt)
        try:
            tr.remove_response(inventory=inventory(), output="VEL",
                               pre_filt=pf, water_level=60, zero_mean=False,
                               taper=False)
        except Exception:
            return None

    tr.stats.coverage = float(real / expected)
    tr.stats.day = _dkey
    if cache:
        ckey.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(ckey, data=tr.data.astype(np.float32),
                            fs=tr.stats.sampling_rate,
                            t0=str(tr.stats.starttime),
                            coverage=tr.stats.coverage)
    return tr


def require_day(sta: str, day, cha: str = C.CHANNEL, **kw):
    """`load_day`, but it explains itself when there is nothing to load.

    `load_day` returns None for several different reasons and the caller
    cannot tell them apart, so a missing station surfaces three lines later
    as `'NoneType' object has no attribute 'stats'` -- which says nothing
    about what went wrong or what to do. This raises instead, and says
    which stations and days actually have waveforms.
    """
    tr = load_day(sta, day, cha=cha, **kw)
    if tr is not None:
        return tr

    if sta not in C.STATIONS:
        raise ValueError(f"{sta!r} is not one of {', '.join(C.STATIONS)}")

    have = available_days(sta, cha)
    if have:
        if len(have) <= 4:
            which = "it does have " + ", ".join(have)
        else:
            which = (f"it does have {len(have)} other days, "
                     f"{have[0]} to {have[-1]}")
        raise FileNotFoundError(
            f"{sta} {cha} has no usable data for {day}; {which}.")

    others = {s: available_days(s, cha) for s in C.STATIONS}
    others = {s: d for s, d in others.items() if d}
    msg = [f"No {cha} waveforms for {sta} at all."]
    if others:
        msg.append("")
        msg.append("The prepared dataset ships raw waveforms for:")
        for s, d in others.items():
            msg.append(f"    {s}: {len(d)} days, {d[0]} to {d[-1]}")
        msg.append("")
        msg.append(
            "Precomputed products (band amplitudes, PPSDs, H/V,\n"
            "autocorrelations) cover ALL four stations -- it is only the raw\n"
            "waveforms that are limited, to keep the download small. To get\n"
            "the rest, run notebook 00.")
    raise FileNotFoundError("\n".join(msg))


def available_days(sta: str, cha: str = C.CHANNEL) -> list[str]:
    """Days this station has, from the raw archive or the packaged cache.

    Students have only the cache, so both sources are searched and the union
    returned.
    """
    days = set()
    d = C.RAW / sta
    if d.exists():
        days |= {f.stem.split(".")[2] for y in d.iterdir() if y.is_dir()
                 for f in y.glob(f"{sta}.{cha}.*.mseed") if f.stat().st_size > 0}
    out = {f"{s[:4]}-{s[4:6]}-{s[6:]}" for s in days}

    cdir = C.DATA / "cache"
    if cdir.exists():
        for f in cdir.glob(f"{sta}.{cha}.*.npz"):
            parts = f.name.split(".")
            if len(parts) > 2 and len(parts[2]) == 10:
                out.add(parts[2])
    return sorted(out)


def load_hvsr_window(sta: str):
    """The 3-component HVSR window for one station, as (Z, N, E) traces.

    Reads the packaged `processed/hvsr_ready.npz` if it is there, and falls
    back to the raw archive otherwise. Students have the first and not the
    second; whoever rebuilt the dataset has both.
    """
    ready = C.PROC / "hvsr_ready.npz"
    if ready.exists():
        z = np.load(ready)
        if f"{sta}_Z" in z:
            fs = float(z[f"{sta}_fs"][0])
            out = []
            for lab, cha in zip("ZNE", C.CHANNELS_3C):
                tr = obspy.Trace(z[f"{sta}_{lab}"].astype(np.float64))
                tr.stats.network, tr.stats.station = C.NETWORK, sta
                tr.stats.channel = cha
                tr.stats.sampling_rate = fs
                tr.stats.starttime = UTCDateTime(C.HVSR_WINDOW[0])
                out.append(tr)
            return tuple(out)

    t0 = UTCDateTime(C.HVSR_WINDOW[0])
    day = datetime.fromisoformat(C.HVSR_WINDOW[0]).date()
    out = {}
    for cha in C.CHANNELS_3C:
        files = [f for f in chunk_files(sta, day, cha) if f.stat().st_size > 0]
        if not files:
            return None
        st = obspy.Stream()
        for f in files:
            st += obspy.read(str(f), format="MSEED")
        st.merge(method=1, fill_value=0)
        tr = st[0]
        tr.detrend("demean")
        tr.detrend("linear")
        tr.taper(0.02)
        # NOTE: no response removal. H/V is a ratio of channels from the SAME
        # three-component sensor, so the response cancels to the accuracy of
        # the channel-to-channel match. Removing it changes nothing and can
        # introduce low-frequency artefacts.
        out[cha[-1]] = tr
    return out["Z"], out["N"], out["E"]
