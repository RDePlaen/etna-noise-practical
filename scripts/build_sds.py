#!/usr/bin/env python3
"""Convert the flat chunk archive into an SDS tree for MSNoise.

MSNoise discovers data by walking a structured archive. The download writes
flat 3-hour chunks, so they are merged per station-day-channel and written as

    <SDS>/<YEAR>/<NET>/<STA>/<CHAN>.D/<NET>.<STA>.<LOC>.<CHAN>.D.<YEAR>.<DOY>

Nothing is resampled or deconvolved here: MSNoise does its own preprocessing
and expects raw counts.

Usage:
    python scripts/build_sds.py                 # every station/channel/day
    python scripts/build_sds.py --day 2026-02-15 --channels HHZ,HHN,HHE
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import obspy

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from etnalib import config as C, io  # noqa: E402

SDS = C.DATA / "sds"


def sds_path(net, sta, cha, t):
    """SDS file for the day `t` BELONGS TO.

    Take the day from the day being converted, never from the first trace's
    starttime: a chunk can begin a few seconds before midnight, which would
    file a whole day under the previous one and leave MSNoise seeing two
    half days.
    """
    p = (SDS / f"{t.year}" / net / sta / f"{cha}.D" /
         f"{net}.{sta}..{cha}.D.{t.year}.{t.julday:03d}")
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def convert(stations, channels, days=None, force=False, hours=None):
    """hours=(h0,h1) trims every day to the same UTC window.

    This is not a space optimisation, it is a correctness one. The archive
    holds 24 h for the days around the three events and 00:00-06:00 for the
    rest. Correlating a 24 h day and a 6 h day and then comparing them means
    comparing different mixtures of the noise field: the ambient field has a
    strong daily cycle, so the two are simply different waveforms. Measured
    earlier in this project, mixing them drops the day-to-day correlation
    coefficient from ~0.9 to ~0.4 and looks exactly like bad data.
    """
    n_written = n_skipped = 0
    for sta in stations:
        for cha in channels:
            dd = days or sorted({f.stem.split(".")[2] for y in
                                 (C.RAW / sta).iterdir() if y.is_dir()
                                 for f in y.glob(f"{sta}.{cha}.*.mseed")
                                 if f.stat().st_size > 0}) if (C.RAW / sta).exists() else []
            dd = [f"{d[:4]}-{d[4:6]}-{d[6:]}" if "-" not in d else d for d in dd]
            for day in dd:
                files = [f for f in io.chunk_files(sta, day, cha)
                         if f.stat().st_size > 0]
                if not files:
                    continue
                st = obspy.Stream()
                for f in files:
                    try:
                        st += obspy.read(str(f), format="MSEED")
                    except Exception:
                        continue
                if len(st) == 0:
                    continue
                t0 = obspy.UTCDateTime(f"{day}T00:00:00")
                out = sds_path(C.NETWORK, sta, cha, t0)
                if out.exists() and not force:
                    n_skipped += 1
                    continue
                # Trim to the day this file is named for, then split: merging
                # across a gap produces a MASKED array, which miniSEED cannot
                # store. split() writes each contiguous segment as its own
                # records, which is what an SDS day file should contain.
                st.merge(method=1, fill_value=None)
                a, b = (t0, t0 + 86400.0) if hours is None else \
                       (t0 + hours[0] * 3600.0, t0 + hours[1] * 3600.0)
                st = st.slice(a, b, nearest_sample=False).split()
                st = obspy.Stream([tr for tr in st if tr.stats.npts > 0])
                if len(st) == 0:
                    continue
                for tr in st:
                    if tr.data.dtype.kind == "f":
                        tr.data = tr.data.astype("int32")
                st.write(str(out), format="MSEED", encoding="STEIM2",
                         reclen=4096)
                n_written += 1
    return n_written, n_skipped


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", default=",".join(C.STATIONS))
    ap.add_argument("--channels", default="HHZ,HHN,HHE")
    ap.add_argument("--day", default=None, help="one day, YYYY-MM-DD")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--hours", default="0,6",
                    help="UTC window kept for every day, e.g. 0,6; "
                         "'all' keeps whole days")
    a = ap.parse_args()
    hrs = None if a.hours == "all" else tuple(int(x) for x in a.hours.split(","))
    w, s = convert(a.stations.split(","), a.channels.split(","),
                   [a.day] if a.day else None, a.force, hours=hrs)
    print(f"SDS at {SDS}: {w} files written, {s} already there")
