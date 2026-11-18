#!/usr/bin/env python3
"""A small SDS archive so students can run MSNoise themselves in the session.

The full archive is 2.4 GB. This is one station, one channel, already
decimated to 20 Hz, night hours only -- about 40 MB, which fits in the data
download and still runs the complete MSNoise chain.

Decimating before MSNoise sees the data is a deliberate shortcut and worth
being honest about in the notebook: MSNoise would normally resample from
100 Hz itself. Here it finds the data already at its target rate and the
resampling step becomes a no-op.

**10 Hz, not 20.** The practical works in 1-2 Hz, so a 5 Hz Nyquist is
ample, and it halves the archive (213 MB -> 107 MB). Measured against an
identical 20 Hz run of the whole MSNoise chain on the same station
(`scripts/diagnostics/test_10hz.py`): 65 stacks each, dv/v correlation
0.965, median absolute difference 0.018 % against a signal standard
deviation of 0.083 % and MWCS errors around 0.13 %. The rate is not what
limits this measurement.

The one thing that MUST change with it is `preprocess_lowpass`, which has to
stay below Nyquist: 8 Hz is fine at 20 Hz sampling and illegal at 10 Hz.
Notebook 03 sets it to 4 Hz. The cost is that no band above ~4 Hz can be
explored from this archive; the full 100 Hz data is in `data/raw/`.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import obspy

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from etnalib import config as C, io  # noqa: E402

OUT = C.DATA / "sds_teaching"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--station", default="ECPN")
    ap.add_argument("--channel", default="HHZ")
    ap.add_argument("--fs", type=float, default=10.0)
    ap.add_argument("--hours", default="0,6")
    a = ap.parse_args()
    h0, h1 = (int(x) for x in a.hours.split(","))

    days = io.available_days(a.station, a.channel)
    n = 0
    for d in days:
        # remove_response=False: MSNoise is given counts, as it expects, and
        # for a single-station autocorrelation the response cancels anyway.
        tr = io.load_day(a.station, d, cha=a.channel, fs_out=a.fs,
                         remove_response=False, min_coverage=0.3)
        if tr is None:
            continue
        t0 = obspy.UTCDateTime(f"{d}T00:00:00")
        seg = tr.slice(t0 + h0 * 3600.0, t0 + h1 * 3600.0)
        if seg.stats.npts == 0:
            continue
        seg = seg.copy()
        seg.data = np.nan_to_num(seg.data).astype(np.int32)
        p = (OUT / f"{t0.year}" / C.NETWORK / a.station / f"{a.channel}.D" /
             f"{C.NETWORK}.{a.station}..{a.channel}.D.{t0.year}.{t0.julday:03d}")
        p.parent.mkdir(parents=True, exist_ok=True)
        obspy.Stream([seg]).write(str(p), format="MSEED", encoding="STEIM2",
                                  reclen=4096)
        n += 1
    size = sum(f.stat().st_size for f in OUT.rglob("*") if f.is_file())
    print(f"{n} day files -> {OUT}  ({size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
