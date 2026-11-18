"""Does MSNoise's single-station ZZ agree with etnalib.acorr?

Two independent implementations of the same thing is the cheapest real test
there is. They will not be identical -- MSNoise whitens and windsorizes with
its own parameters -- but the coda should correlate strongly and the
arrivals must sit at the same lags.
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import obspy
from etnalib import io, acorr, config as C

PROJ = sys.argv[1] if len(sys.argv) > 1 else "/home/claude/etna/msproj"
DAY = sys.argv[2] if len(sys.argv) > 2 else "2026-02-15"
COMP = sys.argv[3] if len(sys.argv) > 3 else "ZZ"


def msnoise_trace(sta, comp=COMP, day=DAY, filt="01"):
    p = glob.glob(f"{PROJ}/STACKS/{filt}/001_DAYS/{comp}/"
                  f"IV_{sta}_IV_{sta}/{day}.MSEED")
    if not p:
        return None, None
    tr = obspy.read(p[0])[0]
    n = tr.stats.npts
    fs = tr.stats.sampling_rate
    lags = (np.arange(n) - (n - 1) // 2) / fs      # symmetric, zero centred
    return lags, tr.data.astype(float)


print(f"{'sta':6s} {'npts':>6s} {'fs':>5s} {'maxlag':>7s} {'sym err':>9s} "
      f"{'cc vs etnalib':>14s} {'peak lag':>9s}")
for sta in C.STATIONS:
    lags, d = msnoise_trace(sta)
    if d is None:
        print(f"{sta:6s}  no MSNoise output")
        continue
    fs = 1.0 / (lags[1] - lags[0])
    # an autocorrelation must be an even function of lag
    half = (len(d) - 1) // 2
    pos, neg = d[half:], d[:half + 1][::-1]
    m = min(len(pos), len(neg))
    sym = np.max(np.abs(pos[:m] - neg[:m])) / np.max(np.abs(d))

    if COMP == "ZZ":
        tr = io.load_day(sta, DAY, fs_out=C.ACORR_FS, min_coverage=0.3,
                         remove_response=False)
        if tr is None:
            print(f"{sta:6s}  no etnalib trace")
            continue
        # match MSNoise's temporal normalisation so we compare like with
        # like: the production run uses windsorizing=-1, i.e. one-bit.
        lg, mine, nw = acorr.daily_autocorr(
            tr, C.ACORR_BAND, C.ACORR_MAXLAG, window_s=C.ACORR_WINDOW,
            norm=os.environ.get("ETNALIB_NORM", "onebit"),
            whiten=True, hours=C.PARTIAL_HOURS)
        if mine is None:
            print(f"{sta:6s}  etnalib produced nothing")
            continue
        theirs = np.interp(lg, lags[half:], pos)
        w = (lg >= 1.0) & (lg <= 15.0)
        cc = np.corrcoef(mine[w], theirs[w])[0, 1]
        pk = lg[w][int(np.argmax(np.abs(theirs[w])))]
        print(f"{sta:6s} {len(d):6d} {fs:5.1f} {lags[-1]:7.1f} {sym:9.2e} "
              f"{cc:14.3f} {pk:9.2f}")
    else:
        print(f"{sta:6s} {len(d):6d} {fs:5.1f} {lags[-1]:7.1f} {sym:9.2e} "
              f"{'n/a':>14s} {lags[int(np.argmax(np.abs(d)))]:9.2f}")
