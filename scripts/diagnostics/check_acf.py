import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from etnalib import io, acorr, dvv, config as C

sta = sys.argv[1] if len(sys.argv) > 1 else "EMSG"
days = io.available_days(sta)[:8]
print(f"{sta}: testing on {len(days)} days")

acfs, used = [], []
for d in days:
    t0 = time.time()
    tr = io.load_day(sta, d, fs_out=C.ACORR_FS, min_coverage=0.5,
                     remove_response=False)
    if tr is None:
        print(f"  {d}: rejected")
        continue
    lags, c, nw = acorr.daily_autocorr(tr, C.ACORR_BAND, C.ACORR_MAXLAG,
                                       window_s=C.ACORR_WINDOW,
                                       estimator="ccgn", norm="ram",
                                       stack="linear")
    if c is None:
        print(f"  {d}: no ACF")
        continue
    acfs.append(c)
    used.append(d)
    coda = c[(lags > 3) & (lags < 15)]
    print(f"  {d}: {nw:3d} windows, {time.time()-t0:5.1f}s, "
          f"coda rms {np.sqrt(np.mean(coda**2)):.4f}, "
          f"max|coda| {np.max(np.abs(coda)):.4f}")

acfs = np.vstack(acfs)
ref = acfs.mean(axis=0)
print("\nday-to-day repeatability (correlation of each day with the mean, "
      "3-15 s coda):")
m = (lags > 3) & (lags < 15)
for d, a in zip(used, acfs):
    cc = np.corrcoef(a[m], ref[m])[0, 1]
    print(f"  {d}: cc = {cc:.4f}")

out, _, _ = dvv.dvv_series(acfs, used, lags, (3.0, 15.0))
print("\ndv/v against the mean of these days:")
for d, v, cc in out:
    print(f"  {d}: dv/v = {v:+7.4f} %   cc = {cc:.4f}")
