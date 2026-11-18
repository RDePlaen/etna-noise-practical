import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from etnalib import io, hvsr, config as C

for sta in C.STATIONS:
    tr = io.load_hvsr_window(sta)
    if tr is None:
        print(sta, "none")
        continue
    z, n, e = tr
    r = hvsr.compute_hvsr(z, n, e, win_s=60.0, fmin=0.2, fmax=25.0, nf=160)
    if r is None:
        print(sta, "no result")
        continue
    f0, A0 = r["f0"], r["A0"]
    d300 = hvsr.depth_from_f0(f0, 300.0) if np.isfinite(f0) else float("nan")
    print(f"{sta:6s} elev{C.STATIONS[sta]['elev']:5.0f}m  nwin={r['nwin']:4d}  "
          f"f0={f0:7.3f} Hz  A0={A0:5.2f}  "
          f"SESAME {r['sesame']['n_sesame_passed']}/5   h(Vs=300)={d300:7.1f} m")
