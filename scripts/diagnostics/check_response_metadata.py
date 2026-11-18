"""Cross-check each response's stated sensitivity against the product of its
stage gains. A mismatch means the StationXML is internally inconsistent, and
anything that uses the scalar sensitivity will be wrong by that factor.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from obspy import UTCDateTime
from etnalib import io, config as C

inv = io.inventory()
t = UTCDateTime("2026-01-20T00:00:00")
print(f"{'sta':6s} {'stated sens':>13s} {'prod of stages':>15s} {'ratio':>7s} "
      f"{'nstages':>7s}  sensor")
for sta in C.STATIONS:
    try:
        r = inv.get_response(f"{C.NETWORK}.{sta}..{C.CHANNEL}", t)
    except Exception as e:
        print(f"{sta}: {e}")
        continue
    stated = r.instrument_sensitivity.value
    prod = 1.0
    for s in r.response_stages:
        if s.stage_gain is not None:
            prod *= s.stage_gain
    desc = ""
    for net in inv:
        for st in net:
            if st.code == sta:
                for ch in st:
                    if ch.code == C.CHANNEL and ch.sensor is not None:
                        desc = ch.sensor.description or ""
    print(f"{sta:6s} {stated:13.4e} {prod:15.4e} {prod/stated:7.4f} "
          f"{len(r.response_stages):7d}  {desc[:30]}")
    if abs(prod / stated - 1.0) > 0.01:
        print(f"       ^ MISMATCH {100*(prod/stated-1):+.1f} %  "
              f"stage gains: {[s.stage_gain for s in r.response_stages]}")
        print(f"       sensitivity quoted at {r.instrument_sensitivity.frequency} Hz")
