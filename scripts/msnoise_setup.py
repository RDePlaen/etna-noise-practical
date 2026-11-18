#!/usr/bin/env python3
"""Configure an MSNoise 1.6.5 project for single-station AC and SC.

MSNoise 1.6.5 only exposes filter creation through its web admin, so this
does it through `msnoise.api` instead, which keeps the whole run scriptable
and reproducible.

What this sets up
-----------------
* no inter-station correlations at all (`components_to_compute` empty)
* `components_to_compute_single_station = ZZ,ZN,ZE,NE`
      ZZ -> autocorrelation of the vertical
      ZN, ZE, NE -> single-station cross-components
* one filter covering the teaching band

Run from inside the MSNoise project directory.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from etnalib import config as C  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sds", default=str(C.DATA / "sds"))
    ap.add_argument("--response", default=str(C.PROC))
    ap.add_argument("--start", default=C.START)
    ap.add_argument("--end", default=C.END)
    ap.add_argument("--components", default="ZZ,ZN,ZE,NE")
    ap.add_argument("--cc-type", default="CC", choices=["CC", "PCC"],
                    help="autocorrelation estimator. PCC is the "
                         "amplitude-unbiased phase cross-correlation of "
                         "Schimmel (1999); it needs the `phasecorr` package, "
                         "which MSNoise does not declare as a dependency. "
                         "Note 1.6.5 implements PCC for AC only, not SC.")
    ap.add_argument("--windsorizing", default="3.0",
                    help="temporal normalisation: 3.0 = clip at 3 sigma, "
                         "-1 = one-bit, 0 = none")
    ap.add_argument("--fmin", type=float, default=C.ACORR_BAND[0])
    ap.add_argument("--fmax", type=float, default=C.ACORR_BAND[1])
    a = ap.parse_args()

    from msnoise.api import connect, update_config, update_filter, get_filters

    db = connect()

    cfg = {
        # ---- archive
        "data_folder": a.sds,
        "data_structure": "SDS",
        "response_path": a.response,
        "startdate": a.start,
        "enddate": a.end,
        "channels": "HH*",

        # ---- preprocessing: match etnalib exactly so the two are comparable
        "cc_sampling_rate": str(C.ACORR_FS),
        "resampling_method": "Lanczos",
        "preprocess_lowpass": "8.0",
        "preprocess_highpass": "0.05",
        # Response removal OFF. For a single station the instrument transfer
        # function cancels in the autocorrelation, C(f)=|H|^2 S(f), and for
        # the cross-components both channels share one sensor and one
        # digitiser, so the phase cancels there too.
        "remove_response": "N",

        # ---- what to correlate
        "components_to_compute": "",            # nothing between stations
        "components_to_compute_single_station": a.components,
        "cc_type_single_station_AC": a.cc_type,
        "cc_type_single_station_SC": "CC",
        "autocorr": "N",                        # deprecated flag, keep off

        # ---- correlation
        "maxlag": str(C.ACORR_MAXLAG),
        "corr_duration": str(C.ACORR_WINDOW),
        "overlap": "0.5",
        "windsorizing": a.windsorizing,
        "whitening": "A",
        "whitening_type": "B",
        # NOTE: cc_normalisation does not exist in 1.6.5 (it is a
        # master-era parameter); temporal normalisation there is controlled
        # by `windsorizing` -- 3.0 clips at 3 sigma, 0 = none, -1 = 1-bit.
        "keep_days": "Y",
        "keep_all": "N",

        # ---- reference stack window
        "ref_begin": a.start,
        "ref_end": a.end,
        "stack_method": "linear",
        "mov_stack": "1",
    }

    for k, v in cfg.items():
        try:
            update_config(db, k, v)
        except Exception as e:
            print(f"  ! could not set {k}: {e}")

    # one filter. mwcs_* are required columns even though we do not run MWCS.
    update_filter(db, ref=1, low=a.fmin, high=a.fmax,
                  mwcs_low=a.fmin, mwcs_high=a.fmax,
                  rms_threshold=0, mwcs_wlen=12, mwcs_step=4, used=True)

    print("config written. filters:")
    for f in get_filters(db, all=True):
        print(f"  filter {f.ref}: {f.low}-{f.high} Hz  used={f.used}")
    print(f"single-station components: {a.components}")
    print(f"AC estimator: {a.cc_type}   windsorizing: {a.windsorizing}")
    db.close()


if __name__ == "__main__":
    main()
