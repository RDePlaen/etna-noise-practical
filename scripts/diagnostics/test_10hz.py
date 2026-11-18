"""Is 10 Hz enough for the 1-2 Hz practical?

Runs the identical MSNoise chain on the same station at 20 Hz and at 10 Hz
and compares the dv/v, stack by stack. Everything else is held fixed:
1-2 Hz filter, maxlag 40 s, 2-day stacks, one-bit, MWCS 5-35 s.

The only configuration that MUST change with the rate is
preprocess_lowpass, which has to stay below the new Nyquist: 8 Hz is fine
at 20 Hz sampling (Nyquist 10) and illegal at 10 Hz sampling (Nyquist 5).
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
MS = os.environ.get("MSNOISE", "msnoise")
PY = sys.executable
NOISY = ("pkg_resources", "SyntaxWarning", "invalid escape", "default[")


def sh(cmd, cwd):
    r = subprocess.run(cmd, cwd=str(cwd), shell=True, capture_output=True,
                       text=True)
    out = [l for l in (r.stdout + r.stderr).splitlines()
           if not any(n in l for n in NOISY)]
    return "\n".join(out[-3:])


def run(tag, sds, fs, lowpass):
    proj = Path("/tmp") / f"ms_{tag}"
    if proj.exists():
        shutil.rmtree(proj)
    proj.mkdir(parents=True)
    sh(f"{MS} db init --tech 1", proj)
    sh(f"{PY} {REPO}/scripts/msnoise_setup.py --sds {sds} "
       f"--response {REPO}/data/processed --components ZZ "
       f"--windsorizing=-1 --fmin 1.0 --fmax 2.0", proj)
    for k, v in (("cc_sampling_rate", fs), ("preprocess_lowpass", lowpass),
                 ("maxlag", 40), ("mov_stack", 2), ("dtt_mincoh", 0.6)):
        sh(f"{MS} config set {k}={v}", proj)
    sh(f"{MS} populate", proj)
    sh(f"{MS} scan_archive --init", proj)
    sh(f"{MS} new_jobs --init", proj)
    sh(f"{MS} -t 4 compute_cc", proj)
    sh(f"{MS} stack -r", proj)
    sh(f"{MS} reset STACK", proj)
    sh(f"{MS} -t 4 stack -m", proj)
    sh(f"{MS} -t 4 compute_mwcs", proj)
    sh(f"{MS} compute_dtt", proj)
    return proj


def read_dtt(proj):
    out = {}
    for f in sorted((proj / "DTT" / "01" / "002_DAYS" / "ZZ").glob("*.txt")):
        try:
            import pandas as pd
            df = pd.read_csv(f, index_col=0)
            if "M" in df.columns and len(df):
                out[f.stem] = float(df["M"].iloc[0]) * 100.0
        except Exception:
            continue
    return out


if __name__ == "__main__":
    sds20 = sys.argv[1] if len(sys.argv) > 1 else str(REPO / "data/sds_teaching")
    sds10 = sys.argv[2] if len(sys.argv) > 2 else "/tmp/sds10"
    print("running 20 Hz ...", flush=True)
    p20 = run("20hz", sds20, 20, 8.0)
    print("running 10 Hz ...", flush=True)
    p10 = run("10hz", sds10, 10, 4.0)

    a, b = read_dtt(p20), read_dtt(p10)
    common = sorted(set(a) & set(b))
    print(f"\n20 Hz: {len(a)} stacks   10 Hz: {len(b)} stacks   "
          f"common: {len(common)}")
    if len(common) < 5:
        sys.exit("too few common stacks to compare")
    va = np.array([a[d] for d in common])
    vb = np.array([b[d] for d in common])
    k = np.isfinite(va) & np.isfinite(vb)
    print(f"correlation      : {np.corrcoef(va[k], vb[k])[0,1]:.4f}")
    print(f"median |20 - 10| : {np.median(np.abs(va[k]-vb[k])):.4f} %")
    print(f"rms  difference  : {np.sqrt(np.mean((va[k]-vb[k])**2)):.4f} %")
    print(f"std of 20 Hz     : {np.std(va[k]):.4f} %")
    print(f"std of 10 Hz     : {np.std(vb[k]):.4f} %")
    print(f"\n{'stack':12s} {'20 Hz':>9s} {'10 Hz':>9s} {'diff':>8s}")
    for d in common[:12]:
        print(f"{d:12s} {a[d]:+9.3f} {b[d]:+9.3f} {a[d]-b[d]:+8.3f}")
