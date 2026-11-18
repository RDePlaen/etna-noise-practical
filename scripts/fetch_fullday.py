"""Fill every day out to 24 h.

The original download took whole days only inside the three event blocks and
00:00-06:00 elsewhere, and the horizontals only ever for the night window.
Diurnal variation is part of the discussion, so every day needs every hour.

Phase 1 is HHZ (what dv/v needs); phase 2 is HHN/HHE (what continuous HVSR
needs). Phase 1 first, so an interrupted run still leaves dv/v complete.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from etnalib import config as C, download as D   # noqa: E402

phase = sys.argv[1] if len(sys.argv) > 1 else "z"
stations = sys.argv[2].split(",") if len(sys.argv) > 2 else list(C.STATIONS)
chans = ("HHZ",) if phase == "z" else ("HHN", "HHE")

jobs = []
for day in D._day_range(C.START, C.END):
    for sta in stations:
        for cha in chans:
            for h in range(0, 24, 3):
                p = D.chunk_path(sta, cha, day, h)
                if p.exists():          # already there (or known empty)
                    continue
                jobs.append((sta, cha, day, h, 3))

print(f"phase {phase}: {len(jobs)} chunks missing, ~{len(jobs)*2.6/1000:.1f} GB, "
      f"~{len(jobs)*2.6/0.75/3600:.1f} h", flush=True)
D.run(jobs, workers=12, label=f"fullday-{phase}")
