"""Download HHN and HHE for the whole window, night hours only.

Single-station cross-components (EN, EZ, NZ) need the horizontals, and the
original download took HHZ only. The autocorrelation work uses 00:00-06:00
UTC on every day for comparability, so the horizontals are only fetched for
those hours -- a third of the volume of fetching whole days, and nothing
downstream can use the rest.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from etnalib import config as C, download as D  # noqa: E402

jobs = []
for day in D._day_range(C.START, C.END):
    for sta in C.STATIONS:
        for cha in ("HHN", "HHE"):
            for h in range(C.PARTIAL_HOURS[0], C.PARTIAL_HOURS[1], 3):
                jobs.append((sta, cha, day, h, 3))

print(f"{len(jobs)} chunks, ~{len(jobs)*2.6/1000:.1f} GB, "
      f"~{len(jobs)*2.6/0.7/60:.0f} min at 0.7 MB/s", flush=True)
D.run(jobs, workers=12, label="HHN/HHE")
