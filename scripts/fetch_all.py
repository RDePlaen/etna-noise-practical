"""Priority-ordered fetch: the most valuable data first, so an interrupted
run still leaves a usable practical."""
import sys, time
sys.path.insert(0, "/home/claude/etna/repo")
from etnalib import config as C, download as D

t0 = time.time()
D.download_inventory()
D.download_hvsr_window()

jobs = D.build_joblist()
def prio(j):
    d = j[2].date().isoformat()
    if "2026-01-14" <= d <= "2026-01-27": return 0   # Harry
    if "2025-12-28" <= d <= "2026-01-03": return 1   # New Year fissure
    if "2026-02-04" <= d <= "2026-02-09": return 2   # February eruption
    return 3                                         # night-only fill
jobs.sort(key=prio)
for p, lab in enumerate(["harry", "newyear", "february", "fill"]):
    sub = [j for j in jobs if prio(j) == p]
    print(f"\n=== phase {p}: {lab}  ({len(sub)} chunks) ===", flush=True)
    D.run(sub, workers=12, label=lab)
print(f"\nTOTAL {(time.time()-t0)/60:.1f} min", flush=True)
