"""Chunked, restartable, polite downloads from the INGV FDSN node.

Why this file exists at all
---------------------------
INGV's dataselect service delivers roughly 0.09 MB/s on a single connection
(measured 2026-11-16). One station-day of HHZ at 100 Hz is about 21 MB, so a
single stream needs ~4 minutes per station-day and a 70-day, 4-station
request would run for most of a day. Running a dozen connections in parallel
recovers most of that; 12 is polite and was measured at ~0.5 MB/s overall.

Consequences, which are also the lesson:
  * ask for one station and a few hours at a time, never one big request --
    the server times out on large ones;
  * write each chunk to disk as it arrives and never re-request a chunk that
    is already there, so an interrupted run can simply be restarted;
  * do this once, offline, and ship the result. Do not do it in a classroom.
"""
from __future__ import annotations

import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
from pathlib import Path

from . import config as C

UA = "etna-noise-practical/1.0 (teaching material; obspy-compatible)"


def _day_range(start: str, end: str):
    d0 = datetime.strptime(start, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    d1 = datetime.strptime(end, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    while d0 < d1:
        yield d0
        d0 += timedelta(days=1)


def full_days() -> set:
    out = set()
    for a, b in C.FULL_DAY_BLOCKS:
        out |= {d.date() for d in _day_range(a, b)}
    return out


def chunk_path(sta: str, cha: str, day: datetime, hour: int) -> Path:
    p = C.RAW / sta / f"{day:%Y}" / f"{sta}.{cha}.{day:%Y%m%d}.{hour:02d}.mseed"
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def build_joblist(stations=None, channels=(C.CHANNEL,),
                  start=C.START, end=C.END, block_hours=3):
    """Every (station, channel, day, hour, nhours) chunk the practical needs."""
    stations = list(stations or C.STATIONS)
    fd = full_days()
    jobs = []
    for day in _day_range(start, end):
        hours = (range(0, 24, block_hours) if day.date() in fd
                 else range(C.PARTIAL_HOURS[0], C.PARTIAL_HOURS[1], block_hours))
        for sta in stations:
            for cha in channels:
                for h in hours:
                    jobs.append((sta, cha, day, h, block_hours))
    return jobs


def _fetch(job, retries=3, timeout=600):
    sta, cha, day, hour, nh = job
    out = chunk_path(sta, cha, day, hour)
    if out.exists() and out.stat().st_size > 0:
        return ("skip", out, out.stat().st_size)
    t0 = day + timedelta(hours=hour)
    t1 = t0 + timedelta(hours=nh)
    url = (f"{C.FDSN_BASE}/dataselect/1/query?net={C.NETWORK}&sta={sta}"
           f"&cha={cha}&start={t0:%Y-%m-%dT%H:%M:%S}&end={t1:%Y-%m-%dT%H:%M:%S}")
    for k in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = r.read()
            if not data:
                out.write_bytes(b"")          # mark "server says nothing here"
                return ("empty", out, 0)
            tmp = out.with_suffix(".part")
            tmp.write_bytes(data)
            tmp.rename(out)
            return ("ok", out, len(data))
        except urllib.error.HTTPError as e:
            if e.code == 204:                 # no data: legitimate, remember it
                out.write_bytes(b"")
                return ("empty", out, 0)
            if e.code in (400, 404):
                return ("fail", out, 0)
            time.sleep(2 ** k)
        except Exception:
            time.sleep(2 ** k)
    return ("fail", out, 0)


def run(jobs, workers=12, label="download"):
    t0 = time.time()
    done = {"ok": 0, "skip": 0, "empty": 0, "fail": 0}
    nbytes = 0
    with ThreadPoolExecutor(workers) as ex:
        futs = [ex.submit(_fetch, j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            st, _p, n = f.result()
            done[st] += 1
            nbytes += n
            if i % 25 == 0 or i == len(futs):
                el = time.time() - t0
                rate = nbytes / 1e6 / max(el, 1e-9)
                eta = (len(futs) - i) * el / i
                print(f"[{label}] {i}/{len(futs)}  {nbytes/1e6:8.1f} MB  "
                      f"{rate:5.2f} MB/s  ETA {eta/60:5.1f} min  {done}",
                      flush=True)
    print(f"[{label}] finished in {(time.time()-t0)/60:.1f} min: {done}",
          flush=True)
    return done


def download_hvsr_window(stations=None):
    """Two quiet hours, three components, for the HVSR notebook."""
    stations = list(stations or C.STATIONS)
    t0 = datetime.fromisoformat(C.HVSR_WINDOW[0]).replace(tzinfo=timezone.utc)
    t1 = datetime.fromisoformat(C.HVSR_WINDOW[1]).replace(tzinfo=timezone.utc)
    nh = int((t1 - t0).total_seconds() // 3600)
    jobs = [(s, c, t0.replace(hour=0, minute=0), t0.hour, nh)
            for s in stations for c in C.CHANNELS_3C]
    return run(jobs, workers=8, label="hvsr-3c")


def download_inventory(stations=None, path=None):
    stations = ",".join(stations or C.STATIONS)
    path = Path(path or C.PROC / "inventory.xml")
    url = (f"{C.FDSN_BASE}/station/1/query?net={C.NETWORK}&sta={stations}"
           f"&cha=HH?&level=response&format=xml")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=300) as r:
        path.write_bytes(r.read())
    print(f"inventory -> {path} ({path.stat().st_size/1e6:.2f} MB)")
    return path


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "inventory"):
        download_inventory()
    if which in ("all", "hvsr"):
        download_hvsr_window()
    if which in ("all", "waveforms"):
        run(build_joblist(), workers=12, label="HHZ")
