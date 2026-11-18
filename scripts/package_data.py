#!/usr/bin/env python3
"""Build the student data archive.

The raw archive is ~3 GB of 100 Hz miniSEED. Students do not need it: the
notebooks read precomputed products plus a small amount of decimated
waveform data for the cells where looking at an actual trace is the point.

What goes in:
  processed/            inventory, PPSDs, daily ACFs, band series, HVSR window
  cache/<EFIU, 20 Hz>   the days notebook 01 plots as waveforms/spectrograms

Everything else is reproducible with notebook 00.
"""
import sys
import tarfile
import hashlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from etnalib import config as C  # noqa: E402

SPECTROGRAM_STATION = "EFIU"
SPECTROGRAM_DAYS = [f"2026-01-{d:02d}" for d in range(14, 27)] + ["2026-01-20"]
# Days of waveforms at EVERY station, so that changing STA in notebook 01
# works instead of returning None. ~6 MB per station-day.
#   2026-01-20 -- storm peak, the notebook default; missing at ESLN
#   2026-01-15 -- present at all four stations
COMMON_DAYS = ("2026-01-20", "2026-01-15")
OUT = C.ROOT / "etna_practical_data.tar.gz"


def members():
    got = []
    for p in sorted((C.PROC).glob("*")):
        if p.is_file():
            got.append(p)
    sdst = C.DATA / "sds_teaching"
    if sdst.exists():
        got += [f for f in sdst.rglob("*") if f.is_file()]
    cache = C.DATA / "cache"
    if cache.exists():
        for day in sorted(set(SPECTROGRAM_DAYS)):
            for p in cache.glob(f"{SPECTROGRAM_STATION}.HHZ.{day}.20.resp.*.npz"):
                got.append(p)
        for sta in C.STATIONS:
            for day in COMMON_DAYS:
                got += list(cache.glob(f"{sta}.HHZ.{day}.20.resp.*.npz"))
    return got


def main():
    files = sorted(set(members()))
    if not files:
        sys.exit("nothing to package -- run build_products.py first")
    total = sum(f.stat().st_size for f in files)
    print(f"{len(files)} files, {total/1e6:.1f} MB uncompressed")
    with tarfile.open(OUT, "w:gz", compresslevel=6) as t:
        for f in files:
            arc = Path("data") / f.relative_to(C.DATA)
            t.add(f, arcname=str(arc))
            print(f"  {arc}  ({f.stat().st_size/1e6:.1f} MB)")
    size = OUT.stat().st_size
    sha = hashlib.sha256(OUT.read_bytes()).hexdigest()
    print(f"\n{OUT.name}: {size/1e6:.1f} MB")
    print(f"sha256: {sha}")
    print("\nUpload to Zenodo, then put the record URL and this checksum into"
          " scripts/get_data.py")


if __name__ == "__main__":
    main()
