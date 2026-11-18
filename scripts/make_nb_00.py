import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nbtools import md, code, build, BOOTSTRAP, COLAB_BADGE

# ====================================================== 00 download
c = []
c.append(md(f"""# 00 — Rebuilding the dataset from the FDSN archive

{COLAB_BADGE.format(name="00_download_data.ipynb")}

> **Do not run this during the practical.** It takes about two hours. The
> prepared archive is what you install in SETUP.md Part C.

This notebook exists so that the dataset is *reproducible*: every byte in
`data/` can be traced back to a request to a public server, and you can
extend the window or add stations later.

## What the dataset is

Four stations of the Italian National Seismic Network (`IV`), on and around
Etna, between **2025-12-10** and **2026-02-20**:

| station | elevation | role |
|---|---|---|
| `ECPN` | 3038 m | summit (Cratere del Piano) — volcanic tremor, the monitoring target |
| `ESLN` | 1787 m | Serra La Nave — the classic Etna reference site |
| `EMSG` | 1435 m | NW flank — the cleanest continuous record in this window |
| `EFIU` |   97 m | Fiumefreddo, 3 km from the Ionian coast — strongest microseism |

The window was chosen because **three different things happen inside it**:

* **24 Dec 2025** — Etna enters a continuous eruptive phase
* **1 Jan 2026** — a fissure opens at ~2100 m in the Valle del Bove
* **16–23 Jan 2026** — Cyclone Harry: the largest secondary-microseism
  amplitudes recorded in the Ionian since digital recording began in 1990
* **6 Feb 2026** — a new eruption begins at the NE base of the SE Crater

So anything you measure in January has at least two candidate causes, one
volcanic and one meteorological. Separating them is the point of the whole
practical.

Data are open, distributed by INGV through EIDA and FDSN
(DOI [10.13127/sd/x0fxnh7qfy](https://doi.org/10.13127/sd/x0fxnh7qfy)).
"""))

c.append(md("""## Why the download is organised the way it is

INGV's `dataselect` service delivers about **0.09 MB/s on one connection**
(measured 2026-11-16). One station-day of 100 Hz data is ~21 MB, so a single
stream needs about four minutes per station-day.

Three consequences, and they are the real lesson of this notebook:

1. **Ask for small chunks.** Three hours of one channel at a time. Large
   requests time out on the server side, and when they do you lose
   everything you were waiting for.
2. **Write each chunk to disk immediately, and never re-request one you
   already have.** Then an interrupted run is restarted, not repeated.
3. **Run a few connections in parallel, but only a few.** Twelve workers
   measured ~0.7 MB/s overall here. More would be rude to a public service
   that other people depend on, including for earthquake response.

Everything below is in `etnalib/download.py`, which you can read."""))

c.append(code(BOOTSTRAP))

c.append(code('''from etnalib import download as D

jobs = D.build_joblist()
full = D.full_days()
print(f"{len(jobs)} chunks to request")
print(f"{len(full)} days downloaded in full (24 h)")
print(f"other days: {C.PARTIAL_HOURS[0]:02d}:00-{C.PARTIAL_HOURS[1]:02d}:00 UTC only")
print(f"\\nestimated volume: {len(jobs) * 2.6 / 1000:.1f} GB")
print(f"estimated time at 0.7 MB/s: {len(jobs) * 2.6 / 0.7 / 3600:.1f} hours")'''))

c.append(md("""### Whole days, and why that was not the first answer

Every day is now downloaded in full, all 24 hours, all three components —
about 14 GB.

The first version of this dataset took only 00:00–06:00 UTC on most days.
That is a defensible choice and standard practice in monitoring: night-time
hours contain much less anthropogenic noise, and using the *same* hours every
day means that when the measurement changes it is not because the mix of
sources changed with the clock.

It is also a choice that hides something you may want to see. The diurnal
cycle is real, it is large, and it differs completely between bands — in this
dataset the 0.5–5 Hz "tremor" band swings by about ±20 % over the day even at
the summit station, while the primary microseism does not. If the daily cycle
is part of what you are discussing, you need the whole day.

Both are available: `etnalib/config.py` has `ACORR_HOURS`, which is `None`
for whole days and can be set back to `(0, 6)` to reproduce the night-only
results. Notebook 02 compares the two."""))

c.append(code('''# This is the cell that takes two hours. It is safe to interrupt and re-run:
# chunks already on disk are skipped.
#
# D.download_inventory()
# D.download_hvsr_window()
# D.run(D.build_joblist(), workers=12)
print("left commented out on purpose -- see the warning at the top")'''))

c.append(md("""## Step 2: the correlations, with MSNoise

Once the waveforms are on disk, the autocorrelations and single-station
cross-components are computed with **MSNoise 1.6.5** — the community package
for ambient-noise correlation — in its single-station mode:

```bash
bash scripts/run_msnoise.sh
```

That script does the whole chain and is safe to re-run:

1. `build_sds.py` — rewrite the flat chunks as an SDS archive, which is what
   MSNoise knows how to walk
2. `msnoise db init --tech 1` — a SQLite project
3. `msnoise_setup.py` — write the configuration, including
   `components_to_compute_single_station = ZZ,ZN,ZE,NE`, and create the
   1–3 Hz filter (MSNoise 1.6.5 only exposes filter creation through its web
   admin, so this does it through `msnoise.api`)
4. `msnoise populate` **then** `msnoise scan_archive --init`
5. `msnoise new_jobs --init`
6. `msnoise -t 4 compute_cc`
7. `export_msnoise.py` — convert `STACKS/.../*.MSEED` into the `.npz` the
   notebooks read

### Four things that will bite you if you rebuild this

**Order: `populate` before `scan_archive`.** `scan_archive` walks the archive
looking for stations that are *already in the Station table*. Run it first
and it prints "Scanning 0 directories", inserts nothing, and **exits 0**. It
looks completely healthy. This is the reverse of the order in most MSNoise
recipes.

**`sqlalchemy` must be < 2.** MSNoise 1.6.5 calls `engine.execute()`, removed
in SQLAlchemy 2.0, so `msnoise new_jobs` dies with `'Engine' object has no
attribute 'execute'`.

**`setuptools` must be < 81.** Every `msnoise` command imports
`pkg_resources`, which setuptools 81 removed.

**`cc_normalisation` does not exist in 1.6.5.** It is a master-era parameter;
temporal normalisation here is `windsorizing` (3.0 = clip at 3σ, −1 = one-bit,
0 = none).

Both version pins are in `environment.yml` with the failure they prevent
written next to them. `run_msnoise.sh` checks them before it does anything,
because both failures are silent until the stage that hits them."""))

c.append(md("""## Using a different node, or a different volcano

`etnalib/config.py` holds the station list, the window and the bands. To point
this at another network, change `FDSN_BASE` and `STATIONS`. ObsPy's routing
client is the general solution when you do not know which node holds the data:

```python
from obspy.clients.fdsn import RoutingClient
client = RoutingClient("eida-routing")
```

Two other open Etna datasets, if you want to go further — both too heavy for
a one-hour practical, both interesting:

* **`4S` (2022–2023)** — Serra La Nave distributed acoustic sensing, a fibre
  in a 200 m borehole, same site as `ESLN`.
* **`7I` (2025)** — one week of rotational-sensor data on the same concrete
  plate as `ESLN`."""))

build(c, "notebooks/00_download_data.ipynb")
print("wrote 00")
