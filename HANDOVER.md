# What to do in the next two days

## Tomorrow morning — the one thing that is time-critical

1. **Create the GitHub repo** and push this directory.
2. **Upload `etna_practical_data.tar.gz` (251 MB) to Zenodo**, get the record
   URL, and put it into `scripts/get_data.py`:
   ```python
   URL = "https://zenodo.org/records/<ID>/files/etna_practical_data.tar.gz"
   SHA256 = "5622ad0b007e057e3ac47c8fb240ced2b85491e755b44a75e159f82b8072afe3"
   ```
   (That checksum is the archive as built on 2026-11-16. If you rebuild it,
   `scripts/package_data.py` prints a new one.)
3. **Replace `ORG` with your GitHub org/user** in two places — the Colab
   badge in `scripts/nbtools.py` and the clone line in `SETUP.md` Part D —
   then regenerate the notebooks:
   ```
   for n in 00 01 02 03 04; do python scripts/make_nb_$n.py; done
   ```
4. **Send the students `SETUP.md`.** That is the deliverable for tomorrow.
   Everything else can slip; this cannot.

## Before the session

* Run `conda env create -f environment.yml` **on a machine that is not this
  one**, ideally Windows, and check `python tests/test_acorr.py` passes. The
  environment is pinned and the packages were verified present on
  conda-forge for linux-64, osx-64, osx-arm64 and win-64, but an actual
  solve on a student-like machine is worth twenty minutes.
* Open one notebook in **Colab** from the pushed repo and run the first cell.
  Colab's base image changes without notice; this is the one check that
  cannot be done in advance from here.
* Read `TEACHING.md`. It has the minute-by-minute plan, the group
  assignments for the notebook 03 exercise, and the answers to the questions
  students actually ask.

## MSNoise

The autocorrelations and single-station cross-components are computed with
**MSNoise 1.6.5** (the current release — not master). `scripts/run_msnoise.sh`
runs the whole chain; `scripts/export_msnoise.py` converts its output into the
`.npz` the notebooks read.

**One-bit normalisation, not the 3σ default.** This was the single biggest
lever in the whole build, measured on this dataset:

| `windsorizing` | usable days (cc>0.7) | median cc | saturated dv/v |
|---|---|---|---|
| `3.0` (clip at 3σ, the usual default) | 21/269 | 0.24–0.50 | 23 |
| `-1` (one-bit) | **128/269** | **0.65–0.78** | **0** |

Your own `etna-local-msnoise-run.md` §8 predicted exactly this: Etna's tremor
is a persistent, strongly localised source at 0.5–5 Hz radiating from the
summit craters, which violates the diffuse-field assumption, and aggressive
temporal normalisation is the standard defence. `run_msnoise.sh` defaults to
`-1`; override with `WINDSORIZING=3.0 bash scripts/run_msnoise.sh` to
reproduce the comparison. `scripts/diagnostics/elect_normalisation.py`
prints the table.

Note that MSNoise 1.6.5 offers no running-absolute-mean option — only
clip-at-Nσ, one-bit and none. `etnalib`'s RAM normalisation still produced
somewhat more usable days than one-bit on this dataset, so if you ever want
the last few percent, RAM is worth adding upstream.

**Two version pins are load-bearing and both failures are silent until you
hit them:**

| pin | what happens without it |
|---|---|
| `sqlalchemy<2` | `msnoise new_jobs` dies: `'Engine' object has no attribute 'execute'` (api.py:788 uses the 1.x API) |
| `setuptools<81` | every `msnoise` command dies at import: `No module named 'pkg_resources'` |

Both are in `environment.yml` with the failure written beside them, and
`run_msnoise.sh` refuses to start without them.

**Not master, deliberately.** Master has the stages your own
`etna-msnoise-handoff.md` documents as broken, and a git commit cannot be
pinned in a file a first-time conda user has to run. 1.6.5 supports exactly
what is needed here: `components_to_compute_single_station = ZZ,ZN,ZE,NE`.

**`environment-minimal.yml` is the fallback** if the solve goes badly on a
student machine — same environment without MSNoise. Nothing in notebooks
01–04 imports it, so the practical is unchanged.

## Whole days (17 Nov)

The dataset originally held 00:00-06:00 UTC on most days. It now holds **all
24 hours of every day** for HHZ, HHN and HHE (~14 GB raw, 12 GB re-downloaded
in about 2 hours), and dv/v, H/V, the band series and the PPSD were all
rebuilt from it. `ACORR_HOURS` in `etnalib/config.py` switches back to
`(0, 6)` to reproduce the night-only results.

One bug this turned up, and it is the kind that does not announce itself: the
day cache keyed on station/day/rate but **not on how much data was on disk
when it was written**, so after re-downloading, every product silently
rebuilt from the cached six hours while reporting success. The key now
carries the source-chunk count, so new data invalidates it. If you ever add
hours to an existing day, check that `load_day` returns the hours you expect
before trusting anything downstream.

The student archive is **251 MB**, most of it the 24-hour `sds_teaching/`
that notebook 03 runs MSNoise on. That archive is decimated to **10 Hz**,
not 20: the practical works in 1-2 Hz, so a 5 Hz Nyquist is ample, and it
halves the download. Verified by running the identical MSNoise chain at both
rates on the same station (`scripts/diagnostics/test_10hz.py`): 65 stacks
each, dv/v correlation **0.965**, median absolute difference **0.018 %**
against a signal standard deviation of 0.083 % and MWCS errors near 0.13 %.
The one setting that must follow the rate is `preprocess_lowpass`, which
notebook 03 sets to 4 Hz; above ~4 Hz nothing can be explored from this
archive, and the full 100 Hz data stays in `data/raw/`.

## Known gaps, honestly

* **`scripts/get_data.py` has a placeholder URL.** Until you fill it in,
  nobody can fetch the data. This is the single blocking item.
* **The archive has not been downloaded through `get_data.py` end to end**,
  because the Zenodo record does not exist yet. The unpacking half is
  tested: the archive was extracted into a clean copy of the repo and all
  four notebooks ran green against it.
* **Single-station cross-components (EN, EZ, NZ) are not in the current
  `daily_acf.npz`.** The SDS was rebuilt with HHZ only to stay inside the
  disk budget, so the MSNoise run computed ZZ alone. The horizontals are in
  the raw archive; `scripts/build_sds.py --channels HHZ,HHN,HHE --hours all`
  restores them, at roughly 9 GB of SDS.
* **The MSNoise conda solve has not been run.** MSNoise was installed and
  exercised here with **pip** in a venv, and that is how every number above
  was produced. The conda-forge package is `msnoise 1.6.5 pyhd8ed1ab_1`
  (noarch) and its dependency `scandir` has cp312 binaries for linux-64,
  osx-64, osx-arm64 and win-64, so the solve should work — but should is not
  the same as did. **This is the thing to test first tomorrow**, and
  `environment-minimal.yml` is the fallback if it misbehaves.
* **No Windows or macOS has run this.** Everything was built and tested on
  Linux with Python 3.13 + obspy 1.5.1. The pinned environment asks for
  Python 3.12, which is the better-tested combination, but that exact
  combination has not been executed.
* **Colab has not been exercised.** The bootstrap cell is written and the
  logic is guarded, but it has not run in Colab.
* `EFIU`'s band amplitudes computed as counts/sensitivity differ from the
  fully deconvolved values by a flat 10 %, in every band, while its stated
  sensitivity matches the product of its stage gains exactly. The pipeline
  uses full deconvolution so nothing here depends on it, but if you ever use
  the scalar sensitivity at that station, know that it is off. `ECPN`,
  `ESLN` and `EMSG` agree to within 1 %. See
  `scripts/diagnostics/check_response_metadata.py`.

## Things that were found the hard way

These are all written up in the code and in `TEACHING.md`, but they are worth
knowing before you stand in front of anyone.

1. **The reference stack dominates everything.** Using the pre-eruption
   fortnight as reference drops median cc from 0.79 to 0.43 and usable days
   from 44 to 15. The grand mean is the default for that reason. This is
   now group C's exercise and it is the best thing in the practical.
2. **The autocorrelations genuinely decorrelate from their December shape
   starting 24–27 December**, the eruption onset, while adjacent days stay
   correlated at 0.85–0.96 throughout. So the measurement is not degrading;
   the waveform is really changing. Whether that is the medium or the tremor
   source is the honest open question, and 1–3 Hz is Etna's tremor band.
3. **Full spectral whitening destroys a single-station reflector** — a 4.00 s
   trough moves to 2.05 s — because the reverberations *are* the spectral
   comb. The default is running-mean whitening. `tests/test_acorr.py`
   measures both.
4. **`remove_response` costs 175 s per station-day at 100 Hz** and 21 s if
   you decimate first. For autocorrelation it is not needed at all, because
   the instrument phase cancels exactly for a single station.
5. **`tr.taper(0.02)` on a day-long trace tapers 29 minutes at each end**,
   which put a downward spike at every midnight in the band-amplitude series.
   Now a fixed 60 s.
6. **`layer_transfer` in the original `belgium-rf/hvsr.py` has the
   attenuation sign backwards** — `(1 + 1j/(2Q))` instead of `(1 - 1j/(2Q))`.
   The denominator passes through zero and the transfer function blows up
   with frequency, so overtones come out larger than the fundamental (A=50
   at the 7th overtone for a 30 m / 300 m-s layer). Fixed here. **Worth
   checking wherever else that function was copied**, since in the Belgium
   work it is a test function and a wrong test is worse than no test.
7. **MSNoise: `populate` must run BEFORE `scan_archive`.** `scan_archive`
   only walks directories for stations already in the Station table. Run it
   first and it prints "Scanning 0 directories", inserts nothing and **exits
   0** — it looks entirely healthy. This is the reverse of the order in most
   recipes, including `run_belgium.sh` (which gets away with it by using
   `populate --fromDA`).
8. **`cc_normalisation` does not exist in MSNoise 1.6.5** — it is master-era.
   Temporal normalisation there is `windsorizing`: 3.0 clips at 3σ, −1 is
   one-bit, 0 is none.
9. **MSNoise agrees with `etnalib/acorr.py` to 0.84–0.90** on the 1–15 s coda,
   and its output is even in lag to ~1e-7. Two independent implementations,
   different whitening and normalisation. Reproduce with
   `scripts/diagnostics/compare_msnoise_etnalib.py`.
10. **INGV dataselect runs at ~0.09 MB/s per connection**, ~0.7 MB/s with 12
   in parallel. The whole 3 GB download took 42 minutes. Students must never
   do this live.
