# Ambient seismic noise at Etna — seminar practical

A one-hour practical on ambient-noise methods for **monitoring** and
**imaging**, built around three months of open data from Mount Etna.

**New here? Go to [SETUP.md](SETUP.md) first.** It takes about 20 minutes and
you must do it before the session.

---

## The idea

Between the earthquakes, a seismometer records a continuous signal made by
the ocean, the weather, human activity and — on a volcano — the volcano
itself. Twenty years ago that signal was a nuisance. Today it is used to
image the subsurface and to monitor it continuously, without ever waiting for
an earthquake.

This practical uses a window of data, **10 December 2025 to 20 February
2026**, in which three different things happen at once:

| when | what |
|---|---|
| 24 Dec 2025 | Etna enters a continuous eruptive phase |
| 1 Jan 2026 | a fissure opens at ~2100 m in the Valle del Bove |
| **16–23 Jan 2026** | **Cyclone Harry** — record secondary-microseism amplitudes in the Ionian |
| 6 Feb 2026 | a new eruption begins at the NE base of the SE Crater |

So a change measured in January has at least two candidate causes, one
volcanic and one meteorological. Telling them apart is the point.

## The notebooks

| notebook | time | what you do |
|---|---|---|
| [`01_what_is_ambient_noise`](notebooks/01_what_is_ambient_noise.ipynb) | ~10 min | spectrograms, band amplitudes, PPSD. Find Cyclone Harry in the data. |
| [`02_hvsr`](notebooks/02_hvsr.ipynb) | ~15 min | H/V spectral ratio at four stations; pick f₀; turn it into a depth. |
| [`03_autocorrelation_monitoring`](notebooks/03_autocorrelation_monitoring.ipynb) | ~20 min | daily autocorrelations → dv/v. Each group changes **one** parameter. |
| [`04_synthesis`](notebooks/04_synthesis.ipynb) | ~10 min | put dv/v, tremor and microseism on one axis and argue. |
| [`00_download_data`](notebooks/00_download_data.ipynb) | **2 h — not in class** | rebuild the whole archive from the FDSN server. |

## The data

Four stations of the Italian National Seismic Network (`IV`), operated by
INGV, distributed openly through EIDA and FDSN.

| station | elevation | distance to coast | role |
|---|---|---|---|
| `ECPN` | 3038 m | 22 km | summit, Cratere del Piano — volcanic tremor |
| `ESLN` | 1787 m | 23 km | Serra La Nave — the classic Etna reference site |
| `EMSG` | 1435 m | 29 km | NW flank — the cleanest continuous record in this window |
| `EFIU` | 97 m | 3 km | Fiumefreddo di Sicilia — strongest microseism, softest ground |

**All 24 hours of every day**, three components, ~14 GB raw. The first
version of this dataset kept only 00:00–06:00 UTC on most days — the usual
monitoring choice, since night hours are quieter and using the same hours
each day removes the daily cycle. Whole days are used instead because the
diurnal cycle is itself part of the discussion: in the 0.5–5 Hz "tremor"
band it is about ±20 % even at the summit station.

`etnalib/config.py` has `ACORR_HOURS` — `None` for whole days, `(0, 6)` to
reproduce the night-only results. Notebook 02 compares them.

**Citation.** INGV Seismological Data Centre (1997), *Rete Sismica Nazionale
(RSN)*, Istituto Nazionale di Geofisica e Vulcanologia (INGV), Italy,
[doi:10.13127/sd/x0fxnh7qfy](https://doi.org/10.13127/sd/x0fxnh7qfy).

## What it looks like

![secondary microseism at four stations](figures/nb01_microseism_stations.png)

Cyclone Harry arriving, 13–27 January 2026. All four stations rise together
by a factor of ~20 — a regional source — with the coastal station highest and
the summit station lowest throughout. The sharp spike on 18 January is an
earthquake, visible on every station: a different kind of signal entirely.

![dv/v at ECPN](figures/nb03_dvv_ECPN.png)

Relative velocity change from MSNoise daily autocorrelations, with the
correlation coefficient underneath. Grey crosses are days we do not trust.
Measured performance, one-bit normalisation, 1–3 Hz, coda 3–15 s,
grand-mean reference:

| station | median cc | days cc>0.7 | dv/v std |
|---|---|---|---|
| `EFIU` | 0.78 | 49/72 | 0.66 % |
| `ECPN` | 0.73 | 31/61 | 0.12 % |
| `EMSG` | 0.67 | 27/72 | 0.28 % |
| `ESLN` | 0.65 | 21/64 | 0.33 % |

## MSNoise version

Notebook 03 runs **MSNoise 1.6.5** live, on a 57 MB single-station archive
shipped in the data download (`data/sds_teaching/`). The whole pipeline —
`compute_cc`, `stack -r`, `stack -m`, `compute_mwcs`, `compute_dtt` and four
plots — takes about **45 seconds**.

> Most MSNoise tutorials online, including Skience 2024, use the development
> version, where the commands are nested (`msnoise cc compute_cc`, `msnoise
> cc dvv compute_mwcs`). In 1.6.5 they are flat (`msnoise compute_cc`,
> `msnoise compute_mwcs`), and `mov_stack` is an integer number of days
> rather than a tuple. The notebook uses the 1.6.5 forms throughout and says
> so, because copying a command from a newer tutorial is the most likely way
> to get stuck.

Defaults worth knowing: `dtt_lag=static`, `dtt_minlag=5.0`, `dtt_width=30.0`
give a **5–35 s** lag window, which is exactly what De Plaen et al. (2019)
used at Etna.

## How the correlations are computed

The autocorrelations (`ZZ`) and the single-station cross-components
(`ZN`, `ZE`, `NE`) are computed with **MSNoise 1.6.5** in its single-station
mode — `components_to_compute_single_station` — and exported to a `.npz` the
notebooks read. `scripts/run_msnoise.sh` does the whole chain; notebook 00
documents it.

`etnalib/acorr.py` is a second, independent implementation of the same thing
in ~40 lines. It exists as a **check**, and notebook 03 shows the comparison:
MSNoise's output is even in lag to one part in 10⁷, as an autocorrelation
must be, and correlates 0.84–0.90 with `etnalib` over the 1–15 s coda despite
different whitening and normalisation.

> MSNoise 1.6.5 needs **`sqlalchemy<2`** and **`setuptools<81`**. Both are
> pinned in `environment.yml` with the exact failure they prevent written
> next to them, and `run_msnoise.sh` refuses to start without them. Neither
> failure is visible until the stage that hits it.

## The code

`etnalib/` is a small library, not a black box — read it.

| module | what |
|---|---|
| `config.py` | stations, window, bands, chronology. Change things here. |
| `download.py` | chunked, restartable, parallel FDSN download |
| `io.py` | one station-day → one clean, response-corrected trace |
| `noise.py` | spectrograms, PPSD, band amplitudes (RSAM / microseism) |
| `hvsr.py` | Nakamura H/V with Konno–Ohmachi smoothing and SESAME criteria |
| `acorr.py` | autocorrelation: CCGN and PCC, linear and phase-weighted stacks. Used as an independent **check** on MSNoise, not as the production path |
| `dvv.py` | stretching, for the relative velocity change |
| `plots.py` | shared styling, and the event chronology on every time axis |

### Nothing here is trusted because it ran

```
python tests/test_acorr.py     # Claerbout sign, whitening, PCC amplitude bias
python -m etnalib.dvv          # plant a known dv/v and recover it
```

Two results from those tests are worth knowing before you start:

* **An impedance increase appears as a *negative* peak** in the
  autocorrelation, at the two-way travel time (Claerbout 1968). The test
  plants a reflector with r = 0.4 and recovers a trough of −0.394.
* **The reference stack is the most consequential parameter**, not the band
  and not the coda window. Taking the quiet fortnight before the eruption as
  the reference — the choice everyone's instinct prefers — drops the median
  correlation coefficient from 0.79 to 0.43 and the usable days from 44 to
  15. Taking the mean of every day fixes it. Notebook 03 makes one group
  discover this.
* **Full spectral whitening destroys the thing you are imaging.** The
  reverberations of a layer *are* a comb in the spectrum, and flattening the
  spectrum erases them — the test measures a 4.00 s reflector moving to
  2.05 s. The default is therefore a running-mean whitening, which removes
  the colour of the noise while keeping the structure.

## Further reading

**Methods**
Nakamura (1989), QR RTRI 30, 25–33 — H/V.
SESAME (2004), *European guidelines for H/V*, deliverable D23.12.
Konno & Ohmachi (1998), BSSA 88, 228–241 — the smoothing.
Claerbout (1968), Geophysics 33, 264–269 — autocorrelation and reflectivity.
Sens-Schönfelder & Wegler (2006), GRL 33 — stretching, and dv/v at a volcano.
Schimmel (1999), IEEE TGRS 37 — phase cross-correlation.
Romero & Schimmel (2018), JGR 123, 5052–5067 — PCC autocorrelation imaging.
McNamara & Buland (2004), BSSA 94, 1517–1527 — PPSD.
Longuet-Higgins (1950), Phil. Trans. R. Soc. A 243 — why the secondary
microseism is at twice the ocean-wave frequency.

**This dataset**
Scicchitano, Cannata et al. (2026), *Unprecedented storm impact in the
Mediterranean revealed by models and seismic noise during Cyclone Harry*,
Communications Earth & Environment.
Cannata et al. (2022), Sci. Rep. 12 — microseism monitoring of Medicane Apollo.
INGV Osservatorio Etneo weekly bulletins, and the Smithsonian Global
Volcanism Program weekly reports, for the eruptive chronology.

**Other open Etna data, if you want to go further**
`4S` (2022–2023) — Serra La Nave distributed acoustic sensing, 200 m borehole.
`9N` (2018) — DAS, seismological and infrasonic data (Jousset et al. 2022).
`7I` (2025) — one week of rotational data on the same plate as `ESLN`.
`15` (2022–2024) — OGS microzonation array, 18 short-period stations in 18
municipalities around Etna. Purpose-built for site-effect work: the obvious
extension of notebook 02.
