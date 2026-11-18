import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nbtools import md, code, build, BOOTSTRAP, COLAB_BADGE

c = []
c.append(md(f"""# 01 — What *is* ambient seismic noise?

{COLAB_BADGE.format(name="01_what_is_ambient_noise.ipynb")}

**~10 minutes.** A seismometer never records nothing. Between the
earthquakes there is a continuous signal, and for a long time it was treated
as a nuisance to be averaged away. This notebook is about what that signal
actually contains — and the answer is: the ocean, the weather, human
activity, and the volcano, each in its own frequency band.

By the end you should be able to look at a spectrogram of an Etna station
and say what each horizontal stripe is."""))

c.append(code(BOOTSTRAP))

c.append(md("""## 1. One day, one station, looked at three ways

Start with the coastal station `EFIU`, 3 km from the Ionian Sea, on
**20 January 2026**. Hold on to that date."""))

c.append(code('''from etnalib import noise

STA, DAY = "EFIU", "2026-01-20"
# Raw waveforms ship for 2026-01-15 and 2026-01-20 at every station
# (ESLN was down on the 20th) and for the whole fortnight at EFIU.
# Everything else in this notebook is precomputed for all four stations.
# If you change STA or DAY, require_day tells you what is available.
print("raw waveform days:",
      {s: len(io.available_days(s)) for s in C.STATIONS})

tr = io.require_day(STA, DAY)
print(tr)
print(f"coverage of the day: {tr.stats.coverage:.0%}")

fig, ax = plt.subplots(figsize=(11, 2.4))
t = np.arange(tr.stats.npts) / tr.stats.sampling_rate / 3600
ax.plot(t, tr.data * 1e6, lw=0.3, color="#333")
ax.set_xlabel("hours UTC on " + DAY); ax.set_ylabel("velocity  [µm/s]")
ax.set_title(f"{STA} — {DAY} — vertical ground velocity"); ax.set_xlim(0, 24)
plt.show()'''))

c.append(md("""The trace is not flat and it is not an earthquake. It wobbles continuously,
and the amplitude changes slowly through the day. To see *why*, we need
frequency."""))

c.append(code('''from scipy.signal import spectrogram

f, tt, S = spectrogram(tr.data, fs=tr.stats.sampling_rate,
                       nperseg=4096, noverlap=2048, scaling="density")
keep = (f > 0.02) & (f < 9)

fig, ax = plt.subplots(figsize=(11, 4))
m = ax.pcolormesh(tt / 3600, f[keep], 10 * np.log10(S[keep]),
                  shading="auto", cmap="magma", vmin=-180, vmax=-110)
ax.set_yscale("log"); ax.set_ylabel("frequency  [Hz]")
ax.set_xlabel("hours UTC on " + DAY); ax.set_xlim(0, 24)
plt.colorbar(m, ax=ax, label="power  [dB rel. 1 (m/s)²/Hz]", pad=0.01)

for name, (f1, f2) in C.BANDS.items():
    if f2 > 9: continue
    ax.axhline(f1, color="w", lw=0.5, alpha=0.4)
    ax.text(0.3, np.sqrt(f1 * f2), name, color="w", fontsize=8, va="center")
ax.set_title(f"{STA} — {DAY}")
plt.show()'''))

c.append(md("""### What you are looking at

Reading from the bottom up:

* **0.05–0.1 Hz — the primary microseism (10–20 s).** Ocean swell reaching
  shallow water and pushing directly on the sea floor. Its frequency is the
  frequency of the ocean waves themselves.
* **0.1–0.35 Hz — the secondary microseism (3–10 s).** The loudest signal on
  almost every seismometer on Earth. Two trains of ocean waves travelling in
  *opposite* directions interfere and produce a standing pressure
  fluctuation at **twice** the wave frequency — which is why this band sits
  at double the primary. Longuet-Higgins worked this out in 1950.
* **0.5–5 Hz — at Etna, volcanic tremor.** Continuous, sustained, radiating
  from the summit craters. Not an earthquake: it has no onset and no end.
* **above ~5 Hz — human beings.** Traffic, industry, machinery.

**Question for you:** which of these four should change between day and
night, and which should not? Look back at the spectrogram before reading on."""))

c.append(md("""## 2. The daily cycle, and what it tells you about sources

Compress the spectrogram into one number per band per half-hour: the RMS
ground velocity in that band. This is exactly the quantity behind **RSAM**
in volcano monitoring and behind the microseism amplitude used in storm
studies — same arithmetic, different band, two different literatures."""))

c.append(code('''fig, ax = plt.subplots(figsize=(11, 3.6))
for name, band in C.BANDS.items():
    if band[1] >= 10: continue
    tt_, a = noise.band_amplitude(tr, band)
    hrs = (tt_ - tt_[0]) / np.timedelta64(1, "h")
    ax.semilogy(hrs, a * 1e6, label=f"{name}  ({band[0]}–{band[1]} Hz)", lw=1.6)
ax.set_xlabel("hours UTC"); ax.set_ylabel("RMS velocity  [µm/s]")
ax.set_xlim(0, 24); ax.legend(fontsize=8, ncol=3)
ax.set_title(f"{STA} — {DAY} — band amplitudes")
plt.show()'''))

c.append(md("""The anthropogenic band has a clear daily cycle; the microseism bands do not.
That is the signature of the source: people sleep, the ocean does not.

> **This is why notebook 00 downloads only 00:00–06:00 UTC on most days.**
> Comparing the same night-time hours every day removes the one thing that
> is guaranteed to change for reasons that have nothing to do with the
> ground."""))

c.append(md("""## 3. Now look at a week — and meet Cyclone Harry

Between **16 and 23 January 2026**, an extratropical cyclone sat over the
Ionian Sea. Ionian stations recorded the largest secondary-microseism
amplitudes **since digital recording began in 1990**, and the highest
individual wave ever measured in the Mediterranean (16.6 m, in the Strait of
Sicily on 21 January). The work was led from Bari and **Catania** —
Scicchitano, Cannata and colleagues, *Communications Earth & Environment*,
2026.

You have the seismic data. Go and find the storm."""))

c.append(code('''days = [f"2026-01-{d:02d}" for d in range(14, 27)]
days = [d for d in days if d in io.available_days(STA)]
print(len(days), "days available")

tt_, f_, S_ = noise.long_spectrogram(STA, days, nperseg=4096, noverlap=2048,
                                     fmax=5.0)
fig, ax = plt.subplots(figsize=(12, 4.2))
m = ax.pcolormesh(tt_, f_[f_ > 0.03], 10 * np.log10(S_[f_ > 0.03]),
                  shading="auto", cmap="magma", vmin=-175, vmax=-105)
ax.set_yscale("log"); ax.set_ylabel("frequency [Hz]")
plots.date_axis(ax); plots.mark_events(ax)
plt.colorbar(m, ax=ax, label="power [dB]", pad=0.01)
ax.set_title(f"{STA} — the arrival of Cyclone Harry")
plt.show()'''))

c.append(md("""**Read the figure.** The secondary microseism band brightens by orders of
magnitude, over days rather than seconds, and fades again. An earthquake
would be a vertical line. A storm is a slow swelling, because it is a
weather system, hundreds of kilometres across, moving at weather speed.

Note also what *does not* change: the tremor band. Different source,
different physics."""))

c.append(md("""## 4. The same storm, seen from four stations

`EFIU` is 3 km from the coast; `ECPN` is on the summit, 22 km inland and
3 km above sea level. If the secondary microseism really comes from the sea,
the two should not record it equally."""))

c.append(code('''# Precomputed, because deconvolving four stations x three months live would
# take an hour. The recipe is noise.band_series() -- read it if you want.
B = np.load(C.PROC / "band_series.npz")

fig, ax = plt.subplots(figsize=(12, 3.8))
for sta in C.STATIONS:
    tk, vk = f"{sta}_secondary_t", f"{sta}_secondary_v"
    if tk not in B: continue
    t_, a_ = B[tk].astype("datetime64[s]"), B[vk]
    w = (t_ >= np.datetime64("2026-01-13")) & (t_ <= np.datetime64("2026-01-28"))
    if w.sum() == 0: continue
    tg, ag = plots.break_gaps(t_[w], a_[w] * 1e6)
    ax.semilogy(tg, ag, lw=1.2, color=plots.station_color(sta),
                label=f"{sta} ({C.STATIONS[sta]['elev']:.0f} m, "
                      f"{C.STATIONS[sta]['coast_km']:.0f} km from coast)")
ax.set_ylabel("secondary microseism\\nRMS velocity [µm/s]")
plots.date_axis(ax); plots.mark_events(ax, kinds=("storm",))
ax.legend(fontsize=8); ax.set_title("0.1–0.35 Hz at four Etna stations")
plt.show()'''))

c.append(md("""**Discuss.** The curves move together — the storm is a regional source, so
every station sees the same event at the same time — but they are offset
from each other. What sets the offset: distance to the coast, elevation, or
the rock under the station? You cannot separate those three with four
stations, which is itself worth knowing.

### Where this goes next
This is the measurement behind an active collaboration between Catania and
Bari: seismic stations as a **permanent, all-weather sea-state sensor**.
Buoys break in big storms, and the biggest storms are the ones you most want
to measure. Seismometers do not break, and Italy has hundreds of them that
have been running since the 1990s."""))

c.append(md("""## 5. The PPSD: is this station *normally* like this?

A spectrogram shows what happened. A **probability power spectral density**
(McNamara & Buland 2004) throws time away and keeps the statistics: for each
frequency, how often is the power at each level? It is how you tell "this
station is noisy today" from "this station is broken"."""))

c.append(code('''import pickle
with open(C.PROC / "ppsd_summary.pkl", "rb") as fh:
    P = pickle.load(fh)

stas = [s for s in C.STATIONS if s in P]
vmax = plots.ppsd_vmax([P[s]["hist"] for s in stas])

fig, axes = plt.subplots(1, len(stas), figsize=(3.3 * len(stas), 4),
                         sharey=True, sharex=True)
for ax, sta in zip(np.atleast_1d(axes), stas):
    m = plots.ppsd_panel(ax, P[sta], vmax)
    ax.set_title(f"{sta} — {plots.ppsd_hours(P[sta]):.0f} h", fontsize=10)
np.atleast_1d(axes)[0].set_ylabel("power [dB rel. 1 (m/s)²/Hz]")
plt.colorbar(m, ax=list(np.atleast_1d(axes)), label="probability", pad=0.01)
plt.show()'''))

c.append(md("""The dashed lines are the Peterson **new low- and high-noise models** — the
quietest and noisiest sites ever recorded, the envelope that most stations
on Earth live inside. A station far outside them is usually a metadata
error rather than a discovery.

The dotted lines mark the two bands from the spectrogram. The hump at
3–10 s is the secondary microseism, and it is there in the *statistics*,
not only on stormy days — it is the permanent background of the planet.

**Two things to notice, and one of them is a trap.**

All four stations sit **above the high-noise model** at periods shorter than
about 1 s, `ECPN` most of all. On a volcano in eruption that is not a fault:
the noise models describe ordinary sites, and continuous volcanic tremor is
not ordinary. Knowing when being "off the scale" is real and when it means a
broken response is a genuinely useful instinct.

Compare the four panels at short period and you have already ranked the
stations by how much tremor they see — summit loudest, coast quietest. That
ranking is what notebook 03 exploits."""))

c.append(md("""---

## Before you move on

1. Which band would you monitor to detect a **storm**? To detect an
   **eruption**? Which band would you avoid entirely, and why?
2. The secondary microseism sits at twice the ocean-wave frequency. If a
   swell has a 10 s period, where is its secondary microseism?
3. `ECPN` sits at 3038 m on a volcano and `EFIU` at 97 m on the coastal
   plain. Before notebook 02: which do you expect to amplify ground motion
   more, and why?"""))

build(c, "notebooks/01_what_is_ambient_noise.ipynb")
print("wrote 01")
