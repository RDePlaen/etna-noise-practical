# Demonstrator's guide — running the practical

One hour, bachelor students, most of whom have not used Python. This is the
run-of-show, the things that go wrong, and the answers to the questions in
the notebooks.

## The one-line version

> Ambient noise is a signal, not a nuisance; one station is enough to do real
> work with it; and the hard part is never the processing, it is deciding
> what your number means.

If they leave with only that, the hour worked.

---

## Timing

Realistically you have **~35 minutes of science**. Protect it.

| min | what | notes |
|---|---|---|
| 0–5 | everyone has a running notebook | see triage below. Do **not** debug one laptop while 25 people wait — move them to Colab and come back at the end |
| 5–15 | **01** — spectrogram, bands, find Cyclone Harry | the hook. Let them find the storm themselves |
| 15–30 | **02** — H/V at four stations | the synthetic check first, always |
| 30–50 | **03** — students run MSNoise themselves, then one knob per group | assign the knobs *before* they start; the pipeline takes ~45 s |
| 50–60 | **04** — collect the answers, argue | the board, not the screen |

If you are running late at minute 30, **cut notebook 02 to a demo you drive
from the front** and protect notebook 03. The dv/v exercise is the one with
the payoff, and it is the one that needs them to actually do something.

---

## Minute 0–5: triage, not debugging

Ask for a show of hands: "who got `obspy 1.5.1` printed out last night?"

* **Hands up:** they're fine.
* **Hands down:** send them straight to the Colab badge in the README. Do not
  investigate. A conda problem takes ten minutes and you do not have ten
  minutes.
* **Pairs:** if more than a third are struggling, pair them up at working
  laptops. This is better than it sounds — the discussion in notebook 03 works
  better in pairs anyway.

The failures you will actually see:

1. **"ModuleNotFoundError: obspy"** — they didn't run `conda activate
   etna-noise` in this terminal. Single most common issue by a wide margin.
2. **Colab shows "Restart session"** — click it, re-run the first cell. Normal.
3. **The conda solve never finished.** MSNoise drags in two downgrades
   (`sqlalchemy<2`, `setuptools<81`). Tell them to use
   `environment-minimal.yml` — identical practical, no MSNoise. Notebooks
   01–04 never import it.

---

## Notebook 01 — what to emphasise

The payoff is the multi-day spectrogram where Harry arrives. Let them look
before you explain.

**The question that makes it land:** *"an earthquake would be a vertical line
on this figure. Why is the storm a slow swelling instead?"* Because it is a
weather system hundreds of kilometres across, moving at weather speed. The
shape of the signal tells you the size of the source.

**The thing they get wrong:** that the band names are source names. They are
not. At `EFIU`, on the coast, the "tremor band" (0.5–5 Hz) rises by a factor
of 22 during the storm — that is not volcanic tremor, it is the sea and the
wind. At `ECPN`, on the summit, the same band is flat through the storm. The
band is a frequency range; what is in it depends on where you are standing.

Numbers from this dataset, if you want them on a slide:

| | 15 Jan (before) | 21 Jan (Harry) | factor |
|---|---|---|---|
| `EFIU` secondary microseism | 0.16 µm/s | 1.84 µm/s | **×12** |
| `ECPN` secondary microseism | 0.18 µm/s | 1.20 µm/s | ×6.7 |
| `ECPN` volcanic tremor | 2.4 µm/s | 1.9 µm/s | ×0.8 (flat) |

**Local hook:** the Cyclone Harry paper was led from Bari and Catania —
Scicchitano, Cannata and colleagues, *Communications Earth & Environment*,
2026. Ionian stations recorded the largest secondary-microseism amplitudes
since digital recording began in 1990, and a 16.6 m wave in the Strait of
Sicily. Several of them will have been in that storm.

---

## Notebook 02 — what to emphasise

**Always run the synthetic first.** It takes thirty seconds and it converts
H/V from a magic black box into a method with a known answer. The exact
transfer function peaks at Vs/4h and the code recovers the 30 m layer to
within 0.3 %.

**The planned surprise.** Ask them at the end of notebook 01 which station
will have the biggest H/V peak. Almost everyone says the coastal one. The
data says the opposite:

| station | f₀ | A₀ | SESAME |
|---|---|---|---|
| `ECPN` 3038 m | 0.89 Hz | 6.0 | 5/5 |
| `ESLN` 1787 m | 0.55 Hz | 4.0 | 5/5 |
| `EMSG` 1435 m | 1.48 Hz | 2.6 | 5/5 |
| `EFIU` 97 m (coast) | — | — | 0/5 |

The summit stations sit on hundreds of metres of fractured lava and loose
volcaniclastic material over denser rock — a perfectly good impedance
contrast. Being high up does not mean being on bedrock.

**Two caveats to raise yourself, because they won't:**

* `ECPN` sits on a steep summit cone, so some of that peak may be
  *topographic* amplification rather than a layer. On a volcano you cannot
  cleanly separate the two from H/V alone.
* `EFIU`'s only bump is at ~0.25 Hz, which is the **secondary microseism
  band** — a source artefact, not a site effect. This is exactly why the
  SESAME criteria exist. It is also a nice callback to notebook 01.

**The honest ending.** h = Vs/(4f₀) has two unknowns and H/V measures one.
Vs from 150 to 800 m/s gives depths differing by a factor of five. H/V tells
you *whether* there is a contrast and *roughly where*; it cannot give a depth
alone. Students find this unsatisfying, which is the correct response.

---

## Notebook 03 — the core, and how to run the exercise

**Say where the data came from, once.** The autocorrelations were computed
with MSNoise, the package the field actually uses, in single-station mode —
not by hand-rolled code. `etnalib/acorr.py` is there as an independent check
and notebook 03 shows the two agreeing (0.84–0.90 on the coda, and MSNoise's
output even in lag to 1e-7). Students should see that a real analysis gets
checked against something, and that two correct implementations do not give
identical numbers.

**Assign the knobs before anyone starts.** Walk round and give each pair a
letter. If you let them choose, half the room picks the same one.

| group | change | what they should find |
|---|---|---|
| A | `LAG_WIN = (1.0, 5.0)` | noisier, less sensitive — early coda has sampled less medium |
| B | `LAG_WIN = (10.0, 20.0)` | more sensitive but lower cc — late coda is better physics, worse data |
| C | `REF_DAYS = ("2025-12-10","2025-12-23")` | **the measurement collapses** — cc 0.79 → 0.43, usable days 44 → 15. The most important group |
| D | `STA = "ECPN"` | summit, noisier, tremor-contaminated |
| E | `STA = "EFIU"` | coast — dominated by the storm, not the volcano |
| F | rebuild in 0.5–1.5 Hz | deeper sensitivity, fewer windows |
| G | rebuild in 2–5 Hz | shallower, more scatter |

Groups F and G need a few minutes to rebuild the ACFs. Give them the letter
early, or skip F/G if you have fewer than six pairs.

**Collect three numbers per group on the board:** largest dv/v excursion,
when it happened, median correlation coefficient. Seven rows. Then the
conversation writes itself.

**What they should conclude:** everybody ran correct code on correct data and
got different numbers, because each made a different *defensible* choice. The
measurement is not "dv/v = x %". It is "dv/v = x %, in this band, in this
coda window, relative to this reference, with this cc". A paper that doesn't
state all five isn't reproducible.

Reference numbers, default settings (grand-mean reference, 1–3 Hz, 3–15 s):

| station | median cc | days cc>0.7 | dv/v std | range |
|---|---|---|---|---|
| `EMSG` | 0.79 | 44/72 | 0.11 % | −0.27 .. +0.23 % |
| `ESLN` | 0.75 | 42/64 | 0.20 % | −0.38 .. +0.35 % |
| `ECPN` | 0.72 | 31/61 | 0.14 % | −0.34 .. +0.17 % |
| `EFIU` | 0.70 | 37/72 | 0.16 % | −0.34 .. +0.46 % |

If a group reports ±5 %, the stretching has hit its search limit and their
correlation coefficient has collapsed — send them to look at the cc panel.

**Group C deserves five minutes of the closing discussion.** Using the quiet
pre-eruption fortnight as the reference is the choice everyone's instinct
says is correct, and it is the one that destroys the measurement. Two causes,
and it is worth making them separate them: (a) circularity — days inside the
reference window fit it partly because they built it, so the early part of
any such plot is flattered; and (b) something real — the autocorrelations
drift away from their December shape starting around 24–27 December, the
onset of the continuous eruptive phase, while *adjacent* days stay correlated
at 0.85–0.96 throughout. The measurement is not degrading. The waveform is
genuinely changing and staying changed.

Whether that is the medium or the noise source is the honest open question,
and it leads straight into notebook 04: Etna's tremor sits at 0.5–5 Hz, which
is the band being correlated, and it radiates from a fixed point at the
summit — exactly what the diffuse-field assumption forbids.

---

## Notebook 04 — the argument

Do this on the board, not on screens.

The structure to hold them to:

1. Which panel responds to which event? (Tremor → volcanic dates. Microseism
   → storm dates. Check it rather than assuming it.)
2. Does dv/v move in the storm window?
3. If it does, **three** explanations compete: the medium really changed
   (loading, pore pressure); the noise *source* moved, so the
   autocorrelation is no longer sampling the same way; or the measurement
   just degraded.
4. How would you tell them apart with the data you already have? — does the
   cc drop at the same time (→ 2 or 3)? Is the excursion coherent across all
   four stations (→ 1) or not (→ 2)? Does it recover instantly (elastic
   loading) or over weeks (damage healing)?

**Do not resolve it for them.** This is a genuinely open kind of question and
the correct student response is "we'd need the meteorological data and GNSS
to be sure". If someone says that, you have won.

The closing point: detection is easy, **attribution is a separate and harder
problem**, and that is why volcano observatories are multi-parameter — and
why the tremor and the dv/v are analysed in different rooms of the same
building at INGV-OE down the road.

---

## Questions you will be asked

**"Why does the autocorrelation start at 1?"**
Zero lag is the trace correlated with itself with no shift — identical by
construction, in every trace, on every day. It carries no information, which
is why the measurement window never includes it.

**"Why is the reflector a trough and not a peak?"**
Claerbout (1968): the autocorrelation of the transmission response equals
*minus* the reflection response at non-zero lag. An impedance increase gives
a negative peak. `python tests/test_acorr.py` plants r = 0.4 and recovers
−0.394.

**"Why not just use earthquakes?"**
You have to wait for them, they come from wherever they come from, and no two
are in the same place. Noise is continuous, free, and — this is the point —
*repeatable*, so you can compare today with yesterday.

**"Is 0.1 % a real change?"**
Only if the correlation coefficient supports it and it is coherent across
stations. Hence the cc panel underneath every dv/v plot.

**"Why don't we remove the instrument response for the autocorrelations?"**
Because for a *single station* it cancels exactly: C(f) = H*(f)H(f)S(f) =
|H(f)|²S(f), so the phase — which is what sets travel time — drops out. It
does *not* cancel for cross-correlation between two stations with different
sensors, and there you must remove it. (It also costs 175 s per station-day
at 100 Hz, which is its own kind of argument.)

---

## If you have spare time

* `python tests/test_acorr.py` live on the projector. Full spectral whitening
  moves a 4.00 s reflector to 2.05 s; running-mean whitening leaves it at
  4.00 s. A one-line parameter change that silently destroys the result is a
  better lesson about scientific software than any lecture about it.
* The `15` OGS microzonation network — 18 stations in 18 municipalities
  around Etna, open on FDSN, purpose-built for site response. Notebook 02
  scaled up to a real hazard map. Good thesis topic, say so out loud.
* The `4S` DAS dataset at Serra La Nave: a fibre in a 200 m borehole at the
  same site as `ESLN`, thousands of channels instead of one.
