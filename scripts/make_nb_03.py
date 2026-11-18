import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nbtools import md, code, build, COLAB_BADGE

c = []
c.append(md(f"""# 03 — Monitoring with MSNoise

{COLAB_BADGE.format(name="03_autocorrelation_monitoring.ipynb")}

**~20 minutes. You run the whole pipeline yourself.**

Notebook 02 asked *what is under the station*. This one asks **is it
changing?** — and does it with [MSNoise](http://www.msnoise.org), the package
the ambient-noise community actually uses, rather than with code written for
this class.

We work on **`ECPN`**, the summit station, in the **1–2 Hz** band — the
station and band of De Plaen et al. (2019), who did this at Etna with the
same method ([Front. Earth Sci. 6:251](https://doi.org/10.3389/feart.2018.00251)).

### Autocorrelation, in one paragraph

Cross-correlate the noise at two stations and you get the wave that would
travel between them. Correlate a station **with itself** and you get the
wave that goes down, reflects, and comes back — the zero-offset reflection
response (Claerbout, 1968). **One station, no pairs.** MSNoise does this
through `components_to_compute_single_station`.

> **Version note.** This runs **MSNoise 1.6.5**. Most tutorials you will find
> online — including Skience 2024 — use the development version, whose
> commands are nested: `msnoise cc compute_cc`, `msnoise cc dvv
> compute_mwcs`. In 1.6.5 they are flat: `msnoise compute_cc`, `msnoise
> compute_mwcs`. If you copy a command from elsewhere and it says "no such
> command", that is why."""))

c.append(code('''import os, sys, shutil, subprocess
from pathlib import Path
%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.abspath("."))
sys.path.insert(0, os.path.abspath(".."))
from etnalib import config as C

REPO = Path(C.ROOT)
SDS = REPO / "data" / "sds_teaching"
PROJ = REPO / "msnoise_live"

assert shutil.which("msnoise"), (
    "msnoise is not on PATH. You are probably in the fallback environment "
    "(etna-noise-minimal). Use the prepared results in notebook 04, or "
    "create the full environment from environment.yml.")
assert SDS.exists(), f"{SDS} missing -- run scripts/get_data.py"
SETUP = REPO / "scripts" / "msnoise_setup.py"
assert SETUP.exists(), (
    f"{SETUP} is missing. The repository is incomplete or out of date -- this "
    "notebook needs the MSNoise version of the code stack. Extract the code "
    "archive into an EMPTY directory and move your data/ folder into it; do "
    "not unpack it on top of an older copy. If you see a nested "
    f"{REPO.name}/{REPO.name}/ directory, that is what happened.")
print("msnoise:", subprocess.run(["msnoise", "--version"], capture_output=True,
                                 text=True).stdout.strip() or "ok")
print("archive:", len(list(SDS.rglob("*.D.*"))), "day files")'''))

c.append(md("""## 1. Create the project

MSNoise keeps its state in a small database — which jobs exist, which have
run, what the parameters are. `db init --tech 1` makes an SQLite one in the
current directory.

Two settings do the real work here:

* `components_to_compute` is **empty** — we want nothing between stations
* `components_to_compute_single_station = ZZ` — the autocorrelation of the
  vertical"""))

c.append(code('''if PROJ.exists():
    shutil.rmtree(PROJ)          # start clean every time you run this
PROJ.mkdir(parents=True)
%cd {PROJ}

!msnoise db init --tech 1 2>&1 | tail -1
# subprocess with check, not "!", so a failure here stops the notebook
# instead of printing an error and leaving an empty project behind.
r = subprocess.run([sys.executable, str(SETUP), "--sds", str(SDS),
                    "--response", str(REPO / "data" / "processed"),
                    "--components", "ZZ", "--windsorizing=-1",
                    "--fmin", "1.0", "--fmax", "2.0"],
                   capture_output=True, text=True)
print("\\n".join(r.stdout.strip().splitlines()[-4:]))
if r.returncode:
    raise SystemExit("msnoise_setup.py failed:\\n" + r.stderr.strip()[-2000:])
!msnoise config set maxlag=40 2>&1 | tail -1
!msnoise config set mov_stack=2 2>&1 | tail -1
# The teaching archive is 10 Hz, so Nyquist is 5 Hz and preprocess_lowpass
# must come down from its 8 Hz default. Everything we do is in 1-2 Hz.
!msnoise config set cc_sampling_rate=10 2>&1 | tail -1
!msnoise config set preprocess_lowpass=4.0 2>&1 | tail -1'''))

c.append(md("""### The parameters, and where they come from

| setting | value | why |
|---|---|---|
| `cc_sampling_rate` | 10 Hz | the archive is already at 10 Hz, so resampling is a no-op |
| `preprocess_lowpass` | 4 Hz | must stay below the 5 Hz Nyquist |
| `maxlag` | 40 s | the measurement needs lags out to 35 s |
| `corr_duration` | 1800 s | 30-minute windows, stacked to a day |
| `windsorizing` | −1 | **one-bit**. Not the default, and it matters — see below |
| `mov_stack` | 2 | 2-day stacks, as in De Plaen et al. (2019) |
| `remove_response` | N | for one station the instrument phase cancels exactly |
| filter 1 | 1–2 Hz | the paper's band |

**On `windsorizing`.** It sets the temporal normalisation: `3.0` clips each
window at 3σ (the usual default), `-1` is one-bit, `0` is none. Etna hosts
continuous volcanic tremor at 0.5–5 Hz radiating from a fixed point at the
summit craters — the opposite of the diffuse, all-directions noise field the
method assumes. Aggressive normalisation is the standard defence, and on
this dataset it took the number of usable days from 21 to 128 out of 269.

**On 10 Hz.** The archive is decimated to 10 Hz, which halves its size. For
a 1–2 Hz measurement that costs nothing measurable: running this identical
chain at 20 Hz and at 10 Hz on the same station gives dv/v curves correlated
at **0.965**, with a median difference of **0.018 %** against a signal
standard deviation of 0.083 % and MWCS errors near 0.13 %. The sampling rate
is not what limits this. It does mean no band above ~4 Hz can be explored
from this archive — the full 100 Hz data is in `data/raw/` if you want it.

**On the lag window.** MSNoise's defaults are `dtt_lag=static`,
`dtt_minlag=5.0`, `dtt_width=30.0` — a **5–35 s** window, which happens to
be exactly what the paper used."""))

c.append(md("""## 2. Scan the archive

`populate` reads the station metadata; `scan_archive` walks the waveform
files and records what exists; `new_jobs` turns that into a to-do list.

> **Order matters and nothing warns you.** `scan_archive` only looks for
> stations that are *already* in the station table. Run it before
> `populate` and it reports "Scanning 0 directories", inserts nothing, and
> exits successfully."""))

c.append(code('''!msnoise populate 2>&1 | tail -2
!msnoise scan_archive --init 2>&1 | tail -1
!msnoise new_jobs --init 2>&1 | tail -1
!msnoise info -j 2>&1 | tail -6'''))

c.append(code('''from msnoise.plots.data_availability import main as plot_da
plot_da(show=True)'''))

c.append(md("""## 3. Compute the correlations

This is the expensive step — a couple of minutes. Each day is cut into
30-minute windows, each window is one-bit normalised and autocorrelated, and
the windows are stacked into one autocorrelation per day."""))

c.append(code('''!msnoise -t 2 compute_cc 2>&1 | grep -Ei "finished|error" | tail -2
!msnoise info -j 2>&1 | tail -5'''))

c.append(md("""## 4. Look at what you just made — the CCF-time plot

This is the plot to spend time on. Every horizontal trace is one day's
autocorrelation; time runs up the page, lag runs across it."""))

c.append(code('''from msnoise.plots.ccftime import main as plot_ccftime
plot_ccftime("IV.ECPN", "IV.ECPN", filterid=1, components="ZZ",
             mov_stack=1, ampli=8, show=True)'''))

c.append(md("""**Read it.**

* **It is symmetric about zero lag.** It has to be: an autocorrelation is an
  even function, because correlating a signal with itself cannot distinguish
  past from future. If yours is not symmetric, something is wrong.
* **The zero-lag spike is identical every day** and carries no information.
  Never measure anything on it — which is why the lag window starts at 5 s.
* **The coda repeats, day after day, for months.** That repeatability is the
  entire basis of the method. Waves that have bounced around for 10 seconds
  have sampled the medium far more thoroughly than a direct arrival, so a
  *tiny* velocity change produces a *measurable* shift out there.

Try `envelope=True`, or `refilter="1:1.5"`, and see what changes."""))

c.append(md("""## 5. Stack: the reference, then the moving stacks

dv/v is always measured **relative to something**. `stack -r` builds the
reference from the whole period; `stack -m` builds the 2-day moving stacks
that get compared against it.

The `reset STACK` in between is not optional — MSNoise deliberately leaves
the stack jobs "in progress" after the reference so you can re-run them for
the moving stacks."""))

c.append(code('''!msnoise stack -r 2>&1 | tail -2
!msnoise reset STACK 2>&1 | tail -1
!msnoise -t 2 stack -m 2>&1 | grep -Ei "finished|error" | tail -1'''))

c.append(code('''from msnoise.plots.interferogram import main as plot_ifg
plot_ifg("IV.ECPN", "IV.ECPN", filterid=1, components="ZZ",
         mov_stack=2, show=True)'''))

c.append(md("""The interferogram is the same data as an image rather than as wiggles.
Vertical stripes that stay put mean a stable medium; a stripe that leans
means the arrival is moving in time, which is what dv/v measures."""))

c.append(md("""## 6. Measure: MWCS, then dt/t

**MWCS** (moving-window cross-spectral analysis) slides a short window along
the coda and, in each one, measures the time shift `dt` between this stack
and the reference from the slope of the cross-spectral phase.

**`compute_dtt`** then regresses those `dt` against lag time `t`. For a
velocity change spread uniformly through the medium that regression is a
straight line through the origin, and

$$\\frac{dv}{v} = -\\frac{dt}{t}$$

A point is thrown away if its error exceeds `dtt_maxerr` (0.1 s) or its
coherence is below `dtt_mincoh`."""))

c.append(code('''!msnoise config set dtt_mincoh=0.6 2>&1 | tail -1
!msnoise -t 2 compute_mwcs 2>&1 | grep -Ei "finished|error" | tail -1
!msnoise compute_dtt 2>&1 | grep -Ei "finished|error" | tail -1'''))

c.append(code('''from msnoise.plots.mwcs import main as plot_mwcs
# mov_stack=2 here, not 1: MWCS and DTT are written per moving stack, and
# we configured 2-day stacks. Asking for mov_stack=1 finds no files and
# fails with a bare IndexError rather than a useful message.
plot_mwcs("IV.ECPN", "IV.ECPN", filterid=1, components="ZZ",
          mov_stack=2, show=True)'''))

c.append(md("""Each panel is one day: `dt` against lag time, with the fitted line. **This
is the actual measurement** — everything else is bookkeeping. A good day is
a tidy straight line through the origin. A day where the points scatter, or
where the line is straight but does not pass through the origin, is telling
you the simple picture does not hold."""))

c.append(code('''from msnoise.plots.dvv import main as plot_dvv
plot_dvv(mov_stack=2, dttname="M", components="ZZ", filterid=1, show=True)'''))

c.append(md("""## 7. ⚠️ The exercise: one knob each

**Work in pairs. Your demonstrator will give each pair ONE parameter.**
Change only that one, re-run from the `compute_dtt` cell (nothing before it
needs redoing — except group F, who must redo `stack -m` and
`compute_mwcs`), and write your answer on the board.

Before re-running `compute_dtt`, clear the old result: `!msnoise reset DTT`.

| group | command | what it changes |
|---|---|---|
| A | `!msnoise config set dtt_minlag=10` | start the window later in the coda |
| B | `!msnoise config set dtt_width=15` | a 5–20 s window instead of 5–35 s |
| C | `!msnoise config set dtt_mincoh=0.8` | demand better-behaved windows |
| D | `!msnoise config set dtt_maxerr=0.05` | demand smaller errors |
| E | `!msnoise config set dtt_sides=left` | one branch only, not both |
| F | `!msnoise config set mov_stack=5` *(re-run `stack -m` too, and plot with `mov_stack=5`)* | 5-day stacks |

Record three numbers: the **largest dv/v excursion**, **when** it happens,
and **how many days survived**.

Predict now, before you run it: will the seven answers agree?

### What this is really about

Everyone here runs correct code on correct data and gets a different number,
because each of you makes a different *defensible* choice. The result is not
"dv/v = x %". It is "dv/v = x %, in this band, in this lag window, relative
to this reference, with this coherence threshold". A paper that does not
state all five cannot be reproduced.

For scale: De Plaen et al. (2019) report dv/v between +0.3 % and −0.2 % at
Etna with median errors near 0.03 %. If your excursion is several per cent,
look at your coherence before you believe it."""))

c.append(md("""## 8. Before you move on

1. Your dv/v is relative to a reference built from the **whole period**,
   including the eruption. What would change if you used only December?
2. The 1–2 Hz band is also Etna's volcanic tremor band, and tremor comes
   from a fixed point at the summit. The method assumes noise arriving from
   all directions. How comfortable are you?
3. Notebook 04 puts this next to the tremor amplitude and the microseism.
   Before you open it: if dv/v changes, how would you decide whether the
   *medium* changed or the *noise source* did?"""))

build(c, "notebooks/03_autocorrelation_monitoring.ipynb")
print("wrote 03")
