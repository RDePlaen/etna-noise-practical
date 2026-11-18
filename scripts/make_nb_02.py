import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nbtools import md, code, build, BOOTSTRAP, COLAB_BADGE

c = []
c.append(md(f"""# 02 — H/V: what is under the station?

{COLAB_BADGE.format(name="02_hvsr.ipynb")}

**~15 minutes.** The first real *imaging* method of the practical, and the
cheapest one in seismology: one three-component sensor, one hour of ambient
noise, no source, no array, no earthquake.

The idea (Nakamura 1989): if a soft layer sits on stiff bedrock, it traps and
amplifies horizontal ground motion at a resonant frequency, while the
vertical component is much less affected. So the **ratio** H/V peaks at that
resonance, and the resonance tells you about the layer.

For a layer of thickness *h* and shear velocity *Vs* over a much stiffer
half-space:

$$ f_0 = \\frac{{V_s}}{{4h}} $$

the quarter-wavelength rule. One peak, one number about the ground."""))

c.append(code(BOOTSTRAP))

c.append(md("""## 1. Does the method work? Test it on a case with a known answer.

Before applying anything to real data, run it on a problem whose answer you
already know. `hvsr.layer_transfer` is the *exact* SH transfer function of
one layer over a half-space — the physics, with no noise and no estimation."""))

c.append(code('''from etnalib import hvsr

h_true, vs_true = 30.0, 300.0          # 30 m of soft stuff, Vs = 300 m/s
f = np.logspace(-1, 1.6, 4000)
A = hvsr.layer_transfer(f, h_true, vs_true, vs_half=1200.0)

f0, A0, i0 = hvsr.pick_f0(f, A, fmin=0.2, fmax=30, min_amp=1.5)

fig, ax = plt.subplots(figsize=(8, 3.4))
ax.loglog(f, A, color="#333", lw=1.6)
ax.axvline(vs_true / (4 * h_true), color="#c1440e", ls="--",
           label=f"theory  Vs/4h = {vs_true/(4*h_true):.2f} Hz")
ax.plot(f0, A0, "o", color="#c1440e", ms=7, label=f"picked f0 = {f0:.2f} Hz")
ax.set_xlabel("frequency [Hz]"); ax.set_ylabel("amplification")
ax.legend(); ax.set_title("one layer over a half-space — the exact answer")
plt.show()

print(f"recovered thickness: {hvsr.depth_from_f0(f0, vs_true):.1f} m "
      f"(true {h_true:.1f} m)")'''))

c.append(md("""Note the smaller peaks at **3×, 5×, 7×** the fundamental. Those are the
overtones of the same layer, not other layers — a standard way to be wrong
about a site.

**Question:** the peak amplitude here is about 4.4. What does it depend on?
Try changing `vs_half` from 1200 to 600 m/s and re-running the cell."""))

c.append(md("""## 2. Real data: two quiet hours at four Etna stations

Now the real thing. Each station gives us Z, N and E, which we cut into
60-second windows, taper, and transform. Three choices matter and they are
all visible in the code:

* horizontals combined as $\\sqrt{(E^2+N^2)/2}$
* **Konno–Ohmachi** smoothing with b = 40 — a smoothing that is constant in
  *log* frequency, so it treats 0.3 Hz and 30 Hz alike
* averaging over windows **geometrically** (mean of the logs), because H/V is
  a ratio and ratios are log-normally distributed

We do **not** remove the instrument response: H/V is a ratio of channels from
the same sensor, so it cancels."""))

c.append(code('''results = {}
for sta in C.STATIONS:
    tr = io.load_hvsr_window(sta)
    if tr is None:
        print(f"{sta}: no data"); continue
    z, n, e = tr
    r = hvsr.compute_hvsr(z, n, e, win_s=60.0, overlap=0.5, b=40.0,
                          fmin=0.2, fmax=25.0, nf=160)
    if r is None:
        print(f"{sta}: not enough clean windows"); continue
    results[sta] = r
    print(f"{sta}: {r['nwin']:3d} windows kept, f0 = {r['f0']:.2f} Hz, "
          f"A0 = {r['A0']:.2f}, SESAME {r['sesame']['n_sesame_passed']}/5")'''))

c.append(code('''fig, axes = plt.subplots(1, len(results), figsize=(3.4 * len(results), 3.8),
                         sharey=True, sharex=True)
axes = np.atleast_1d(axes)
for ax, (sta, r) in zip(axes, results.items()):
    ax.fill_between(r["f"], r["hv"] / r["sd"], r["hv"] * r["sd"],
                    color=plots.station_color(sta), alpha=0.25)
    ax.loglog(r["f"], r["hv"], color=plots.station_color(sta), lw=1.8)
    ax.axhline(1.0, color="k", lw=0.6, ls=":")
    if np.isfinite(r["f0"]):
        ax.plot(r["f0"], r["A0"], "v", color="#c1440e", ms=8)
        ax.annotate(f"{r['f0']:.2f} Hz", (r["f0"], r["A0"]),
                    textcoords="offset points", xytext=(4, 6),
                    color="#c1440e", fontsize=9)
    ax.set_title(f"{sta}  ({C.STATIONS[sta]['elev']:.0f} m)", fontsize=10)
    ax.set_xlabel("frequency [Hz]")
axes[0].set_ylabel("H / V"); axes[0].set_ylim(0.2, 20)
plt.suptitle("H/V, geometric mean and geometric standard deviation", y=1.02)
plt.show()'''))

c.append(md("""## 3. Read the four curves — and be surprised

Before you convert anything to a depth, look. And compare what you see with
the prediction you made at the end of notebook 01.

Most people predict that the **coastal** station, lowest and furthest from
the volcano, will have the softest ground and the biggest peak. The data says
the opposite: the three mountain stations have clean, isolated peaks, and the
coastal station has essentially none.

Two things to take from that.

**A flat H/V near 1 is not a failed measurement.** It is the correct answer
for a station with no strong impedance contrast beneath it in the frequency
band you are looking at — either because it sits on competent rock, or
because its contrast is too deep (f₀ below your band) or too shallow (f₀
above it). "No peak" is a result; "no peak, therefore the method failed" is
not.

**On a volcano, a peak is not automatically a sedimentary layer.** The
summit stations sit on hundreds of metres of fractured lava and loose
volcaniclastic material over denser rock, which is a perfectly good impedance
contrast. And a station on a steep summit cone also experiences *topographic*
amplification, which can mimic or distort a site resonance. Keep both in mind
before you quote a layer thickness.

### Look closely at the coastal station

Its H/V is not perfectly flat — there is a weak bump near 0.25 Hz. Compare
that frequency with the bands you identified in notebook 01.

It sits in the middle of the **secondary microseism band**. This is a
well-known trap: when one direction of arriving energy dominates the noise
field, H/V can develop a peak that reflects the *source* rather than the
site. This is precisely why published criteria exist, and why the criterion
about scatter between windows matters.

### The SESAME criteria

A peak is only a *reliable* peak if it satisfies published criteria (SESAME
2004) on its amplitude, its sharpness, and the scatter between windows. The
code reports them per station rather than silently dropping anything, so you
can see **which** criterion fails."""))

c.append(code("""# Plain stdlib table -- no pandas, so this cell cannot fail on a
# workstation whose Jupyter kernel is not the etna-noise environment.
def tick(b):
    return " yes" if b else "  no"

hdr = (f'{"station":<8}{"elev":>6}{"f0 [Hz]":>9}{"A0":>7}{"wins":>6}'
       f'{"resolv":>7}{"A0>2":>7}{"low":>7}{"high":>7}{"sigma":>7}{"passed":>8}')
print(hdr)
print("-" * len(hdr))
for sta, r in results.items():
    s = r["sesame"]
    f0 = f'{s["f0"]:9.3f}' if np.isfinite(s["f0"]) else f'{"--":>9}'
    A0 = f'{s["A0"]:7.2f}' if np.isfinite(s["A0"]) else f'{"--":>7}'
    print(f'{sta:<8}{C.STATIONS[sta]["elev"]:>6}{f0}{A0}'
          f'{s["n_windows"]:>6}'
          f'{tick(s["ok_f0_window"]):>7}{tick(s["ok_amplitude"]):>7}'
          f'{tick(s["ok_lowside"]):>7}{tick(s["ok_highside"]):>7}'
          f'{tick(s["ok_sigma"]):>7}'
          f'{str(s["n_sesame_passed"]) + "/5":>8}')
"""))

c.append(md("""## 4. From a frequency to a depth — and the honest caveat

$h = V_s / (4 f_0)$ needs $V_s$, and **we have not measured $V_s$**. H/V alone
gives one number and the equation has two unknowns. So the depth you get is
only as good as the velocity you assumed.

Rather than pick one value and pretend, look at the whole range."""))

c.append(code('''vs_range = np.array([150, 200, 300, 400, 600, 800])
fig, ax = plt.subplots(figsize=(7.5, 4))
for sta, r in results.items():
    if not np.isfinite(r["f0"]):
        continue
    ax.plot(vs_range, hvsr.depth_from_f0(r["f0"], vs_range), "o-",
            color=plots.station_color(sta), label=f"{sta} (f0={r['f0']:.2f} Hz)")
ax.set_xlabel("assumed Vs of the soft layer  [m/s]")
ax.set_ylabel("implied thickness h = Vs / (4 f0)  [m]")
ax.legend(fontsize=8)
ax.set_title("the same measurement, six different assumptions")
plt.show()

for sta, r in results.items():
    if np.isfinite(r["f0"]):
        lo = hvsr.depth_from_f0(r["f0"], 150)
        hi = hvsr.depth_from_f0(r["f0"], 800)
        print(f"{sta}: f0 = {r['f0']:5.2f} Hz  ->  h between {lo:6.1f} and {hi:6.1f} m")'''))

c.append(md("""**This is the honest result**, and it is worth sitting with: a factor of ~5
in assumed velocity gives a factor of ~5 in depth. H/V is superb at telling
you *whether* there is a contrast and *roughly where*, and it cannot give you
a depth on its own.

To break the ambiguity you need an independent constraint — a borehole, an
active-source survey, or a dispersion curve from array noise measurements.
That last one is where **cross-correlation** comes in, and it is why real
studies combine the two."""))

c.append(md("""## 5. Why this matters in Catania

H/V is not an academic exercise here. Site amplification is why two buildings
in the same earthquake, a few hundred metres apart, can have completely
different fates — the one on sediment shakes harder and for longer than the
one on rock.

The towns on Etna's flanks sit on a patchwork of lava flows, volcaniclastic
deposits and alluvium, and they are exposed both to regional earthquakes and
to shallow, local volcano-tectonic ones. That is exactly why OGS deployed the
`15` network — 18 short-period stations in 18 municipalities around Etna, from
November 2022 to October 2024 — purely for microzonation.

**That network is open and on FDSN.** Extending this notebook to it, and
mapping f₀ across the municipalities, is a real project and a good thesis
topic."""))

c.append(md("""## 6. Bonus: does f0 change over three months?

Everything above used two quiet hours. With complete days for all three
components, H/V can be computed **per day** across the whole window — and
**per hour**, which resolves the daily cycle.

```
python scripts/make_overview_figures.py hvsr
```

![continuous H/V](../figures/overview_hvsr_continuous.png)

Measured, with the peak tracked inside a per-station band and refined to
sub-bin precision:

| station | days with a pick | median f0 | inter-quartile range |
|---|---|---|---|
| `EMSG` | 72/72 | 1.484 Hz | 0.016 Hz (**1.1 %**) |
| `ESLN` | 59/65 | 0.561 Hz | 0.046 Hz (8 %) |
| `ECPN` | 54/62 | 1.128 Hz | 0.258 Hz (23 %) |
| `EFIU` | 0/72 | — | no resonance in band |

**The answer is no.** At `EMSG`, the cleanest station, f0 stays within about
1 % of its own median for three months, and does not step at the eruption
onset or during the storm. The site resonance is a property of the
*structure*, and the structure did not change.

Now put that next to notebook 03, where dv/v moves by **0.1–0.4 %**. Those
velocity changes are three to twenty times smaller than the daily scatter of
f0 — so H/V could not possibly have detected them, which is the answer to
question 3 below. Two methods, same station, same data, sensitivities an
order of magnitude apart.

### And over the day?

![H/V by hour of day](../figures/diurnal_hvsr.png)

Stacking the hourly curves by **hour of day** answers a different question:
does the resonance move with the clock? At `EMSG`, the station whose peak is
cleanest, **f0 varies by 0.8 % over the 24 hours** — flat. The *amplitude* of
H/V at high frequency clearly brightens during the day, because the noise
field changes, but the resonant frequency does not.

That is the expected and reassuring answer: f0 is a property of the
structure, and the structure does not care what time it is. Treat the large
numbers quoted for `ECPN` and `ESLN` with suspicion rather than excitement —
look at where their picks drop out, and see the trap note below.

> **Two traps this exercise walks into, both worth knowing.** Picking f0 day
> after day across the full band makes the pick hop between the fundamental
> and an overtone whenever their amplitudes cross — the first attempt gave
> `ECPN` an apparent 89 % scatter that was entirely the picker. And a peak
> picked as a *grid index* can only take discrete values: on this 140-point
> log grid that is 3.5 % quantisation, which makes a stability plot show the
> grid rather than the ground. `hvsr.refine_peak` fixes the second; the
> per-station tracking band fixes the first."""))

c.append(md("""---

## Before you move on

1. One station probably has no peak at all. Write down, in one sentence, what
   that tells you about the ground beneath it.
2. If you had to pick **one** extra measurement to turn your f₀ into a real
   depth, what would it be?
3. Notebook 03 measures a change of a *few tenths of a percent* in seismic
   velocity. Look again at the scatter bands on your H/V curves. Do you think
   H/V could detect a 0.1 % velocity change? What would you need instead?"""))

build(c, "notebooks/02_hvsr.ipynb")
print("wrote 02")
