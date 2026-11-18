# Setting up — read this before the practical

You need **two things working before you walk in**: a Python environment, and
the data. Allow 20 minutes. If you have never installed Python before, follow
part A exactly; it is written for that case.

If anything fails, do not spend the evening fighting it — go to **Part D
(Google Colab)**, which needs nothing but a browser and a Google account, and
tell the demonstrator at the start of the session.

---

## Part A — Install Miniforge (Python + conda)

Miniforge is a small installer that gives you Python and the `conda` package
manager, already pointed at the scientific-software repository we need. It does
not interfere with any other Python on your machine.

**Windows**

1. Download the installer: <https://conda-forge.org/download/> → *Windows x86_64*
2. Run it. Accept the defaults. **Install for "Just Me"**, and do not tick
   "Add to PATH" if it is offered.
3. When it finishes, open **Miniforge Prompt** from the Start menu. Every
   command below goes in that window, not in PowerShell or cmd.

**macOS**

1. Download the installer for your chip from <https://conda-forge.org/download/>
   — *macOS arm64* for M1/M2/M3/M4, *macOS x86_64* for an Intel Mac. If unsure:
   Apple menu → About This Mac.
2. Open **Terminal** and run, replacing the filename with what you downloaded:
   ```
   bash ~/Downloads/Miniforge3-MacOSX-arm64.sh
   ```
3. Answer `yes` to the licence and `yes` when it offers to initialise. Then
   **close Terminal and open a new one**.

**Linux**

```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```
Answer `yes` to initialise, then open a new terminal.

**Check it worked.** In a new terminal (Miniforge Prompt on Windows):
```
conda --version
```
You should see something like `conda 24.x.x`. If you see "command not found",
the shell was not initialised — close and reopen the terminal first, and if it
still fails, go to Part D.

---

## Part B — Create the environment

Download or clone this repository, then from inside its folder:

```
conda env create -f environment.yml
conda activate etna-noise
```

The first command downloads about 600 MB and takes 3–12 minutes. It includes
MSNoise, which is used to compute the correlations.

> **If that command fails to solve, or runs for more than ~15 minutes**, use
> the fallback instead:
> ```
> conda env create -f environment-minimal.yml
> conda activate etna-noise-minimal
> ```
> It is the same environment without MSNoise. Nothing in notebooks 01–04
> imports MSNoise — the correlations were computed with it offline and ship
> as a data file — so the practical is identical. Tell the demonstrator which
> one you used. **You must
run `conda activate etna-noise` every time you open a new terminal** — this is
the single most common thing that goes wrong.

**Check it worked:**
```
python -c "import obspy, numpy, scipy, matplotlib; print('obspy', obspy.__version__)"
```
Expected: `obspy 1.5.1` (any 1.5.x is fine).

Then run the self-test, which takes a few seconds and checks the physics, not
just the imports:
```
python tests/test_acorr.py
python -m etnalib.dvv
```
Every line should end in `OK`.

---

## Part C — Get the data

The data archive is **not** in this repository — it is too large for git.

```
# from inside the repository folder
python scripts/get_data.py
```

This downloads one compressed archive (**~175 MB**) and unpacks it into
`data/`. It takes a minute or two on a normal connection.

**Check it worked:**
```
python -c "from etnalib import io; print(len(io.available_days('EFIU')), 'days of EFIU waveforms')"
```
Expected: `14 days of EFIU waveforms`. Notebook 03 also needs
`data/sds_teaching/` (61 day files), which is in the same archive. The archive
carries precomputed products for all four stations across all 72 days, plus
raw waveforms only for the days notebook 01 actually plots — that is why this
number is small and nothing is wrong.

Raw waveforms, by station:

| station | days shipped |
|---|---|
| `EFIU` | 2026-01-14 … 2026-01-26 (the storm fortnight) |
| `ECPN`, `EMSG` | 2026-01-15 and 2026-01-20 |
| `ESLN` | 2026-01-15 (the station was down on the 20th) |

So in notebook 01 you can set `STA` to any of the four, as long as `DAY` is
one the station has. If it is not, `io.require_day` says so and lists the
days that are there, rather than failing later with `'NoneType' object has
no attribute 'stats'`.

> **Do not run `notebooks/00_download_data.ipynb` before the session.** That
> notebook rebuilds the archive from scratch, straight from the INGV FDSN
> server, and it takes about **two hours**. It is there so you can reproduce
> or extend the dataset later, and so you can see exactly where every byte
> came from. It is not a setup step.

---

## Part D — Plan B: Google Colab

Needs only a browser and a Google account. Open
`notebooks/01_what_is_ambient_noise.ipynb` on GitHub and click the
**"Open in Colab"** badge at the top, then run the first cell:

```python
!pip install -q obspy==1.5.1
!git clone -q https://github.com/<ORG>/etna-noise-practical.git
%cd etna-noise-practical
!python scripts/get_data.py
```

Notes:
* If Colab shows a **"Restart session"** banner after the install, click it,
  then re-run the cell. This is normal and harmless.
* Colab forgets everything if you leave it idle for a while. If your variables
  vanish, re-run the first cell and carry on.
* Everything works in Colab except interactive matplotlib widgets, which we
  do not use.

---

## `ModuleNotFoundError` even though the package is installed

This means the notebook is running on a **different Python** from the one
`conda install` put the package into — the environment is fine, the kernel is
wrong. It happens when JupyterLab was launched from `base` instead of
`etna-noise`.

The first cell of every notebook prints `python:` followed by the interpreter
in use. If that path does not contain `etna-noise`, that is the problem.

Two ways to fix it:

```bash
# 1. launch Jupyter from inside the environment (preferred)
conda activate etna-noise
jupyter lab
```

```bash
# 2. or register the environment as a kernel, once, then pick it in
#    Kernel > Change kernel... inside the notebook
conda activate etna-noise
python -m ipykernel install --user --name etna-noise --display-name "etna-noise"
```

## If you get stuck

Bring your laptop anyway. A broken environment is a five-minute fix with
someone looking at it, and Colab is always there as a fallback. Please do not
skip the session because the install failed.
