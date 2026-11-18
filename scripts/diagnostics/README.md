# Diagnostics

Scripts used to design and validate the practical. Not needed to run it, but
they document how the defaults were chosen — and each one found something.

| script | what it establishes |
|---|---|
| `check_hvsr.py` | f0 and SESAME score per station |
| `check_acf.py` | per-day autocorrelation timing and day-to-day repeatability |
| `check_stability_time.py` | **where** in the window the ACFs decorrelate (24–27 Dec) and that adjacent days stay coherent throughout |
| `sweep_acf.py` | adjacent-day correlation across bands, normalisations and lag windows |
| `check_lagwin.py` | how coda window trades correlation against dv/v stability |
| `check_dvv_full.py` | dv/v and cc per station over the whole window |
| `check_response_metadata.py` | stated sensitivity vs product of stage gains, per station |

Run any of them from the repository root, e.g.

    python scripts/diagnostics/check_dvv_full.py
