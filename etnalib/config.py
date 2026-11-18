"""Configuration for the Etna ambient-noise practical.

Everything that a student or a demonstrator might want to change lives here,
so that the notebooks themselves stay readable.

Window
------
2025-12-10 to 2026-02-20 was chosen because three independent things happen
inside it, and the practical is about telling them apart:

  * 24 Dec 2025  Etna enters a continuous eruptive phase
  *  1 Jan 2026  a fissure opens at ~2100 m in the Valle del Bove (effusive)
  * 16-23 Jan 2026  Cyclone Harry: record secondary-microseism amplitudes in
                    the Ionian Sea (Scicchitano, Cannata et al., Comm. Earth
                    Environ., 2026)
  *  6 Feb 2026  a new eruption begins at the NE base of the SE Crater

So a velocity change measured in January has at least two candidate causes,
one volcanic and one meteorological, and the students have to argue about it.
"""
from __future__ import annotations
import os
from pathlib import Path

# ----------------------------------------------------------------- paths
ROOT = Path(os.environ.get("ETNA_ROOT", Path(__file__).resolve().parent.parent))
DATA = ROOT / "data"
RAW = DATA / "raw"
PROC = DATA / "processed"
FIG = ROOT / "figures"
for _d in (RAW, PROC, FIG):
    _d.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------- FDSN
FDSN_BASE = "https://webservices.ingv.it/fdsnws"
FDSN_CLIENT = "INGV"          # for obspy.clients.fdsn.Client
NETWORK = "IV"

# ----------------------------------------------------------------- window
START = "2025-12-10"
END = "2026-02-20"

# Days downloaded in full (24 h). Everything else is 00:00-06:00 UTC only.
FULL_DAY_BLOCKS = [
    ("2025-12-28", "2026-01-03"),   # eruption onset / New Year fissure
    ("2026-01-14", "2026-01-27"),   # Cyclone Harry and the storm train
    ("2026-02-04", "2026-02-09"),   # 6 February eruption
]
PARTIAL_HOURS = (0, 6)              # UTC, the night window

# Hours used to build the daily autocorrelations and the HVSR curves.
# None = the whole day.
#
# This started as (0, 6) because the archive only held night hours on most
# days, and mixing a 24 h day with a 6 h day is not a like-for-like
# comparison: the ambient field has a strong daily cycle, so the two stack
# different mixtures of sources and the day-to-day correlation falls from
# ~0.9 to ~0.4. Now that every day is complete, the whole day is used and
# the diurnal cycle becomes something to look at rather than something to
# avoid. Set it back to (0, 6) to reproduce the night-only results.
ACORR_HOURS = None

# ----------------------------------------------------------------- stations
# Four stations, four jobs. Elevation in m, distance to the Ionian coast in km
# (straight line, for the microseism discussion).
STATIONS = {
    "ECPN": dict(lat=37.7437, lon=14.9865, elev=3038, coast_km=22.0,
                 site="Etna Cratere del Piano",
                 role="summit: volcanic tremor, the monitoring target"),
    "ESLN": dict(lat=37.6934, lon=14.9744, elev=1787, coast_km=23.0,
                 site="Serra La Nave",
                 role="mid flank: the classic Etna reference site"),
    "EMSG": dict(lat=37.8208, lon=14.9498, elev=1435, coast_km=29.0,
                 site="EMSG, NW flank",
                 role="NW flank: the cleanest record in the window"),
    "EFIU": dict(lat=37.7896, lon=15.2103, elev=97, coast_km=3.0,
                 site="Fiumefreddo di Sicilia",
                 role="coastal plain: strongest microseism, softest ground"),
}
CHANNEL = "HHZ"
CHANNELS_3C = ["HHZ", "HHN", "HHE"]

# Two quiet hours used for the HVSR exercise (3 components, all stations).
HVSR_WINDOW = ("2026-02-15T01:00:00", "2026-02-15T03:00:00")

# ----------------------------------------------------------------- bands
BANDS = {
    "primary":   (0.05, 0.10),   # primary microseism,   10-20 s
    "secondary": (0.10, 0.35),   # secondary microseism,  3-10 s  <- Harry
    "tremor":    (0.50, 5.00),   # Etna volcanic tremor
    "anthro":    (5.0, 20.0),    # cultural / anthropogenic noise
}
# Band used for the autocorrelation monitoring. Deliberately a parameter:
# notebook 03 asks different groups to change it and compare.
ACORR_BAND = (1.0, 3.0)
ACORR_MAXLAG = 20.0       # s
ACORR_FS = 20.0           # Hz, after decimation
ACORR_WINDOW = 1800.0     # s, length of the sub-windows stacked into a day

# Stretching search for dv/v
DVV_EPS = 0.05            # +/- 5 %
DVV_NEPS = 201
DVV_LAG = (3.0, 15.0)     # s, coda window used for the measurement

# ----------------------------------------------------------------- chronology
# Dates used to annotate every time-series figure. Sources are listed in
# README.md; these are the published INGV / Global Volcanism Program dates.
# (date, kind, short label for figures, full description)
CHRONOLOGY = [
    ("2025-12-24", "volcanic", "eruptive phase",
     "Etna enters a continuous eruptive phase"),
    ("2026-01-01", "volcanic", "VdB fissure",
     "Fissure opens at ~2100 m in the Valle del Bove"),
    ("2026-01-16", "storm", "Harry",
     "Cyclone Harry reaches the Ionian Sea"),
    ("2026-01-21", "storm", "16.6 m wave",
     "Record 16.6 m individual wave, Strait of Sicily"),
    ("2026-01-23", "storm", "Harry ends",
     "Cyclone Harry decays"),
    ("2026-01-27", "storm", "Kristin",
     "Storm Kristin, part of the Jan-Feb 2026 storm train"),
    ("2026-02-06", "volcanic", "new eruption",
     "New eruption at the NE base of the SE Crater"),
]

COLORS = dict(volcanic="#c1440e", storm="#1f6f8b", neutral="#444444")
