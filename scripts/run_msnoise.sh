#!/usr/bin/env bash
# Run MSNoise 1.6.5 for single-station autocorrelation (ZZ) and
# cross-components (EN, EZ, NZ) over the whole teaching window.
#
# Run once, offline. Students get the exported .npz, not this.
#
#   bash scripts/run_msnoise.sh [PROJECT_DIR]
#
# Everything is idempotent: re-running skips what is already done. To start
# over, delete the project directory.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJ="${1:-$REPO/msnoise_project}"
SDS="$REPO/data/sds"
RESP="$REPO/data/processed"
MSNOISE="${MSNOISE:-msnoise}"
PYTHON="${PYTHON:-python3}"
THREADS="${THREADS:-4}"

say() { printf '\n=== %s ===\n' "$*"; }
# MSNoise 1.6.5 prints a pkg_resources deprecation warning on every single
# command. Filter it so real messages are visible.
quiet() { grep -vEi 'pkg_resources|SyntaxWarning|invalid escape|^ *default\[' || true; }

say "0. preflight"
"$PYTHON" - <<'EOF'
import sqlalchemy, sys
major = int(sqlalchemy.__version__.split('.')[0])
if major >= 2:
    sys.exit(f"sqlalchemy {sqlalchemy.__version__} will break `msnoise new_jobs` "
             "('Engine' object has no attribute 'execute'). Need < 2.")
try:
    import pkg_resources  # noqa: F401
except ModuleNotFoundError:
    sys.exit("pkg_resources is missing: setuptools >= 81 removed it and every "
             "msnoise command imports it. Need setuptools < 81.")
print(f"sqlalchemy {sqlalchemy.__version__}, pkg_resources present -- OK")
EOF

say "1. build the SDS archive"
"$PYTHON" "$REPO/scripts/build_sds.py" --channels HHZ,HHN,HHE

say "2. create the project"
mkdir -p "$PROJ"
cd "$PROJ"
if [ ! -f db.ini ]; then
    $MSNOISE db init --tech 1 2>&1 | quiet
else
    echo "db.ini exists, keeping the project"
fi

say "3. configure for single-station AC + SC"
# windsorizing=-1 is one-bit normalisation, and it is NOT the obvious default.
# Measured on this dataset, against the usual 3-sigma clip:
#   usable days (cc>0.7)   21/269  ->  128/269
#   median cc           0.24-0.50  ->  0.65-0.78
#   saturated dv/v             23  ->  0
# Etna hosts persistent volcanic tremor at 0.5-5 Hz radiating from a fixed
# point at the summit -- the opposite of the diffuse field the method
# assumes -- and aggressive temporal normalisation is the standard defence.
WINDSORIZING="${WINDSORIZING:--1}"
"$PYTHON" "$REPO/scripts/msnoise_setup.py" --sds "$SDS" --response "$RESP" \
    --windsorizing="$WINDSORIZING" 2>&1 | quiet

# ORDER MATTERS AND NOTHING WARNS YOU.
# `scan_archive` walks the archive for stations that are already in the
# Station table. Run it before `populate` and it reports "Scanning 0
# directories", inserts nothing, and exits 0 looking perfectly healthy.
say "4. populate the station table (BEFORE scanning)"
$MSNOISE populate 2>&1 | quiet

say "5. scan the archive"
$MSNOISE scan_archive --init 2>&1 | quiet | tail -3

say "6. create the CC jobs"
$MSNOISE new_jobs --init 2>&1 | quiet | tail -3

say "7. compute the correlations"
$MSNOISE -t "$THREADS" compute_cc 2>&1 | quiet | grep -Ei 'job|finish|error' | tail -20

say "8. export to the format the notebooks read"
"$PYTHON" "$REPO/scripts/export_msnoise.py" --project "$PROJ"

say "done"
$MSNOISE info 2>&1 | quiet | sed -n '/^Jobs:/,$p'
