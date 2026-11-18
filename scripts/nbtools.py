"""Tiny helper to build .ipynb files from (kind, source) pairs."""
import json
from pathlib import Path


def md(s):
    return ("markdown", s)


def code(s):
    return ("code", s)


def build(cells, path, kernel="python3", title=None):
    out = []
    for kind, src in cells:
        lines = src.strip("\n").split("\n")
        lines = [l + "\n" for l in lines[:-1]] + [lines[-1]]
        cid = f"c{len(out):03d}"
        if kind == "markdown":
            out.append({"cell_type": "markdown", "id": cid, "metadata": {},
                        "source": lines})
        else:
            out.append({"cell_type": "code", "id": cid, "metadata": {},
                        "execution_count": None, "outputs": [],
                        "source": lines})
    nb = {"cells": out,
          "metadata": {"kernelspec": {"display_name": "Python 3",
                                      "language": "python", "name": kernel},
                       "language_info": {"name": "python",
                                         "version": "3.12"}},
          "nbformat": 4, "nbformat_minor": 5}
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(nb, indent=1))
    return path


COLAB_BADGE = (
    "[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]"
    "(https://colab.research.google.com/github/ORG/etna-noise-practical/blob/"
    "main/notebooks/{name})"
)

BOOTSTRAP = '''# --- run me first ---------------------------------------------------
# Works both in the conda environment and in Google Colab.
import sys, os, importlib.util
if importlib.util.find_spec("obspy") is None:          # Colab
    !pip install -q obspy==1.5.1
if not os.path.exists("etnalib"):                      # Colab
    !git clone -q https://github.com/ORG/etna-noise-practical.git
    %cd etna-noise-practical
    !python scripts/get_data.py
sys.path.insert(0, os.path.abspath("."))
sys.path.insert(0, os.path.abspath(".."))

import numpy as np, matplotlib.pyplot as plt
from etnalib import config as C, io, plots
plots.style()
print("ready -- stations:", ", ".join(C.STATIONS))
# If a later cell says "No module named ...", it is almost always this line
# pointing somewhere other than the etna-noise environment: the notebook is
# running on a different kernel. Kernel > Change kernel, or see SETUP.md.
print("python:", sys.executable)'''
