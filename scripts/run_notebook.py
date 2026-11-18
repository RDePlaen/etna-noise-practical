#!/usr/bin/env python3
"""Execute a notebook end to end and report cell-by-cell timing and errors.

Usage:  python scripts/run_notebook.py notebooks/03_....ipynb [--save]

The bootstrap cell contains Colab-only IPython magics (`!pip`, `%cd`) guarded
by an `if`, which the plain kernel parses fine but we strip anyway so the
test matches what a conda user actually runs.
"""
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

ROOT = Path(__file__).resolve().parent.parent


def main():
    path = Path(sys.argv[1])
    save = "--save" in sys.argv
    nb = nbformat.read(path, as_version=4)

    # drop the Colab-install lines; in a real conda session they are skipped
    for cell in nb.cells:
        if cell.cell_type == "code" and "find_spec(\"obspy\")" in cell.source:
            # %matplotlib inline, NOT Agg: under Agg a notebook's plt.show()
            # draws nothing and nbclient captures no images, so a run can
            # pass with zero figures and look fine.
            cell.source = (
                'import sys, os\n'
                'sys.path.insert(0, os.path.abspath("."))\n'
                'sys.path.insert(0, os.path.abspath(".."))\n'
                'get_ipython().run_line_magic("matplotlib", "inline")\n'
                'import numpy as np\n'
                'import matplotlib.pyplot as plt\n'
                'from etnalib import config as C, io, plots\n'
                'plots.style()\n'
                'print("ready -- stations:", ", ".join(C.STATIONS))'
            )

    client = NotebookClient(nb, timeout=900, kernel_name="python3",
                            resources={"metadata": {"path": str(ROOT)}},
                            allow_errors=True)
    t0 = time.time()
    client.execute()
    total = time.time() - t0

    nerr = 0
    for i, cell in enumerate(nb.cells):
        if cell.cell_type != "code":
            continue
        errs = [o for o in cell.get("outputs", [])
                if o.get("output_type") == "error"]
        if errs:
            nerr += 1
            print(f"\n*** CELL {i} FAILED ***")
            print(cell.source[:400])
            print("---")
            print("\n".join(errs[0].get("traceback", []))[-1800:])
    print(f"\n{path.name}: {total:.1f}s total, "
          f"{nerr} failing cell(s) out of "
          f"{sum(1 for c in nb.cells if c.cell_type=='code')}")
    if save:
        nbformat.write(nb, path.with_suffix(".executed.ipynb"))
    return 1 if nerr else 0


if __name__ == "__main__":
    sys.exit(main())
