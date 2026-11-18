#!/usr/bin/env python3
"""Fetch and unpack the prepared data archive.

Deliberately dependency-free (standard library only) so that it runs even if
the conda environment is not quite right yet.
"""
import hashlib
import sys
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# TODO: replace with the Zenodo record once minted.
URL = "https://zenodo.org/records/XXXXXXX/files/etna_practical_data.tar.gz"
SHA256 = None          # fill in when the archive is uploaded
DEST = ROOT / "data"
ARCHIVE = ROOT / "etna_practical_data.tar.gz"


def main():
    if (DEST / "processed" / "inventory.xml").exists() and "--force" not in sys.argv:
        print("data/ already looks populated; pass --force to download again")
        return
    if not ARCHIVE.exists():
        print(f"downloading {URL}")
        with urllib.request.urlopen(URL, timeout=300) as r, \
                open(ARCHIVE, "wb") as f:
            total = int(r.headers.get("Content-Length", 0))
            got = 0
            while chunk := r.read(1 << 20):
                f.write(chunk)
                got += len(chunk)
                if total:
                    print(f"\r  {got/1e6:7.1f} / {total/1e6:.1f} MB", end="")
            print()
    if SHA256:
        h = hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
        if h != SHA256:
            sys.exit(f"checksum mismatch:\n  got      {h}\n  expected {SHA256}")
        print("checksum OK")
    print(f"unpacking into {DEST}")
    with tarfile.open(ARCHIVE) as t:
        t.extractall(ROOT, filter="data")
    print("done. Check with:")
    print("  python -c \"from etnalib import io; print(io.available_days('EMSG')[:3])\"")


if __name__ == "__main__":
    main()
