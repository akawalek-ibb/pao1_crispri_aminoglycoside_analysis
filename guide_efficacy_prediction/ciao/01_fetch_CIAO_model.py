#!/usr/bin/env python3
"""Download the two released CIAO prediction files from a fixed upstream commit.

Run from the project root: python guide_efficacy_prediction/ciao/01_fetch_CIAO_model.py
Uses only the Python standard library. Files are saved beside the predictor in
saved_model/. Existing matching files are reused; differing files are not replaced.
"""
import argparse
import hashlib
from pathlib import Path
import sys
import urllib.request

COMMIT = "2dad9153ea61db49001e2cb48095fbb2eac46003"
BASE = "https://raw.githubusercontent.com/BarquistLab/CRISPRi_guide_efficiency_bacteria/" + COMMIT + "/5_gRNA_design/saved_model/"
FILES = {
    "CRISPRi_headers.sav": "0e4cf56a486705ad7a140ab20f6a476edbdfb067358e5bbf6f0a7b157415abbe",
    "CRISPRi_model.sav": "04bf43713d1ed8f89b016d1de16e9be11b259e2af8b7beebf6e6e694a8ad453d",
}

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-only", action="store_true", help="Check local files without downloading.")
    args = parser.parse_args()
    dest = Path(__file__).resolve().parent / "saved_model"
    if not args.check_only:
        dest.mkdir(parents=True, exist_ok=True)
    failures = []
    for name, expected in FILES.items():
        path = dest / name
        if path.exists():
            if sha256(path) != expected:
                failures.append(name + ": checksum differs; keep your copy and choose the intended model explicitly.")
            else:
                print("Verified:", path)
            continue
        if args.check_only:
            failures.append(name + ": missing; run this script without --check-only.")
            continue
        temp = path.with_suffix(path.suffix + ".part")
        try:
            print("Downloading:", name, flush=True)
            with urllib.request.urlopen(BASE + name, timeout=120) as src, temp.open("wb") as out:
                for chunk in iter(lambda: src.read(1024 * 1024), b""):
                    out.write(chunk)
            if sha256(temp) != expected:
                raise ValueError("download checksum mismatch")
            temp.replace(path)
            print("Verified:", path)
        except Exception as exc:
            if temp.exists():
                temp.unlink()
            failures.append(name + ": " + str(exc))
    for message in failures:
        print("ERROR:", message, file=sys.stderr)
    return 1 if failures else 0

if __name__ == "__main__":
    raise SystemExit(main())
