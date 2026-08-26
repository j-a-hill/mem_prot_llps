"""
STEP 1 -- get the source database files into minimal/raw/.

For each of the six sources in config.SOURCES it tries the URL first, and if that
does not work it copies the local snapshot instead. Either way you end up with the
same set of filenames in minimal/raw/, so step 2 does not care which happened.

It writes minimal/raw/PROVENANCE.csv recording, per source, whether the file came
from the web or from disk, and how big it is. That file is your record of what the
run was actually built from.

Run:  python 01_download.py
"""

import json
import shutil
import sys
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))   # find config.py
import config as C

TIMEOUT = 120


def count_records(path, kind):
    """
    Count the data records in a downloaded file. Raises if the file is not the
    format we expected.

    This is the important part of this script. A download can "succeed" -- HTTP 200,
    plausible size -- and still be useless:
      * PhaSepDB's download URL serves an HTML page, not a table.
      * CD-CODE's API is paginated 25 records at a time, so one GET gives you 25 of
        11,144 and looks perfectly valid.
    Checking the size cannot catch either. Parsing it can.
    """
    if kind == "json_records":
        obj = json.loads(path.read_text())
        # CD-CODE wraps its records in {count, current_page, data, pages}; the stored
        # snapshot is a bare list. Accept both.
        if isinstance(obj, dict) and "data" in obj:
            return len(obj["data"])
        return len(obj)
    if kind == "json_by_key":
        return len(json.loads(path.read_text()))
    if kind == "table":
        # sep=None lets pandas sniff comma vs tab -- UniProt serves TSV, the stored
        # snapshot is CSV, and DrLLPS is tab-separated under a .csv name.
        return len(pd.read_csv(path, sep=None, engine="python", on_bad_lines="skip"))
    if kind == "excel":
        return len(pd.read_excel(path))
    raise ValueError(f"unknown kind {kind}")


def fetch(url, dest, kind, min_records):
    """
    Download url to dest, and keep it only if it parses to at least min_records
    records. Returns a short status string for the provenance log.
    """
    try:
        r = requests.get(url, timeout=TIMEOUT)
        r.raise_for_status()
        dest.write_bytes(r.content)
    except Exception as e:
        return f"download failed ({type(e).__name__})"

    try:
        n = count_records(dest, kind)
    except Exception as e:
        return f"downloaded but unparseable ({type(e).__name__})"

    if n < min_records:
        return f"downloaded but only {n} records (expected >={min_records})"
    return "downloaded"


rows = []

for name, src in C.SOURCES.items():
    dest = C.RAW / src["out"]

    status = fetch(src["url"], dest, src["kind"], src["min_records"])

    if status == "downloaded":
        origin, why = "downloaded", ""
    else:
        # Anything other than a clean, complete download falls back to the snapshot.
        shutil.copy(src["local"], dest)
        origin, why = "local snapshot", status

    rows.append({
        "source": name,
        "file": src["out"],
        "origin": origin,
        "n_records": count_records(dest, src["kind"]),
        "size_kb": round(dest.stat().st_size / 1024, 1),
        "fallback_reason": why,
        "url": src["url"],
        "note": src["note"],
    })

prov = pd.DataFrame(rows)
prov.to_csv(C.RAW / "PROVENANCE.csv", index=False)

print(prov[["source", "origin", "n_records", "size_kb"]].to_string(index=False))

fell_back = prov[prov.fallback_reason != ""]
if len(fell_back):
    print("\nfell back to the local snapshot:")
    for _, r in fell_back.iterrows():
        print(f"  {r.source:<11s} {r.fallback_reason}")

print(f"\nwrote {C.RAW / 'PROVENANCE.csv'}")
