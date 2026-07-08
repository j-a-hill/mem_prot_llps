"""
integrate_pdb_segments.py

Appends 'segment' approach rows from PICNIC and PSPire pdb_segment runs into
their respective *_all_approaches.csv files, then rebuilds topology_scores_master.csv.

Run after run_picnic_pdb_segments.py and run_pspire_pdb_segments.py:
  python3 integrate_pdb_segments.py
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT    = Path(__file__).parent
RAW_DIR = ROOT / "output/topology_scores_raw"
MASTER  = ROOT / "build_topology_scores_master.py"

TOOLS = {
    "PICNIC": RAW_DIR / "PICNIC_pdb_segments_raw.csv",
    "PSPire": RAW_DIR / "PSPire_pdb_segments_raw.csv",
}


def integrate_tool(tool: str, src: Path) -> None:
    dest = RAW_DIR / f"{tool}_all_approaches.csv"
    if not src.exists():
        print(f"[{tool}] segment file not found: {src} — skipping")
        return
    if not dest.exists():
        print(f"[{tool}] all_approaches file not found: {dest} — skipping")
        return

    raw = pd.read_csv(src)
    print(f"\n[{tool}] segment rows loaded: {len(raw)}")
    print(raw.groupby("region")["score"].describe().round(3).to_string())

    # Build rows matching *_all_approaches schema
    new_rows = raw[["UniProt_ID", "region", "seg_idx", "score"]].copy()
    new_rows["approach"] = "segment"
    new_rows = new_rows[["UniProt_ID", "approach", "region", "seg_idx", "score"]]

    existing = pd.read_csv(dest)
    existing = existing[existing["approach"] != "segment"]

    combined = pd.concat([existing, new_rows], ignore_index=True)
    combined.to_csv(dest, index=False)
    print(f"Updated {dest.name}: {len(combined)} rows")
    print(combined["approach"].value_counts().to_string())


for tool, src in TOOLS.items():
    integrate_tool(tool, src)

print("\nRebuilding topology_scores_master.csv...")
result = subprocess.run([sys.executable, str(MASTER)], capture_output=True, text=True)
print(result.stdout)
if result.returncode != 0:
    print("STDERR:", result.stderr, file=sys.stderr)
    sys.exit(result.returncode)

print("Done.")
