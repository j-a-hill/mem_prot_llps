"""
integrate_picnic_pdb_region.py

Merges PICNIC pdb_region scores into:
  1. output/topology_scores_raw/PICNIC_all_approaches.csv   (appends pdb_region rows)
  2. output/topology_scores_master.csv                       (rebuilds from build script)

Run after run_picnic_pdb_region.py completes:
  python3 integrate_picnic_pdb_region.py
"""
from pathlib import Path
import pandas as pd
import subprocess, sys

ROOT    = Path(__file__).parent
RAW_DIR = ROOT / "output/topology_scores_raw"
SRC     = RAW_DIR / "PICNIC_pdb_region_raw.csv"
DEST    = RAW_DIR / "PICNIC_all_approaches.csv"
MASTER  = ROOT / "build_topology_scores_master.py"

if not SRC.exists():
    sys.exit(f"Not found: {SRC} — run run_picnic_pdb_region.py first")

# ── Load pdb_region scores ─────────────────────────────────────────────────
new = pd.read_csv(SRC)
print(f"pdb_region rows loaded: {len(new)}")
print(new.groupby("region")["score"].describe().round(3).to_string())

# ── Build new rows in the all_approaches schema ────────────────────────────
new_rows = new.rename(columns={"score": "score"})[["UniProt_ID", "region", "score"]].copy()
new_rows["approach"] = "pdb_region"
new_rows["seg_idx"]  = None
new_rows = new_rows[["UniProt_ID", "approach", "region", "seg_idx", "score"]]

# ── Append to PICNIC_all_approaches.csv ────────────────────────────────────
existing = pd.read_csv(DEST)
print(f"\nExisting rows: {len(existing)}")

# Remove any stale pdb_region rows (allow re-run)
existing = existing[existing["approach"] != "pdb_region"]

combined = pd.concat([existing, new_rows], ignore_index=True)
combined.to_csv(DEST, index=False)
print(f"Updated PICNIC_all_approaches.csv: {len(combined)} rows")
print(combined["approach"].value_counts().to_string())

# ── Rebuild topology_scores_master.csv ─────────────────────────────────────
print("\nRebuilding topology_scores_master.csv...")
result = subprocess.run([sys.executable, str(MASTER)], capture_output=True, text=True)
print(result.stdout)
if result.returncode != 0:
    print("STDERR:", result.stderr, file=sys.stderr)
    sys.exit(result.returncode)

print("Done.")
