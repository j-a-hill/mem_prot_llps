"""
run_picnic_pdb_region.py

Downloads AlphaFold2 PDB files for the 60 topology-annotated membrane proteins,
slices each to Cytoplasmic / Transmembrane / Extracellular+Lumenal regions
(ATOM records only — no SEQRES, so PICNIC extracts region sequence from ATOM
records and runs IUPred on the isolated-region context), then scores each slice
with PICNIC manual mode (programmatic Python API, no subprocess).

Output: output/topology_scores_raw/PICNIC_pdb_region_raw.csv
  UniProt_ID, region, n_residues, score

Run with the PSPHunter conda env and IUPred2A on PYTHONPATH:

  PYTHONPATH=/home/jake/iupred2a_lib \\
  PATH="/home/jake/miniconda3/envs/PSPHunter/bin:$PATH" \\
  /home/jake/miniconda3/envs/PSPHunter/bin/python3 run_picnic_pdb_region.py

Approach label: 'pdb_region' — parallel to 'concatenated' for sequence tools.
The score reflects the region treated as a standalone protein; pLDDT values from
the AlphaFold B-factor column are preserved from the full-protein model.
STRIDE secondary structure and IUPred disorder are computed on the isolated
region (boundary effects expected at TM junctions; minimal for large regions).
TM region results are expected to be uninformative (hydrophobic helices fully
exposed in isolation — negative control only).
"""

import os
import re
import sys
import time
import warnings
from pathlib import Path

import pandas as pd
import requests

ROOT       = Path(__file__).parent
TOPO_CACHE = ROOT / "data/uniprot_topology_cache.csv"
PDB_DIR    = ROOT / "data/alphafold_pdbs"
REG_DIR    = PDB_DIR / "regions"
OUT_RAW    = ROOT / "output/topology_scores_raw"
OUT_CSV    = OUT_RAW / "PICNIC_pdb_region_raw.csv"

for d in [PDB_DIR, REG_DIR, OUT_RAW]:
    d.mkdir(parents=True, exist_ok=True)

# AF2 API — try v6 first (matching PICNIC's own download logic), then fall back
AF_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_{ver}.pdb"

# Topology region labels and the UniProt annotation keywords that map to each
REGION_KEYWORDS = {
    "Cytoplasmic":           ["cytoplasmic", "nuclear"],
    "Transmembrane":         None,   # handled separately from TRANSMEM column
    "Extracellular/Lumenal": ["extracellular", "lumenal", "perinuclear space"],
}

# Known proteins with no AF2 model (HTTP 404 from EBI API)
NO_AF_MODEL = {"O75445", "Q8WXG9"}

MIN_RESIDUES = 10   # skip region slices shorter than this


# ── Topology parsing ────────────────────────────────────────────────────────

def parse_topology(topo_str, tm_str):
    """Return residue-number sets for each region, or empty sets if unannotated."""
    regions = {"Cytoplasmic": set(), "Transmembrane": set(), "Extracellular/Lumenal": set()}

    if isinstance(topo_str, str):
        for m in re.finditer(r'TOPO_DOM (\d+)\.\.(\d+).*?/note="([^"]*)"', topo_str):
            start, end, note = int(m.group(1)), int(m.group(2)), m.group(3).lower()
            residues = set(range(start, end + 1))
            if any(kw in note for kw in REGION_KEYWORDS["Cytoplasmic"]):
                regions["Cytoplasmic"] |= residues
            elif any(kw in note for kw in REGION_KEYWORDS["Extracellular/Lumenal"]):
                regions["Extracellular/Lumenal"] |= residues

    if isinstance(tm_str, str):
        for m in re.finditer(r'TRANSMEM (\d+)\.\.(\d+)', tm_str):
            regions["Transmembrane"] |= set(range(int(m.group(1)), int(m.group(2)) + 1))

    return regions


# ── PDB download ─────────────────────────────────────────────────────────────

def download_pdb(uid: str) -> Path | None:
    """Download the latest AF2 PDB for uid; cache under PDB_DIR.

    Returns the local Path on success, None if no model exists.
    Checks for any already-cached version before fetching.
    """
    existing = sorted(PDB_DIR.glob(f"AF-{uid}-F1-model_v*.pdb"))
    if existing:
        return existing[-1]

    for ver in ["v6", "v7", "v8", "v4", "v3"]:
        url  = AF_URL_TEMPLATE.format(uid=uid, ver=ver)
        dest = PDB_DIR / f"AF-{uid}-F1-model_{ver}.pdb"
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                dest.write_bytes(r.content)
                print(f"  downloaded {uid} ({ver})")
                return dest
            elif r.status_code == 404:
                continue
        except requests.RequestException as e:
            warnings.warn(f"  {uid}: {e}")
    print(f"  {uid}: no AF2 model found")
    return None


# ── PDB slicing ───────────────────────────────────────────────────────────────

def slice_pdb(src: Path, dest: Path, residue_numbers: set) -> int:
    """Write ATOM records for residue_numbers from src to dest.

    Returns the number of residues written.  SEQRES records are deliberately
    omitted so PICNIC falls back to ATOM-based sequence extraction, ensuring
    IUPred2A runs on the isolated-region sequence.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    written = set()
    lines = []
    with open(src) as fh:
        for line in fh:
            rec = line[:6].rstrip()
            if rec in ("ATOM", "HETATM"):
                try:
                    res_num = int(line[22:26])
                except ValueError:
                    continue
                if res_num in residue_numbers:
                    lines.append(line)
                    written.add(res_num)
            elif rec == "TER":
                # Only include TER if we wrote something before it
                if lines and not lines[-1].startswith("TER"):
                    lines.append(line)
    lines.append("END\n")
    dest.write_text("".join(lines))
    return len(written)


# ── PICNIC scoring ────────────────────────────────────────────────────────────

def score_with_picnic(pdb_file: Path, uid: str, tmp_dir: Path) -> float | None:
    """Call PICNIC manual pipeline on pdb_file; return score or None on error."""
    tmp_dir.mkdir(parents=True, exist_ok=True)
    try:
        from picnic_bio.prediction.inference_model import inference_model_manual_without_go_one
        score, _ = inference_model_manual_without_go_one(
            str(pdb_file), uid, output_path=str(tmp_dir)
        )
        return float(score)
    except Exception as e:
        warnings.warn(f"  PICNIC failed for {pdb_file.name}: {e}")
        return None


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    topo_df = pd.read_csv(TOPO_CACHE)
    uids = topo_df["Entry"].tolist()
    print(f"Processing {len(uids)} proteins")

    records = []

    for uid in uids:
        if uid in NO_AF_MODEL:
            print(f"[skip] {uid}: no AF2 model")
            continue

        # ── download full PDB ─────────────────────────────────────────────
        full_pdb = download_pdb(uid)
        if full_pdb is None:
            print(f"[skip] {uid}: download failed")
            continue
        time.sleep(0.3)   # polite rate-limiting

        # ── parse topology ────────────────────────────────────────────────
        row = topo_df[topo_df["Entry"] == uid].iloc[0]
        regions = parse_topology(row.get("Topological domain"), row.get("Transmembrane"))

        for region, res_set in regions.items():
            if len(res_set) < MIN_RESIDUES:
                print(f"  {uid} [{region}]: only {len(res_set)} residues — skip")
                continue

            # ── slice PDB ─────────────────────────────────────────────────
            # Keep same filename as the full PDB so PICNIC's AF-<uid>-F<i>-v<j> regex works
            reg_label  = region.replace("/", "_")
            region_dir = REG_DIR / reg_label
            region_pdb = region_dir / full_pdb.name
            tmp_dir    = region_dir / f"{uid}_picnic_tmp"

            n_written = slice_pdb(full_pdb, region_pdb, res_set)
            if n_written < MIN_RESIDUES:
                print(f"  {uid} [{region}]: only {n_written} residues in PDB — skip")
                continue

            # ── score ────────────────────────────────────────────────────
            print(f"  {uid} [{region}]: {n_written} residues — scoring...", end=" ", flush=True)
            score = score_with_picnic(region_pdb, uid, tmp_dir)
            if score is not None:
                print(f"{score:.4f}")
                records.append({
                    "UniProt_ID": uid,
                    "region":     region,
                    "n_residues": n_written,
                    "score":      score,
                })
            else:
                print("FAILED")

    df = pd.DataFrame(records)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved {len(df)} rows → {OUT_CSV}")
    print(df.groupby("region")["score"].describe().to_string())


if __name__ == "__main__":
    main()
