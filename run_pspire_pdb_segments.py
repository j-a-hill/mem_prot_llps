"""
run_pspire_pdb_segments.py

Scores individual topology segment PDB files with PSPire.

For each region, batches all segment PDBs into one PSPire run (one CSV per
region), then splits results back by UniProt_ID + seg_idx using the manifest.

PSPire requires:
  - HEADER record in PDB (prepended from full AF2 PDB)
  - Filenames that embed the UniProt ID in AF-<uid>-F1-model_v<j> format,
    or we supply a project name and parse PSPire's CSV output

Reads:  data/alphafold_pdbs/segment_manifest.csv
Output: output/topology_scores_raw/PSPire_pdb_segments_raw.csv
  UniProt_ID, region, seg_idx, start, end, n_residues, score

Run:
  python3 run_pspire_pdb_segments.py
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT       = Path(__file__).parent
MANIFEST   = ROOT / "data/alphafold_pdbs/segment_manifest.csv"
PDB_DIR    = ROOT / "data/alphafold_pdbs"
PSPIRE_PY  = ROOT / "external_tools/PSPire/PSPire.py"
OUT_CSV    = ROOT / "output/topology_scores_raw/PSPire_pdb_segments_raw.csv"
PYTHON     = "/home/jake/miniconda3/envs/PSPire/bin/python3"


def pdb_header_lines(full_pdb: Path) -> list[str]:
    """Return HEADER/TITLE/COMPND/SOURCE lines from the full AF2 PDB."""
    lines = []
    with open(full_pdb) as fh:
        for line in fh:
            if re.match(r"^(HEADER|TITLE|COMPND|SOURCE)", line):
                lines.append(line)
            elif line.startswith("ATOM"):
                break
    return lines


def build_pdb_with_header(full_pdb: Path, seg_pdb: Path, dest: Path) -> None:
    """dest = HEADER lines from full_pdb + ATOM/TER/END from seg_pdb."""
    lines = pdb_header_lines(full_pdb)
    with open(seg_pdb) as fh:
        lines.extend(fh)
    dest.write_text("".join(lines))


def run_pspire_on_dir(pdb_dir: Path, project_name: str,
                      work_dir: Path) -> pd.DataFrame | None:
    out_csv = work_dir / f"{project_name}_scores.csv"
    cmd = [PYTHON, str(PSPIRE_PY),
           "-d", str(pdb_dir), "-o", str(out_csv), "-n", project_name]
    env = os.environ.copy()
    env["PATH"] = "/usr/bin:" + env.get("PATH", "")

    result = subprocess.run(cmd, cwd=str(work_dir),
                            capture_output=True, text=True, env=env)
    if result.returncode != 0:
        print(f"  PSPire failed (exit {result.returncode})", file=sys.stderr)
        if result.stderr.strip():
            print(result.stderr[-2000:], file=sys.stderr)
        return None
    if not out_csv.exists():
        print(f"  PSPire produced no output: {out_csv}", file=sys.stderr)
        return None
    return pd.read_csv(out_csv)


def main():
    if not PSPIRE_PY.exists():
        sys.exit(f"PSPire not found at {PSPIRE_PY}")

    manifest = pd.read_csv(MANIFEST)
    print(f"Segments in manifest: {len(manifest)}")

    # Build a lookup: uid → full PDB path
    pdb_map = {}
    for pdb in sorted(PDB_DIR.glob("AF-*-F1-model_v*.pdb")):
        m = re.search(r"AF-([A-Z0-9]+)-F1-model_v(\d+)\.pdb$", pdb.name)
        if m:
            uid, ver = m.group(1), m.group(2)
            if uid not in pdb_map or int(ver) > int(pdb_map[uid][1]):
                pdb_map[uid] = (pdb, ver)

    all_rows = []

    for region in ["Cytoplasmic", "Transmembrane", "Extracellular/Lumenal"]:
        region_segs = manifest[manifest.region == region].copy()
        if region_segs.empty:
            continue
        print(f"\n{'='*60}")
        print(f"Region: {region}  ({len(region_segs)} segments)")

        with tempfile.TemporaryDirectory(prefix="pspire_seg_") as tmp:
            tmp = Path(tmp)
            input_dir = tmp / "input"
            input_dir.mkdir()
            work_dir  = tmp / "work"
            work_dir.mkdir()

            # Track (filename_stem → manifest row) for matching results back
            stem_to_row = {}

            for _, seg in region_segs.iterrows():
                uid      = seg["UniProt_ID"]
                seg_idx  = int(seg["seg_idx"])
                seg_pdb  = ROOT / seg["pdb_path"]
                full_pdb = pdb_map.get(uid, (None,))[0]

                if not seg_pdb.exists():
                    print(f"  [missing] {seg_pdb.name}")
                    continue
                if full_pdb is None:
                    print(f"  {uid}: no full PDB, skipping")
                    continue

                # PSPire parses UniProt ID from AF-<uid>-F1-model_v<j>.pdb
                # Use canonical segment name so result CSV has identifiable Uniprot_ID
                prefix_map = {"Cytoplasmic": "cyto", "Transmembrane": "tm",
                              "Extracellular/Lumenal": "ec"}
                pfx = prefix_map[region]
                # Embed seg_idx in the "uid" field so PSPire's output is unique per seg
                pseudo_uid = f"{uid}S{seg_idx:02d}"   # e.g. O00254S01
                ver = pdb_map[uid][1]
                dest_name = f"AF-{pseudo_uid}-F1-model_v{ver}.pdb"
                dest = input_dir / dest_name

                build_pdb_with_header(full_pdb, seg_pdb, dest)
                stem_to_row[pseudo_uid] = seg

            n_pdbs = len(list(input_dir.glob("*.pdb")))
            print(f"  Running PSPire on {n_pdbs} segment PDBs...")
            if n_pdbs == 0:
                continue

            proj = f"pspire_seg_{region.replace('/', '_').replace(' ', '_')}"
            scores_df = run_pspire_on_dir(input_dir, proj, work_dir)

            if scores_df is None:
                print(f"  !! PSPire returned no results for {region}")
                continue

            # PSPire output column: Uniprot_ID (which we set to the pseudo_uid)
            for _, r in scores_df.iterrows():
                pseudo_uid = str(r["Uniprot_ID"])
                if pseudo_uid not in stem_to_row:
                    print(f"  [no match] {pseudo_uid}")
                    continue
                seg = stem_to_row[pseudo_uid]
                score = float(r["Score"])
                print(f"  {seg['UniProt_ID']} {region} seg{int(seg['seg_idx']):02d} "
                      f"({seg['start']}-{seg['end']}, n={seg['n_residues']}) → {score:.4f}")
                all_rows.append({
                    "UniProt_ID":  seg["UniProt_ID"],
                    "region":      region,
                    "seg_idx":     int(seg["seg_idx"]),
                    "start":       seg["start"],
                    "end":         seg["end"],
                    "n_residues":  seg["n_residues"],
                    "score":       score,
                })

    if not all_rows:
        sys.exit("No scores collected")

    df = pd.DataFrame(all_rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {len(df)} rows to {OUT_CSV}")
    print(df.groupby("region")["score"].describe().round(3).to_string())


if __name__ == "__main__":
    main()
