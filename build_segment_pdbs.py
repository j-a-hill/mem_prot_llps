"""
build_segment_pdbs.py

Slices each full AF2 PDB into individual topology segments (one PDB file per
contiguous annotated span) and writes a manifest CSV.

Output directory:
  data/alphafold_pdbs/segments/{Cytoplasmic|Transmembrane|Extracellular_Lumenal}/
    AF-{uid}-F1-model_v{ver}_cyto00.pdb   (cytoplasmic segment 0)
    AF-{uid}-F1-model_v{ver}_tm00.pdb     (transmembrane segment 0)
    ...

Output manifest:
  data/alphafold_pdbs/segment_manifest.csv
    UniProt_ID, region, seg_idx, start, end, n_residues, pdb_path

Run:
  python3 build_segment_pdbs.py
"""

import re
from pathlib import Path

import pandas as pd

ROOT       = Path(__file__).parent
TOPO_CACHE = ROOT / "data/uniprot_topology_cache.csv"
PDB_DIR    = ROOT / "data/alphafold_pdbs"
SEG_DIR    = PDB_DIR / "segments"
MANIFEST   = PDB_DIR / "segment_manifest.csv"

MIN_RESIDUES = 5   # skip segments shorter than this

REGION_KEYWORDS = {
    "Cytoplasmic":           ["cytoplasmic", "nuclear"],
    "Extracellular/Lumenal": ["extracellular", "lumenal", "perinuclear space"],
}

REGION_DIR = {
    "Cytoplasmic":           "Cytoplasmic",
    "Transmembrane":         "Transmembrane",
    "Extracellular/Lumenal": "Extracellular_Lumenal",
}

REGION_PREFIX = {
    "Cytoplasmic":           "cyto",
    "Transmembrane":         "tm",
    "Extracellular/Lumenal": "ec",
}


def parse_segments(topo_str, tm_str):
    """Return list of (region, start, end) for each individual annotated span."""
    segs = []

    if isinstance(topo_str, str):
        for m in re.finditer(r'TOPO_DOM (\d+)\.\.(\d+).*?/note="([^"]*)"', topo_str):
            start, end, note = int(m.group(1)), int(m.group(2)), m.group(3).lower()
            for region, kws in REGION_KEYWORDS.items():
                if any(kw in note for kw in kws):
                    segs.append((region, start, end))

    if isinstance(tm_str, str):
        for m in re.finditer(r'TRANSMEM (\d+)\.\.(\d+)', tm_str):
            start, end = int(m.group(1)), int(m.group(2))
            segs.append(("Transmembrane", start, end))

    return segs


def slice_pdb(src: Path, dest: Path, residue_numbers: set) -> int:
    """Write ATOM records for residue_numbers from src to dest. Returns n_residues."""
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
                if lines and not lines[-1].startswith("TER"):
                    lines.append(line)
    lines.append("END\n")
    dest.write_text("".join(lines))
    return len(written)


def main():
    topo = pd.read_csv(TOPO_CACHE)

    # Find all available full PDBs
    pdb_map = {}
    for pdb in sorted(PDB_DIR.glob("AF-*-F1-model_v*.pdb")):
        m = re.search(r"AF-([A-Z0-9]+)-F1-model_v(\d+)\.pdb", pdb.name)
        if m:
            uid, ver = m.group(1), m.group(2)
            # keep highest version
            if uid not in pdb_map or int(ver) > int(pdb_map[uid][1]):
                pdb_map[uid] = (pdb, ver)

    manifest_rows = []
    total_written = 0

    for _, row in topo.iterrows():
        uid = row["Entry"]
        if uid not in pdb_map:
            print(f"  {uid}: no AF2 PDB found, skipping")
            continue

        full_pdb, ver = pdb_map[uid]
        segments = parse_segments(row["Topological domain"], row["Transmembrane"])

        # Count by region to assign seg_idx
        region_counts = {}

        for region, start, end in segments:
            n_expected = end - start + 1
            if n_expected < MIN_RESIDUES:
                print(f"  {uid} {region} {start}-{end}: {n_expected} res < {MIN_RESIDUES}, skipping")
                continue

            seg_idx = region_counts.get(region, 0)
            region_counts[region] = seg_idx + 1

            prefix = REGION_PREFIX[region]
            out_dir = SEG_DIR / REGION_DIR[region]
            out_name = full_pdb.stem + f"_{prefix}{seg_idx:02d}.pdb"
            dest = out_dir / out_name

            res_set = set(range(start, end + 1))
            n_written = slice_pdb(full_pdb, dest, res_set)

            manifest_rows.append({
                "UniProt_ID":  uid,
                "region":      region,
                "seg_idx":     seg_idx,
                "start":       start,
                "end":         end,
                "n_residues":  n_written,
                "pdb_path":    str(dest.relative_to(ROOT)),
            })
            total_written += 1

    manifest = pd.DataFrame(manifest_rows)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(MANIFEST, index=False)

    print(f"\nWrote {total_written} segment PDB files")
    print(manifest.groupby("region")["n_residues"].describe().round(1).to_string())
    print(f"\nManifest: {MANIFEST}")


if __name__ == "__main__":
    main()
