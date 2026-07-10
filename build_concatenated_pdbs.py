"""
build_concatenated_pdbs.py

Builds the structural "pdb_concatenated" analogue of the sequence-level
'concatenated' approach, for the 3 topology-class regions per protein
(Cytoplasmic / Transmembrane / Extracellular/Lumenal).

Input:  data/alphafold_pdbs/regions/{region}/AF-{uid}-F1-model_v{ver}.pdb
        (true, disjoint AlphaFold coordinates for all annotated spans of a
        given topology class in one protein — built by run_picnic_pdb_region.py;
        residues keep their true/original UniProt numbering, so a multi-span
        region like "4 cytoplasmic loops of a GPCR" has numbering gaps and no
        implied chain continuity between spans.)

Output: data/alphafold_pdbs/concatenated/{region}/AF-{uid}-F1-model_v{ver}.pdb
        Same atoms, same 3D coordinates, same atom order (ascending true
        residue number, matching how the region PDB was sliced) — but
        residues are RENUMBERED consecutively (1, 2, 3, ...) and atom serials
        are reset, with no TER between spans. This creates one nominally
        unbroken polypeptide chain per the method the user selected
        (ask_user, 'Structural splicing method', option 1): the seam between
        spatially distant real loops is now invisible to residue numbering,
        exactly mirroring what sequence-level 'concatenated' does when it
        joins disjoint region substrings into one string before scoring.
        3D coordinates are NOT moved — this is a labeling-only manipulation,
        not physical splicing (rigid-body placement), so PICNIC/PSPire (whose
        STRIDE/pLDDT-context and structural features depend on numbering
        continuity and standard atom serials, not on absolute 3D adjacency
        for every feature) will see a single chain while the coordinates
        still reflect each fragment's true position in the whole protein.

Manifest: data/alphafold_pdbs/concatenated_manifest.csv
  UniProt_ID, region, n_residues, n_fragments, pdb_path

Run:
  python3 build_concatenated_pdbs.py
"""

import re
from pathlib import Path

import pandas as pd

ROOT     = Path(__file__).parent
PDB_DIR  = ROOT / "data/alphafold_pdbs"
REG_DIR  = PDB_DIR / "regions"
CAT_DIR  = PDB_DIR / "concatenated"
MANIFEST = PDB_DIR / "concatenated_manifest.csv"

REGIONS = ["Cytoplasmic", "Transmembrane", "Extracellular_Lumenal"]

# Directory name -> canonical region label matching PICNIC_pdb_region_raw.csv /
# PSPire_pdb_region_raw.csv / *_all_approaches.csv conventions.
REGION_LABEL = {
    "Cytoplasmic":            "Cytoplasmic",
    "Transmembrane":          "Transmembrane",
    "Extracellular_Lumenal":  "Extracellular/Lumenal",
}


def renumber_pdb(src: Path, dest: Path):
    """Read ATOM/HETATM records from src (ascending true residue number,
    possibly with numbering gaps between disjoint spans), and write dest
    with residues renumbered consecutively from 1 and atom serials reset,
    coordinates untouched. Returns (n_residues, n_fragments)."""
    dest.parent.mkdir(parents=True, exist_ok=True)

    lines_out = []
    serial = 0
    new_resnum = 0
    last_true_resnum = None
    n_fragments = 0

    with open(src) as fh:
        for line in fh:
            rec = line[:6].rstrip()
            if rec not in ("ATOM", "HETATM"):
                continue
            try:
                true_resnum = int(line[22:26])
            except ValueError:
                continue

            if true_resnum != last_true_resnum:
                new_resnum += 1
                if last_true_resnum is not None and true_resnum != last_true_resnum + 1:
                    n_fragments += 1
                last_true_resnum = true_resnum

            serial += 1
            new_line = (
                line[:6]
                + f"{serial:>5d}"
                + line[11:22]
                + f"{new_resnum:>4d}"
                + line[26:]
            )
            lines_out.append(new_line)

    if not lines_out:
        return 0, 0

    lines_out.append("TER\n")
    lines_out.append("END\n")
    dest.write_text("".join(lines_out))
    return new_resnum, n_fragments + 1  # +1 for the first fragment


def main():
    manifest_rows = []
    total_written = 0

    for region in REGIONS:
        src_dir = REG_DIR / region
        if not src_dir.is_dir():
            print(f"  {region}: source dir not found, skipping")
            continue

        for src in sorted(src_dir.glob("AF-*-F1-model_v*.pdb")):
            m = re.search(r"AF-([A-Z0-9]+)-F1-model_v(\d+)\.pdb", src.name)
            if not m:
                continue
            uid = m.group(1)

            dest = CAT_DIR / region / src.name
            n_res, n_frag = renumber_pdb(src, dest)
            if n_res == 0:
                print(f"  {uid} {region}: no residues written, skipping")
                continue

            manifest_rows.append({
                "UniProt_ID":  uid,
                "region":      REGION_LABEL[region],
                "n_residues":  n_res,
                "n_fragments": n_frag,
                "pdb_path":    str(dest.relative_to(ROOT)),
            })
            total_written += 1

    manifest = pd.DataFrame(manifest_rows)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(MANIFEST, index=False)

    print(f"\nWrote {total_written} concatenated region PDB files")
    print(manifest.groupby("region")[["n_residues", "n_fragments"]].describe().round(1).to_string())
    print(f"\nManifest: {MANIFEST}")


if __name__ == "__main__":
    main()
