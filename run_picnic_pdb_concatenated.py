"""
run_picnic_pdb_concatenated.py

Scores the structural 'pdb_concatenated' PDBs (built by build_concatenated_pdbs.py)
with PICNIC manual mode. These files hold the SAME true 3D coordinates as the
'pdb_region' PDBs for a given protein/topology-class, but with residues
renumbered into one continuous chain (1..N) and no gaps at fragment
boundaries. This isolates the "does slicing-and-splicing itself change the
score" question (comparing pdb_region vs pdb_concatenated, same fold) from
the "does topology class matter" question (comparing across regions).

Output: output/topology_scores_raw/PICNIC_pdb_concatenated_raw.csv
  UniProt_ID, region, n_residues, n_fragments, score

Run with the PSPHunter conda env's IUPred2A on PYTHONPATH and the custom
STRIDE binary on PATH:

  PYTHONPATH=external_tools/iupred2a/iupred2a \\
  PATH="external_tools/stride_bin:$PATH" \\
  python3 run_picnic_pdb_concatenated.py
"""

import warnings
from pathlib import Path

import pandas as pd

ROOT      = Path(__file__).parent
MANIFEST  = ROOT / "data/alphafold_pdbs/concatenated_manifest.csv"
OUT_RAW   = ROOT / "output/topology_scores_raw"
OUT_CSV   = OUT_RAW / "PICNIC_pdb_concatenated_raw.csv"

OUT_RAW.mkdir(parents=True, exist_ok=True)


def score_with_picnic(pdb_file: Path, uid: str, tmp_dir: Path) -> float | None:
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


def main():
    manifest = pd.read_csv(MANIFEST)
    print(f"Scoring {len(manifest)} concatenated region PDBs with PICNIC")

    records = []
    for i, row in manifest.iterrows():
        uid, region = row["UniProt_ID"], row["region"]
        pdb_path = ROOT / row["pdb_path"]
        tmp_dir = pdb_path.parent / f"{uid}_picnic_tmp"

        print(f"  [{i+1}/{len(manifest)}] {uid} [{region}]: "
              f"{row['n_residues']} res, {row['n_fragments']} frag — scoring...",
              end=" ", flush=True)
        score = score_with_picnic(pdb_path, uid, tmp_dir)
        if score is not None:
            print(f"{score:.4f}")
            records.append({
                "UniProt_ID":  uid,
                "region":      region,
                "n_residues":  row["n_residues"],
                "n_fragments": row["n_fragments"],
                "score":       score,
            })
        else:
            print("FAILED")

    df = pd.DataFrame(records)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved {len(df)} rows -> {OUT_CSV}")
    print(df.groupby("region")["score"].describe().to_string())


if __name__ == "__main__":
    main()
