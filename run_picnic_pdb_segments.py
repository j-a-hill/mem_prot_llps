"""
run_picnic_pdb_segments.py

Scores individual topology segment PDB files with PICNIC.

Reads:  data/alphafold_pdbs/segment_manifest.csv  (built by build_segment_pdbs.py)
Output: output/topology_scores_raw/PICNIC_pdb_segments_raw.csv
  UniProt_ID, region, seg_idx, start, end, n_residues, score

Run:
  PYTHONPATH=/home/jake/iupred2a_lib \\
  PATH="/home/jake/miniconda3/envs/PSPHunter/bin:$PATH" \\
  /home/jake/miniconda3/envs/PSPHunter/bin/python3 run_picnic_pdb_segments.py
"""

import shutil
import tempfile
import warnings
from pathlib import Path

import pandas as pd

ROOT     = Path(__file__).parent
MANIFEST = ROOT / "data/alphafold_pdbs/segment_manifest.csv"
OUT_CSV  = ROOT / "output/topology_scores_raw/PICNIC_pdb_segments_raw.csv"


def score_with_picnic(pdb_file: Path, uid: str, tmp_dir: Path) -> float | None:
    tmp_dir.mkdir(parents=True, exist_ok=True)
    try:
        from picnic_bio.prediction.inference_model import inference_model_manual_without_go_one
        score, _ = inference_model_manual_without_go_one(
            str(pdb_file), uid, output_path=str(tmp_dir)
        )
        return float(score)
    except KeyboardInterrupt:
        raise
    except BaseException as e:
        # STRIDE calls sys.exit() on failure (raises SystemExit, not Exception)
        warnings.warn(f"  PICNIC/STRIDE failed for {pdb_file.name}: {type(e).__name__}: {e}")
        return None


def main():
    manifest = pd.read_csv(MANIFEST)
    print(f"Segments to score: {len(manifest)}")
    print(manifest.groupby("region").size().to_string())

    # Resume: skip already-scored rows
    done = set()
    if OUT_CSV.exists():
        existing = pd.read_csv(OUT_CSV)
        for _, r in existing.iterrows():
            done.add((r["UniProt_ID"], r["region"], int(r["seg_idx"])))
        print(f"\nResuming — {len(done)} already done")

    rows = []
    for _, seg in manifest.iterrows():
        uid     = seg["UniProt_ID"]
        region  = seg["region"]
        seg_idx = int(seg["seg_idx"])
        pdb_path = ROOT / seg["pdb_path"]

        if (uid, region, seg_idx) in done:
            continue

        if not pdb_path.exists():
            print(f"  [missing] {pdb_path.name}")
            continue

        # PICNIC needs the filename to be AF-<uid>-F<i>-model_v<j>.pdb to parse uid.
        # Copy to a temp dir with the canonical name before scoring.
        with tempfile.TemporaryDirectory(prefix="picnic_seg_") as tmp:
            tmp = Path(tmp)
            # strip the _cyto00 etc. suffix for PICNIC's filename parser
            stem_base = pdb_path.stem.split("_cyto")[0].split("_tm")[0].split("_ec")[0]
            canonical = tmp / f"{stem_base}.pdb"
            shutil.copy2(pdb_path, canonical)
            picnic_tmp = tmp / "picnic_out"
            score = score_with_picnic(canonical, uid, picnic_tmp)

        status = f"{score:.4f}" if score is not None else "FAILED"
        print(f"  {uid} {region} seg{seg_idx:02d} ({seg['start']}-{seg['end']}, "
              f"n={seg['n_residues']}) → {status}")

        rows.append({
            "UniProt_ID":  uid,
            "region":      region,
            "seg_idx":     seg_idx,
            "start":       seg["start"],
            "end":         seg["end"],
            "n_residues":  seg["n_residues"],
            "score":       score,
        })

    # Append new rows to any existing output
    all_rows = []
    if OUT_CSV.exists():
        all_rows = [pd.read_csv(OUT_CSV)]
    if rows:
        all_rows.append(pd.DataFrame(rows))

    if all_rows:
        out = pd.concat(all_rows, ignore_index=True)
        OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(OUT_CSV, index=False)
        print(f"\nWrote {len(out)} rows to {OUT_CSV}")
        print(out.groupby("region")["score"].describe().round(3).to_string())


if __name__ == "__main__":
    main()
