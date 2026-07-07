"""
Run PSPsPredict on per-segment topology FASTA files.
Mirrors run_pspspredict.py but operates on the _segments FASTAs.

Expected layout (WSL paths):
  ~/PSPsPredict/src/                      -- Predict.py, embedding.py
  ~/PSPsPredict/model/                    -- best-0.901.pth
  ~/PSPsPredict/input/*_segments.fasta    -- placed here before running
  ~/PSPsPredict/embeddings/               -- created here
  ~/PSPsPredict/output/*_segments_scores.csv  -- results
"""

import csv
import subprocess
from pathlib import Path

BASE    = Path.home() / "PSPsPredict"
SRC     = BASE / "src"
INPUT   = BASE / "input"
EMB_DIR = BASE / "embeddings"
OUT_DIR = BASE / "output"
PYTHON  = Path.home() / "miniconda3/envs/PSPsPredict/bin/python"
MODEL_ID = "Rostlab/prot_t5_xl_half_uniref50-enc"

EMB_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

REGIONS = ["Cytoplasmic_segments", "Transmembrane_segments", "Extracellular_Lumenal_segments"]


def read_fasta(path):
    records, header, chunks = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:].strip(), []
        else:
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


for region in REGIONS:
    fasta = INPUT / f"{region}.fasta"
    if not fasta.exists():
        print(f"[skip] {fasta} not found — copy it to {INPUT} first")
        continue

    records = read_fasta(fasta)
    print(f"\n=== {region}: {len(records)} sequences ===")

    out_csv = OUT_DIR / f"{region}_scores.csv"
    if out_csv.exists():
        print(f"  Already done ({out_csv.name}), skipping.")
        continue

    csv_path = BASE / f"{region}.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "seq"])
        for name, seq in records:
            w.writerow([name, seq])

    emb_subdir = EMB_DIR / region
    emb_subdir.mkdir(exist_ok=True)

    print(f"  Generating embeddings -> {emb_subdir}")
    r = subprocess.run(
        [str(PYTHON), str(SRC / "embedding.py"),
         "-i", str(csv_path),
         "-o", str(emb_subdir) + "/",
         "--model", MODEL_ID],
        cwd=str(SRC)
    )
    if r.returncode != 0:
        print(f"  [ERROR] embedding failed for {region}")
        continue

    print(f"  Running predictions -> {out_csv}")
    r = subprocess.run(
        [str(PYTHON), str(SRC / "Predict.py"),
         "-i", str(csv_path),
         "-src", str(emb_subdir) + "/",
         "-o", str(out_csv)],
        cwd=str(SRC)
    )
    if r.returncode != 0:
        print(f"  [ERROR] predict failed for {region}")
        continue

    print(f"  Done -> {out_csv}")

print("\nAll segment regions complete.")
