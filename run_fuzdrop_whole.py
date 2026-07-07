"""
run_fuzdrop_whole.py

Run FuzDrop on the whole-protein FASTA (Whole.fasta) for the 60 membrane
LLPS proteins, using each protein's full ESpritz-DisProt residue array from
the PhaSePred JSON.

Output:
  output/topology_scores_raw/FuzDrop_whole_raw.csv
    columns: UniProt_ID, length, fuzdrop_score
"""

import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

ROOT           = Path(__file__).parent
PRED           = ROOT / "Predictors_whole_genome_sets"
OUT            = ROOT / "output"
RAW_OUT        = OUT / "topology_scores_raw"
FASTA_DIR      = OUT / "topology_fasta"
PHASEPRED_JSON = PRED / "PhaSePred_human_reviewed.json"
FUZDROP_BIN    = ROOT / "external_tools" / "FuzDrop" / "FuzDrop"
WHOLE_FASTA    = FASTA_DIR / "Whole.fasta"
OUT_CSV        = RAW_OUT / "FuzDrop_whole_raw.csv"

RAW_OUT.mkdir(parents=True, exist_ok=True)

ESPRITZ_HEADER = (
    "*" * 108 + "\n"
    "Licensed to: Prof. MONIKA FUXREITER (University of Padova) located in Padova, Italy. \n"
    "This license is for non-commercial use only. Please see LICENSE file for details \n"
    "(https://biocomputingup.it/assets/data/LICENSE.txt)\n\n"
    "Contact silvio.tosatto@unipd.it for commercial licensing details.\n\n\n"
    + "*" * 108 + "\n"
)


def parse_fasta(fasta_path):
    seqs = {}
    uid_cur, chunks = None, []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if uid_cur is not None:
                    seqs[uid_cur] = "".join(chunks)
                uid_cur = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if uid_cur is not None:
        seqs[uid_cur] = "".join(chunks)
    return seqs


def make_espritz_content(labels_iter, scores_iter):
    body = "\n".join(f"{l}\t{s}" for l, s in zip(labels_iter, scores_iter))
    return ESPRITZ_HEADER + body + "\n"


def run_fuzdrop(stem, fasta_content, espritz_content, work_dir):
    fasta_path   = work_dir / f"{stem}.fasta"
    espritz_path = work_dir / f"{stem}.espritz"
    res_path     = work_dir / f"{stem}_res.txt"

    fasta_path.write_text(fasta_content)
    espritz_path.write_text(espritz_content)

    try:
        result = subprocess.run(
            [str(FUZDROP_BIN.resolve()), fasta_path.name, espritz_path.name],
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=120,
        )
    except subprocess.TimeoutExpired:
        return None
    finally:
        fasta_path.unlink(missing_ok=True)
        espritz_path.unlink(missing_ok=True)

    score = None
    if result.returncode == 0 and res_path.exists():
        last_line = res_path.read_text().strip().splitlines()[-1]
        try:
            score = float(last_line.split("=")[-1].strip())
        except ValueError:
            pass
        res_path.unlink(missing_ok=True)

    return score


print("Loading PhaSePred JSON (737 MB) — please wait...")
with open(PHASEPRED_JSON) as f:
    phasepred = json.load(f)
print(f"  Loaded {len(phasepred):,} entries.")

whole_seqs = parse_fasta(WHOLE_FASTA)
print(f"\nWhole.fasta: {len(whole_seqs)} sequences")

work_dir = Path(tempfile.mkdtemp(prefix="fuzdrop_whole_"))
print(f"Working dir: {work_dir}\n")

rows = []
n_ok = n_skip = n_fail = 0

for i, (uid, seq) in enumerate(whole_seqs.items()):
    if i % 10 == 0:
        print(f"  [{i}/{len(whole_seqs)}] {uid}  len={len(seq)}")

    esp = phasepred.get(uid, {}).get("ESpritz-DisProt", {})
    esp_labels = esp.get("label", "").split(",")
    esp_scores = esp.get("residue", "").split(",")
    full_seq = phasepred.get(uid, {}).get("Sequence", "")

    if not full_seq or len(esp_labels) != len(full_seq) or len(esp_scores) != len(full_seq):
        print(f"    [skip] {uid}: no ESpritz data")
        n_skip += 1
        continue

    if seq != full_seq:
        print(f"    [warn] {uid}: Whole.fasta seq != PhaSePred seq "
              f"(fasta={len(seq)}, json={len(full_seq)}) — using JSON seq")
        seq = full_seq

    fasta_content   = f">{uid}\n{seq}\n"
    espritz_content = make_espritz_content(esp_labels, esp_scores)

    score = run_fuzdrop(uid, fasta_content, espritz_content, work_dir)

    if score is not None:
        rows.append({"UniProt_ID": uid, "length": len(seq), "fuzdrop_score": score})
        n_ok += 1
    else:
        print(f"    [FAIL] {uid}")
        n_fail += 1

shutil.rmtree(work_dir, ignore_errors=True)

df = pd.DataFrame(rows)
df.to_csv(OUT_CSV, index=False)
print(f"\nDone: {n_ok} scored, {n_fail} failed, {n_skip} skipped")
print(f"Saved {len(df)} rows to {OUT_CSV.relative_to(ROOT)}")
