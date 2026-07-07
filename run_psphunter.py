#!/usr/bin/env python
"""
Run PSPHunter on topology FASTAs.

PSPHunter uses a word2vec trigram model + 100 trained RandomForest classifiers.
This script replicates the Perl script logic in pure Python, avoiding the
15-char header truncation problem.

Usage:
    python run_psphunter.py

Outputs to output/topology_scores_raw/
"""
import os
import sys
import numpy as np
import joblib
import pandas as pd
from pathlib import Path

# ── paths ────────────────────────────────────────────────────────────────────
PSPHUNTER_DIR = Path("/home/jake/PSPHunter")
WORDVEC_FILE  = PSPHUNTER_DIR / "datasets/wordvec/uniprot_sprot70_size60.txt"
TRAINED_DIR   = PSPHUNTER_DIR / "Trained_model"
N_MODELS      = 100
VEC_SIZE      = 60

FASTA_DIR = Path("/home/jake/mem_prot_llps/mem_prot_llps/output/topology_fasta")
OUT_DIR   = Path("/home/jake/mem_prot_llps/mem_prot_llps/output/topology_scores_raw")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ── helpers ──────────────────────────────────────────────────────────────────

def load_wordvec(path: Path) -> dict:
    """Return {trigram: np.array(60,)} from the PSPHunter word2vec file."""
    wv = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("Uniprot"):
                continue
            parts = [p for p in line.split("\t") if p != ""]
            wv[parts[0]] = np.array(parts[1:], dtype=float)
    print(f"  Loaded {len(wv)} trigrams from word2vec file")
    return wv


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Return list of (header_without_>, sequence) pairs."""
    records = []
    header = None
    seq_parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts).upper()))
    return records


def seq_to_features(seq: str, wv: dict) -> np.ndarray:
    """Convert amino-acid sequence to 60-dim trigram word2vec feature vector."""
    vec = np.zeros(VEC_SIZE, dtype=float)
    for i in range(len(seq) - 2):
        tri = seq[i:i+3]
        if tri in wv:
            vec += wv[tri]
    return vec


def load_models(trained_dir: Path, n: int) -> list:
    """Load all n trained RandomForest models; return as list."""
    models = []
    for i in range(1, n + 1):
        model_path = trained_dir / str(i) / "word2vec70_60" / "train_model.m"
        models.append(joblib.load(model_path))
    return models


def predict_fasta(records: list[tuple[str,str]], wv: dict, models: list) -> list[tuple[str, float]]:
    """
    Returns (header, avg_prob) for each record.
    Empty sequences get score=None.
    """
    headers = [r[0] for r in records]
    feature_matrix = np.array([seq_to_features(r[1], wv) for r in records])  # (N, 60)

    # For empty sequences, features are zero – flag them
    empty_mask = np.array([len(r[1]) == 0 for r in records])

    # Accumulate probabilities across all 100 models
    prob_sum = np.zeros(len(records))
    for model in models:
        probs = model.predict_proba(feature_matrix)[:, 1]
        prob_sum += probs

    avg_probs = prob_sum / len(models)
    # Return None for empty sequences
    results = []
    for i, (h, _) in enumerate(records):
        score = None if empty_mask[i] else float(avg_probs[i])
        results.append((h, score))
    return results


# ── main run ─────────────────────────────────────────────────────────────────

def main():
    print("Loading word2vec model...")
    wv = load_wordvec(WORDVEC_FILE)

    print(f"Loading {N_MODELS} trained RandomForest models...")
    models = load_models(TRAINED_DIR, N_MODELS)
    print("  Models loaded.")

    # All rows for long-form output
    all_rows = []

    # ── Whole ──────────────────────────────────────────────────────────────
    fasta_path = FASTA_DIR / "Whole.fasta"
    print(f"\nProcessing {fasta_path.name}...")
    records = parse_fasta(fasta_path)
    print(f"  {len(records)} sequences")
    results = predict_fasta(records, wv, models)
    df = pd.DataFrame([
        {"UniProt_ID": h, "score": s} for h, s in results
    ])
    out_path = OUT_DIR / "PSPHunter_Whole_raw.csv"
    df.to_csv(out_path, index=False)
    print(f"  Saved -> {out_path}")
    for h, s in results:
        all_rows.append({"UniProt_ID": h, "approach": "concat", "region": "Whole", "seg_idx": np.nan, "score": s})

    # ── Concat (whole-protein per region) ─────────────────────────────────
    concat_regions = ["Cytoplasmic", "Transmembrane", "Extracellular_Lumenal"]
    for region in concat_regions:
        fasta_path = FASTA_DIR / f"{region}.fasta"
        print(f"\nProcessing {fasta_path.name}...")
        records = parse_fasta(fasta_path)
        print(f"  {len(records)} sequences")
        results = predict_fasta(records, wv, models)
        df = pd.DataFrame([{"UniProt_ID": h, "score": s} for h, s in results])
        out_name = f"PSPHunter_{region}_concat_raw.csv"
        out_path = OUT_DIR / out_name
        df.to_csv(out_path, index=False)
        print(f"  Saved -> {out_path}")
        for h, s in results:
            all_rows.append({"UniProt_ID": h, "approach": "concat", "region": region, "seg_idx": np.nan, "score": s})

    # ── Segments ───────────────────────────────────────────────────────────
    segment_regions = ["Cytoplasmic", "Transmembrane", "Extracellular_Lumenal"]
    for region in segment_regions:
        fasta_path = FASTA_DIR / "segments" / f"{region}_segments.fasta"
        print(f"\nProcessing {fasta_path.name}...")
        records = parse_fasta(fasta_path)
        print(f"  {len(records)} sequences")
        results = predict_fasta(records, wv, models)

        rows = []
        for h, s in results:
            parts = h.split("|")
            uid = parts[0]
            seg_idx = int(parts[2].replace("seg", "")) if len(parts) >= 3 else np.nan
            rows.append({"UniProt_ID": uid, "seg_idx": seg_idx, "score": s})
            all_rows.append({
                "UniProt_ID": uid,
                "approach": "segments",
                "region": region,
                "seg_idx": seg_idx,
                "score": s
            })

        df = pd.DataFrame(rows)
        out_name = f"PSPHunter_{region}_segments_raw.csv"
        out_path = OUT_DIR / out_name
        df.to_csv(out_path, index=False)
        print(f"  Saved -> {out_path}")

    # ── Long-form all approaches ───────────────────────────────────────────
    df_all = pd.DataFrame(all_rows)
    out_path = OUT_DIR / "PSPHunter_all_approaches.csv"
    df_all.to_csv(out_path, index=False)
    print(f"\nSaved long-form -> {out_path}")
    print(f"Total rows: {len(df_all)}")


if __name__ == "__main__":
    main()
