"""
run_deephase_topology.py

Run DeePhase on all 7 topology FASTAs (whole, 3 concat, 3 segments).
Must be run from the mem_prot_llps directory with the DeePhase conda env.

Usage:
  /home/jake/miniconda3/bin/conda run -n DeePhase python run_deephase_topology.py

Outputs (topology_scores_raw/):
  DeePhase_Whole_raw.csv
  DeePhase_{Cytoplasmic,Transmembrane,Extracellular_Lumenal}_concat_raw.csv
  DeePhase_{Cytoplasmic,Transmembrane,Extracellular_Lumenal}_segments_raw.csv
  DeePhase_all_approaches.csv
"""

import os
import sys
import pickle
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from pathlib import Path

ROOT           = Path(__file__).parent
DEEPHASE_DIR   = Path("/home/jake/DeePhase/__PREDICT")
RAW_OUT        = ROOT / "output/topology_scores_raw"
FASTA_DIR      = ROOT / "output/topology_fasta"
RAW_OUT.mkdir(parents=True, exist_ok=True)

# DeePhase must run from its own directory (relative model paths)
os.chdir(DEEPHASE_DIR)
sys.path.insert(0, str(DEEPHASE_DIR))

from gensim.models import word2vec as gensim_w2v

# ProtVec must exist in the module scope that pickle was saved from.
# The saved model references __main__.ProtVec, so we define it here
# (this script IS __main__) before calling Word2Vec.load().
class ProtVec(gensim_w2v.Word2Vec):
    """Stub — only needed so pickle can reconstruct the saved model."""
    def __init__(self, fasta_fname=None, corpus=None, n=3, size=100,
                 corpus_fname="corpus.txt", sg=1, window=25, min_count=1, workers=20):
        self.n = n
        self.size = size
        self.fasta_fname = fasta_fname

def load_protvec(model_fname):
    return gensim_w2v.Word2Vec.load(model_fname)

def split_ngrams(seq, n):
    a = zip(*[iter(seq)]*n)
    b = zip(*[iter(seq[1:])]*n)
    c = zip(*[iter(seq[2:])]*n)
    result = []
    for ngrams in [a, b, c]:
        result.append(["".join(ng) for ng in ngrams])
    return result

def get_vector(pv, seq):
    ngram_patterns = split_ngrams(seq, 3)
    vecs = []
    for ngrams in ngram_patterns:
        ngram_sum = None
        for ngram in ngrams:
            try:
                v = pv.wv[ngram]
            except KeyError:
                v = np.zeros(pv.wv.vector_size)
            ngram_sum = v if ngram_sum is None else ngram_sum + v
        if ngram_sum is not None:
            vecs.append(ngram_sum)
    return sum(vecs) if vecs else np.zeros(pv.wv.vector_size)

print("Loading word2vec model...")
pv = load_protvec("tools/Embeddings/swissprot_size200_window25.model")
print(f"  vector_size = {pv.wv.vector_size}")

print("Loading RF models...")
phys_model = pickle.load(open("tools/Models/phys_multi.sav", "rb"))
w2v_model  = pickle.load(open("tools/Models/w2v_multi.sav",  "rb"))

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
import math, subprocess, tempfile

AA_array = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
kd = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,"H":-3.2,
      "I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}

Hydrophobic_AAs = ['A','I','L','M','F','V']
Polar_AAs       = ['S','Q','N','G','C','T','P']
Cation_AAs      = ['K','R','H']
Anion_AAs       = ['D','E']
Arom_AAs        = ['W','Y','F']

def hydrophobicity(seq):
    s = ProteinAnalysis(seq)
    return sum(s.count_amino_acids()[aa] * kd[aa] for aa in AA_array)

def shannon_entropy(seq):
    s = ProteinAnalysis(seq)
    pct = s.get_amino_acids_percent()
    return sum(-math.log2(p) * p for aa in AA_array if (p := pct[aa]) > 0)

def extract_lcr(seq):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(f'>1\n{seq}')
        fname = f.name
    try:
        out = subprocess.run(["segmasker", "-in", fname],
                             capture_output=True, timeout=30)
        tokens = out.stdout.split()[1:]
        starts, ends = [], []
        for i in range(0, len(tokens) // 3):
            starts.append(int(tokens[3*i]))
            ends.append(int(tokens[3*i+2]))
    except Exception:
        starts, ends = [], []
    finally:
        os.unlink(fname)
    residues = set()
    for s, e in zip(starts, ends):
        residues.update(range(s, e+1))
    lcr_seq = "".join(seq[r-1] for r in sorted(residues) if 0 < r <= len(seq))
    return len(residues), lcr_seq

def extract_idr(seq):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(f'>1\n{seq}')
        fname = f.name
    try:
        out = subprocess.run(
            ["python", "tools/iupred2a.py", fname, "long"],
            capture_output=True, timeout=60
        )
        tokens = out.stdout.split()[40:]
        probs = [float(tokens[3*i+2]) for i in range(len(tokens)//3)]
    except Exception:
        probs = [0.0] * len(seq)
    finally:
        os.unlink(fname)
    TH1, TH2 = 0.5, 20
    idr_residues = []
    current = 0
    for t, p in enumerate(probs):
        if p > TH1:
            current += 1
            if t == len(probs) - 1 and current > TH2:
                idr_residues.extend(range(t - current, t + 1))
        else:
            if current > TH2:
                idr_residues.extend(range(t - current, t))
            current = 0
    return len(idr_residues)

def get_aa_count(seq, aa):
    return ProteinAnalysis(seq).count_amino_acids()[aa] if seq and not isinstance(seq, float) else 0


def score_sequence(seq):
    """Return DeePhase score (float) for a single sequence string."""
    seq = seq.upper().replace("*", "").replace("-", "")
    # strip non-standard AAs
    seq_clean = "".join(c for c in seq if c in AA_array)
    if len(seq_clean) < 10:
        return None

    # --- physical features ---
    lcr_len, lcr_seq = extract_lcr(seq_clean)
    idr_len = extract_idr(seq_clean)
    seqlen = len(seq_clean)

    hydro  = hydrophobicity(seq_clean)
    shanno = shannon_entropy(seq_clean)
    lcr_f  = lcr_len / seqlen if seqlen > 0 else 0
    idr_f  = idr_len / seqlen if seqlen > 0 else 0
    arom_f = sum(get_aa_count(seq_clean, aa) for aa in Arom_AAs) / seqlen
    cat_f  = sum(get_aa_count(seq_clean, aa) for aa in Cation_AAs) / seqlen

    # predict_multiclass expects features in alphabetical column order
    phys_row = pd.DataFrame([{
        "Arom_frac": arom_f, "Cation_frac": cat_f,
        "Hydrophobicity": hydro, "IDR_frac": idr_f,
        "LCR_frac": lcr_f, "Shannon_entropy": shanno,
        "sequence_final": seq_clean,
    }])
    phys_row = phys_row.reindex(sorted(phys_row.columns), axis=1)
    phys_X = phys_row.drop(columns=["sequence_final"]).to_numpy()
    phys_score = phys_model.predict_proba(phys_X)[0, 0] + 0.5 * phys_model.predict_proba(phys_X)[0, 1]

    # --- w2v features ---
    vec = get_vector(pv, seq_clean)
    vec_size = pv.wv.vector_size
    w2v_row = pd.DataFrame([dict(zip([str(i) for i in range(vec_size)], vec))])
    w2v_row["sequence_final"] = seq_clean
    w2v_row = w2v_row.reindex(sorted(w2v_row.columns), axis=1)
    w2v_X = w2v_row.drop(columns=["sequence_final"]).to_numpy()
    w2v_score = w2v_model.predict_proba(w2v_X)[0, 0] + 0.5 * w2v_model.predict_proba(w2v_X)[0, 1]

    return 0.5 * (phys_score + w2v_score)


def parse_fasta(fasta_path):
    seqs = {}
    uid, chunks = None, []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if uid:
                    seqs[uid] = "".join(chunks)
                uid = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if uid:
        seqs[uid] = "".join(chunks)
    return seqs


FASTAS = {
    "Whole":                         (FASTA_DIR / "Whole.fasta",                               "whole",        "Whole"),
    "Cytoplasmic_concat":            (FASTA_DIR / "Cytoplasmic.fasta",                         "concatenated", "Cytoplasmic"),
    "Transmembrane_concat":          (FASTA_DIR / "Transmembrane.fasta",                       "concatenated", "Transmembrane"),
    "Extracellular_Lumenal_concat":  (FASTA_DIR / "Extracellular_Lumenal.fasta",               "concatenated", "Extracellular/Lumenal"),
    "Cytoplasmic_segments":          (FASTA_DIR / "segments/Cytoplasmic_segments.fasta",        "segment",      "Cytoplasmic"),
    "Transmembrane_segments":        (FASTA_DIR / "segments/Transmembrane_segments.fasta",      "segment",      "Transmembrane"),
    "Extracellular_Lumenal_segments":(FASTA_DIR / "segments/Extracellular_Lumenal_segments.fasta","segment",   "Extracellular/Lumenal"),
}

all_rows = []

for fasta_key, (fasta_path, approach, region) in FASTAS.items():
    out_stem = f"DeePhase_{fasta_key.replace('/', '_')}_raw.csv"
    out_path = RAW_OUT / out_stem

    if out_path.exists():
        print(f"\n[SKIP] {out_stem} already exists")
        existing = pd.read_csv(out_path)
        for _, r in existing.iterrows():
            seg_idx = int(r["seg_idx"]) if "seg_idx" in r and pd.notna(r.get("seg_idx")) else None
            all_rows.append({"UniProt_ID": r["UniProt_ID"], "approach": approach,
                             "region": region, "seg_idx": seg_idx, "score": r["score"]})
        continue

    seqs = parse_fasta(fasta_path)
    print(f"\n[{fasta_key}] {len(seqs)} sequences from {fasta_path.name}")
    rows = []
    for j, (header, seq) in enumerate(seqs.items()):
        # Extract UniProt_ID and seg_idx from header
        parts = header.split("|")
        uid = parts[0]
        seg_idx = None
        if len(parts) >= 3 and parts[2].startswith("seg"):
            try:
                seg_idx = int(parts[2][3:])
            except ValueError:
                pass

        if j % 10 == 0:
            print(f"  [{j}/{len(seqs)}] {uid}  len={len(seq)}")

        score = score_sequence(seq)
        if score is not None:
            rows.append({"UniProt_ID": uid, "seg_idx": seg_idx, "score": score})
        else:
            print(f"    [skip] {uid} — sequence too short")

    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    print(f"  Saved {len(df)} rows to {out_stem}")
    for _, r in df.iterrows():
        all_rows.append({"UniProt_ID": r["UniProt_ID"], "approach": approach,
                         "region": region, "seg_idx": r.get("seg_idx"), "score": r["score"]})

# Build combined file
combined = pd.DataFrame(all_rows)
combined_path = RAW_OUT / "DeePhase_all_approaches.csv"
combined.to_csv(combined_path, index=False)
print(f"\n=== DeePhase done: {len(combined)} total rows ===")
print(combined.groupby(["approach","region"])["UniProt_ID"].count().to_string())
