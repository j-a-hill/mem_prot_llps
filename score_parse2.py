"""
score_parse2.py

Phase 0d: ParSe v2 (Ibrahim et al., J. Biol. Chem. 2023) ported directly to
Python from the deployed web tool's source (stevewhitten.github.io/Parse_v2_FASTA,
js/parse.js + js/part-b.js). The tool has no trained model weights and no
external dependencies — it's deterministic per-residue biophysical math over
a sliding window (polymer scaling exponent vs. helix propensity, classified
against fixed reference points) — so it can be reproduced exactly rather
than relying on browser-based FASTA upload.

Score: PS potential = "Sigma classifier distance P" — the web tool's own
per-protein summary statistic, i.e. the sum of the per-residue classifier
distance over residues classified as phase-separation-driving ('P').

Validated against the actual JS engine (via Node) on FUS/TDP-43/p53: matches
to float precision (see conversation/PR for the harness).

Constraints carried over from the web tool: sequences <25 or >10,000
residues are not scored (NaN), matching its own stated limits.

Output: output/parse2_scores.csv (UniProt_ID, ParSe2_score)
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
PRED = ROOT / "Predictors_whole_genome_sets"
OUT  = ROOT / "output"

WINDOW = 25
MIN_LEN, MAX_LEN = 25, 10000

# ppii_h, helix, hydr propensities per residue (parse.js `amino_acids`).
# X -> ala defaults, U (selenocysteine) -> cys defaults, matching parse.js.
AA_PROPS = {
    "A": (0.37, 1.42, 0.0728),  "C": (0.25, 0.73, 0.3557),  "D": (0.30, 1.01, -0.0552),
    "E": (0.42, 1.63, -0.0295), "F": (0.17, 1.16, 0.4201),  "G": (0.13, 0.50, -0.0589),
    "H": (0.20, 1.20, 0.0874),  "I": (0.39, 1.12, 0.3805),  "K": (0.56, 1.24, -0.0053),
    "L": (0.24, 1.29, 0.3819),  "M": (0.36, 1.21, 0.1613),  "N": (0.27, 0.71, -0.0390),
    "P": (1.00, 0.65, -0.0492), "Q": (0.53, 1.02, 0.0126),  "R": (0.38, 1.06, 0.0394),
    "S": (0.24, 0.71, -0.0282), "T": (0.32, 0.78, 0.0239),  "V": (0.39, 0.99, 0.2947),
    "W": (0.25, 1.05, 0.4114),  "Y": (0.25, 0.67, 0.3113),
    "X": (0.37, 1.42, 0.0728),
    "U": (0.25, 0.73, 0.3557),
}

HYDR_F_CUTOFF = 0.08280152
C1, C2 = -0.244078945, 0.7885823
M = -1.0 / C1


def _dist_const(c3, c4):
    """Fixed reference-point distance (independent of sequence/window)."""
    b = c3 - M * c4
    x = (b - C2) / (C1 - M)
    y = M * x + b
    return np.sqrt((c4 - x) ** 2 + (c3 - y) ** 2)


PS_DIST = _dist_const(0.5416, 0.9327272)


def parse_fasta(path):
    seqs = {}
    uid, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if uid is not None:
                    seqs[uid] = "".join(chunks)
                parts = line[1:].split("|")
                uid = parts[1].strip() if len(parts) >= 2 else line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if uid is not None:
        seqs[uid] = "".join(chunks)
    return seqs


def score_protein(seq, window=WINDOW):
    """PS potential (Sigma classifier distance P) for one sequence."""
    L = len(seq)
    if L < MIN_LEN or L > MAX_LEN:
        return np.nan
    n_windows = L - window + 1
    if n_windows <= 0:
        return 0.0

    seq = seq.upper()
    ppii_arr   = np.zeros(L)
    helix_arr  = np.zeros(L)
    hydr_arr   = np.zeros(L)
    charge_arr = np.zeros(L)   # +1 for D/E, -1 for K/R (net_charge = |sum|)

    for i, ch in enumerate(seq):
        props = AA_PROPS.get(ch)
        if props is None:
            continue
        ppii_arr[i], helix_arr[i], hydr_arr[i] = props
        if ch in ("D", "E"):
            charge_arr[i] = 1.0
        elif ch in ("K", "R"):
            charge_arr[i] = -1.0

    def window_sums(arr):
        c = np.concatenate(([0.0], np.cumsum(arr)))
        return c[window:] - c[:-window]

    ppii_w       = window_sums(ppii_arr) / window
    helix_w      = window_sums(helix_arr) / window
    hydr_w       = window_sums(hydr_arr) / window
    net_charge_w = np.abs(window_sums(charge_arr))

    is_F = hydr_w >= HYDR_F_CUTOFF

    ppii_eff   = np.where(ppii_w == 1.0, 0.98, ppii_w)
    v_exponent = 0.503 - 0.11 * np.log(1.0 - ppii_eff)
    rh         = 2.16 * (4 * window) ** v_exponent + 0.26 * 4 * net_charge_w - 0.29 * np.sqrt(4 * window)
    nu_model   = np.log(rh / 2.16) / np.log(4 * window)

    b = nu_model - M * helix_w
    x = (b - C2) / (C1 - M)
    y = M * x + b

    is_P = (~is_F) & (((nu_model - C2) / C1) > helix_w)
    dist_norm_P = np.sqrt((helix_w - x) ** 2 + (nu_model - y) ** 2) / PS_DIST

    return float(np.sum(np.where(is_P, dist_norm_P, 0.0)))


def main():
    fasta_path = PRED / "Seq2Phase_swiss_prot_human_220916.fasta"
    seqs = parse_fasta(fasta_path)
    print(f"Loaded {len(seqs)} sequences")

    rows = []
    for i, (uid, seq) in enumerate(seqs.items()):
        rows.append({"UniProt_ID": uid, "ParSe2_score": score_protein(seq)})
        if (i + 1) % 5000 == 0:
            print(f"  {i + 1}/{len(seqs)}")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "parse2_scores.csv", index=False)
    print(f"Scored {out['ParSe2_score'].notna().sum()} / {len(out)} proteins")
    print(f"Saved {OUT / 'parse2_scores.csv'}")


if __name__ == "__main__":
    main()
