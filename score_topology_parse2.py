"""
score_topology_parse2.py

Run ParSe2 scoring on per-segment and concatenated region FASTAs.
Outputs output/topology_scores_raw/ParSe2_all_approaches.csv.

Score used: mean of non-NaN per-residue scores (normalised); also records
raw nansum. Sequences outside MIN_LEN/MAX_LEN limits get NaN scores.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent

# Import ParSe2 constants from score_parse2.py
sys.path.insert(0, str(ROOT))
from score_parse2 import (
    AA_PROPS, HYDR_F_CUTOFF, C1, C2, M, PS_DIST, WINDOW, MIN_LEN, MAX_LEN,
)

OUT_DIR = ROOT / "output" / "topology_scores_raw"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FASTA_DIR = ROOT / "output" / "topology_fasta"
SEG_DIR   = FASTA_DIR / "segments"

# ── FASTA sources ─────────────────────────────────────────────────────────────

SOURCES = [
    # (path, approach, region, is_segment)
    (SEG_DIR / "Cytoplasmic_segments.fasta",       "segment",     "Cytoplasmic",           True),
    (SEG_DIR / "Transmembrane_segments.fasta",      "segment",     "Transmembrane",         True),
    (SEG_DIR / "Extracellular_Lumenal_segments.fasta", "segment", "Extracellular/Lumenal", True),
    (FASTA_DIR / "Cytoplasmic.fasta",               "concatenated", "Cytoplasmic",          False),
    (FASTA_DIR / "Transmembrane.fasta",             "concatenated", "Transmembrane",        False),
    (FASTA_DIR / "Extracellular_Lumenal.fasta",     "concatenated", "Extracellular/Lumenal", False),
    (FASTA_DIR / "Whole.fasta",                     "whole",        "Whole",                False),
]

# Segment slug to region name mapping
SLUG_TO_REGION = {
    "cyto":  "Cytoplasmic",
    "tm":    "Transmembrane",
    "extra": "Extracellular/Lumenal",
}


# ── Per-residue ParSe2 scores ─────────────────────────────────────────────────

def parse2_residue_scores(seq, window=WINDOW):
    """Per-residue dist_norm_P (0.0 where not classified 'P'), NaN outside scored range."""
    L = len(seq)
    out = np.full(L, np.nan)
    if L < MIN_LEN or L > MAX_LEN:
        return out
    n_windows = L - window + 1
    if n_windows <= 0:
        return out

    seq = seq.upper()
    ppii_arr = np.zeros(L)
    helix_arr = np.zeros(L)
    hydr_arr = np.zeros(L)
    charge_arr = np.zeros(L)

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

    is_F       = hydr_w >= HYDR_F_CUTOFF
    ppii_eff   = np.where(ppii_w == 1.0, 0.98, ppii_w)
    v_exponent = 0.503 - 0.11 * np.log(1.0 - ppii_eff)
    rh         = 2.16 * (4 * window) ** v_exponent + 0.26 * 4 * net_charge_w - 0.29 * np.sqrt(4 * window)
    nu_model   = np.log(rh / 2.16) / np.log(4 * window)

    b = nu_model - M * helix_w
    x = (b - C2) / (C1 - M)
    y = M * x + b

    is_P        = (~is_F) & (((nu_model - C2) / C1) > helix_w)
    dist_norm_P = np.sqrt((helix_w - x) ** 2 + (nu_model - y) ** 2) / PS_DIST

    mid_offset = window // 2
    out[mid_offset: mid_offset + n_windows] = np.where(is_P, dist_norm_P, 0.0)
    return out


# ── FASTA parsing ─────────────────────────────────────────────────────────────

def parse_fasta_segment(path):
    """Parse segment FASTA: headers like >{uid}|{slug}|seg{i}.
    Returns list of (uid, region, seg_idx, seq)."""
    records = []
    uid = region = seg_idx = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if uid is not None:
                    records.append((uid, region, seg_idx, "".join(chunks)))
                header = line[1:]
                parts = header.split("|")
                uid     = parts[0].strip()
                slug    = parts[1].strip() if len(parts) > 1 else ""
                seg_str = parts[2].strip() if len(parts) > 2 else "seg0"
                region  = SLUG_TO_REGION.get(slug, slug)
                seg_idx = int(seg_str.replace("seg", "")) if seg_str.startswith("seg") else 0
                chunks = []
            else:
                chunks.append(line)
    if uid is not None:
        records.append((uid, region, seg_idx, "".join(chunks)))
    return records


def parse_fasta_plain(path, approach, region_name):
    """Parse concatenated/whole FASTA: headers like >{uid} or >sp|{uid}|...
    Returns list of (uid, approach, region, seg_idx=NaN, seq)."""
    records = []
    uid = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if uid is not None:
                    records.append((uid, approach, region_name, np.nan, "".join(chunks)))
                header = line[1:]
                parts = header.split("|")
                # Handle >sp|{uid}|... or plain >{uid}
                if len(parts) >= 3 and parts[0] in ("sp", "tr"):
                    uid = parts[1].strip()
                elif len(parts) >= 2:
                    uid = parts[1].strip()
                else:
                    uid = parts[0].strip()
                chunks = []
            else:
                chunks.append(line)
    if uid is not None:
        records.append((uid, approach, region_name, np.nan, "".join(chunks)))
    return records


# ── Main scoring loop ─────────────────────────────────────────────────────────

all_rows = []

for fasta_path, approach, region_name, is_segment in SOURCES:
    if not fasta_path.exists():
        print(f"WARNING: missing {fasta_path} — skipping")
        continue

    if is_segment:
        records = [(uid, approach, reg, idx, seq)
                   for uid, reg, idx, seq in parse_fasta_segment(fasta_path)]
    else:
        records = parse_fasta_plain(fasta_path, approach, region_name)

    n_scored = 0
    n_nan    = 0
    valid_scores = []

    for uid, appr, reg, seg_idx, seq in records:
        L   = len(seq)
        arr = parse2_residue_scores(seq)
        out_of_range = bool(L < MIN_LEN or L > MAX_LEN)

        if out_of_range:
            score_mean       = np.nan
            score_sum        = np.nan
            n_scored_res     = 0
            n_nan += 1
        else:
            non_nan_vals     = arr[~np.isnan(arr)]
            n_scored_res     = int(len(non_nan_vals))
            score_mean       = float(np.mean(non_nan_vals)) if n_scored_res > 0 else np.nan
            score_sum        = float(np.nansum(arr))
            n_scored += 1
            if not np.isnan(score_mean):
                valid_scores.append(score_mean)

        all_rows.append({
            "UniProt_ID":        uid,
            "approach":          appr,
            "region":            reg,
            "seg_idx":           seg_idx,
            "length":            L,
            "score_mean":        score_mean,
            "score_sum":         score_sum,
            "n_scored_residues": n_scored_res,
        })

    score_range = (
        f"{min(valid_scores):.4f} – {max(valid_scores):.4f}"
        if valid_scores else "N/A"
    )
    print(
        f"{fasta_path.name:<45}  "
        f"scored={n_scored:>5}  too_short_or_long={n_nan:>4}  "
        f"score_mean range: {score_range}"
    )

# ── Write output ──────────────────────────────────────────────────────────────

df = pd.DataFrame(all_rows, columns=[
    "UniProt_ID", "approach", "region", "seg_idx",
    "length", "score_mean", "score_sum", "n_scored_residues",
])

out_path = OUT_DIR / "ParSe2_all_approaches.csv"
df.to_csv(out_path, index=False)
print(f"\nWrote {len(df)} rows to {out_path}")
print(f"  Approaches: {df['approach'].unique().tolist()}")
print(f"  Regions:    {df['region'].unique().tolist()}")
print(f"  Total scored (non-NaN): {df['score_mean'].notna().sum()}")
print(f"  Total NaN (out-of-range): {df['score_mean'].isna().sum()}")
