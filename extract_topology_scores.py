"""
extract_topology_scores.py

Does the cytoplasmic / transmembrane / extracellular-or-luminal region of a
membrane protein drive its whole-protein LLPS score — or at least score
higher than the whole-protein value? Supersedes the earlier cytoplasmic-only,
2-tool version (extract_cytoplasmic_scores.py) with a full 3-way topology
split across every tool that actually has per-residue resolution.

For each of the 60 experimentally-confirmed membrane LLPS proteins:
  - Fetches UniProt topological-domain + transmembrane features (cached to
    data/uniprot_topology_cache.csv), bucketed into:
      Cytoplasmic             UniProt "Cytoplasmic"
      Extracellular/Lumenal   UniProt "Extracellular", "Lumenal", etc. —
                              topologically the same "outside" face, so
                              pooled together
      Transmembrane           the TM-helix spans themselves
  - Extracts per-residue arrays for the 6 tools with residue-level
    resolution: catGRANULE, PLAAC, PScore, ESpritz-DisProt (continuous) and
    SEG (binary low-complexity-region label) from PhaSePred_human_reviewed
    .json; ParSe2 computed locally (reusing score_parse2.py's algorithm).
  - Computes each region's mean score and compares it, paired per protein,
    against the whole-protein score.

Run extract_topology_fasta.py afterwards — it adds R+Y and LLPhyScore (both
locally runnable, no per-residue PhaSePred data needed) and regenerates the
delta summary/figure with all 8 tools, plus a manual-retrieval checklist for
the predictors that have neither.

Outputs
-------
  output/topology_domain_scores.csv   one row per protein, wide: whole +
                                       per-region mean per tool
  output/topology_score_deltas.csv    per tool x region paired-comparison summary
  output/figures/topology_score_deltas.png
"""

import json
import re
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

from score_parse2 import (
    AA_PROPS, HYDR_F_CUTOFF, C1, C2, M, PS_DIST, WINDOW, MIN_LEN, MAX_LEN,
)

ROOT       = Path(__file__).parent
PRED       = ROOT / "Predictors_whole_genome_sets"
OUT        = ROOT / "output"
FIG_DIR    = OUT / "figures"
TOPO_CACHE = ROOT / "data" / "uniprot_topology_cache.csv"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"

FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Load the 60 positives + sequences ──────────────────────────────────────

pos_ids = pd.read_csv(ROOT / "exp_db" / "membrane_exp_db_matches.csv")["Entry"].tolist()
print(f"Positive proteins: {len(pos_ids)}")

seqs = {}
uid, chunks = None, []
with open(PRED / "Seq2Phase_swiss_prot_human_220916.fasta") as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if uid is not None:
                seqs[uid] = "".join(chunks)
            parts = line[1:].split("|")
            uid = parts[1].strip() if len(parts) >= 2 else None
            chunks = []
        else:
            chunks.append(line)
    if uid is not None:
        seqs[uid] = "".join(chunks)

# ── Fetch + cache UniProt topology ──────────────────────────────────────────

def fetch_topology_raw(entry_ids, cache_path):
    if cache_path.exists():
        print(f"Loading cache: {cache_path}")
        return pd.read_csv(cache_path)
    print(f"Fetching topology for {len(entry_ids)} proteins...")
    params = {
        "query": "accession:(" + " OR ".join(entry_ids) + ")",
        "fields": "accession,ft_topo_dom,ft_transmem",
        "format": "tsv",
        "size": 500,
    }
    resp = requests.get(UNIPROT_URL, params=params, timeout=30)
    resp.raise_for_status()
    result = pd.read_csv(StringIO(resp.text), sep="\t")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(cache_path, index=False)
    return result

topo_raw = fetch_topology_raw(pos_ids, TOPO_CACHE)
topo_raw = topo_raw.rename(columns={
    "Entry": "UniProt_ID", "Topological domain": "topo_dom_raw", "Transmembrane": "transmem_raw",
})
topo_raw = topo_raw.set_index("UniProt_ID")

FEATURE_RE = re.compile(r'(TOPO_DOM|TRANSMEM)\s+(\d+)\.\.(\d+);\s*/note="([^"]*)"')


def parse_features(raw_text):
    if pd.isna(raw_text):
        return []
    return [(m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))
            for m in FEATURE_RE.finditer(str(raw_text))]


def bucket_topo_dom(note):
    return "Cytoplasmic" if "ytoplasmic" in note else "Extracellular/Lumenal"


def build_region_labels(uid, length):
    """Per-residue region label array (index 0 = residue 1), default None (unannotated)."""
    labels = np.full(length, None, dtype=object)
    if uid not in topo_raw.index:
        return labels
    row = topo_raw.loc[uid]
    for _, s, e, note in parse_features(row["topo_dom_raw"]):
        labels[s - 1:min(e, length)] = bucket_topo_dom(note)
    for _, s, e, note in parse_features(row["transmem_raw"]):
        labels[s - 1:min(e, length)] = "Transmembrane"
    return labels


def build_segments(labels, region):
    """Return a list of (start, end) 0-based, end-exclusive index pairs for each
    contiguous run of `region` in the label array."""
    segments = []
    n = len(labels)
    i = 0
    while i < n:
        if labels[i] == region:
            start = i
            while i < n and labels[i] == region:
                i += 1
            segments.append((start, i))
        else:
            i += 1
    return segments


# ── PhaSePred residue-level scores ──────────────────────────────────────────

print("Loading PhaSePred JSON (737 MB, this may take a moment)...")
with open(PRED / "PhaSePred_human_reviewed.json") as f:
    phasepred = json.load(f)


def parse_residue_array(uid, tool, subkey):
    entry = phasepred.get(uid, {})
    s = entry.get(tool, {}).get(subkey, "")
    if not s:
        return None
    try:
        return np.array([float(x) for x in str(s).split(",") if x.strip() != ""])
    except ValueError:
        return None


def whole_score(uid, tool, subkey):
    try:
        return float(phasepred[uid][tool][subkey])
    except (KeyError, TypeError, ValueError):
        return np.nan


# tool label -> (PhaSePred key, residue subkey, whole-score subkey)
PHASEPRED_TOOLS = {
    "catGRANULE": ("catGRANULE", "residue", "single"),
    "PLAAC":      ("PLAAC", "residue", "NLLR"),
    "PScore":     ("PScore", "residue", "single"),
    "ESpritz":    ("ESpritz-DisProt", "residue", "single"),
    "SEG":        ("SEG", "label", "single"),
}

# ── ParSe2, computed locally (mirrors score_parse2.score_protein) ──────────

def parse2_residue_scores(seq, window=WINDOW):
    """Per-residue dist_norm_P (0.0 where not classified 'P'), NaN outside scoreable range."""
    L = len(seq)
    out = np.full(L, np.nan)
    if L < MIN_LEN or L > MAX_LEN:
        return out
    n_windows = L - window + 1
    if n_windows <= 0:
        return out

    seq = seq.upper()
    ppii_arr, helix_arr, hydr_arr, charge_arr = (np.zeros(L) for _ in range(4))
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

    is_P        = (~is_F) & (((nu_model - C2) / C1) > helix_w)
    dist_norm_P = np.sqrt((helix_w - x) ** 2 + (nu_model - y) ** 2) / PS_DIST

    mid_offset = window // 2
    out[mid_offset: mid_offset + n_windows] = np.where(is_P, dist_norm_P, 0.0)
    return out


# ── Per-protein, per-tool, per-region means ─────────────────────────────────

REGIONS = ["Cytoplasmic", "Transmembrane", "Extracellular/Lumenal"]


def region_mean(arr, labels, bucket=None):
    """Mean over residues labelled `bucket`, or the whole array if bucket is None.
    Always the same aggregation (mean of the per-residue array) so whole-vs-region
    comparisons are apples-to-apples — NOT each tool's own precomputed single-value
    summary, which may be a sum, max, or separately-calibrated statistic."""
    if arr is None:
        return np.nan
    if bucket is None:
        vals = arr
    else:
        n = min(len(arr), len(labels))
        mask = labels[:n] == bucket
        vals = arr[:n][mask]
    vals = vals[~np.isnan(vals.astype(float))]
    return float(np.mean(vals)) if len(vals) else np.nan


def segment_stats_for_arr(arr, segments, n_arr):
    """Compute per-segment score statistics for a score array.

    Parameters
    ----------
    arr      : np.ndarray or None — per-residue scores
    segments : list of (start, end) pairs from build_segments()
    n_arr    : length of arr (used for boundary clipping)

    Returns a dict with keys:
      mean_seg_score    unweighted mean of per-segment means (each segment equal weight)
      longest_seg_score mean score of the single longest segment
      max_seg_score     highest per-segment mean (best segment)
      seg_lengths       comma-separated segment lengths  (diagnostic)
      seg_scores        comma-separated per-segment mean scores (same order)
    """
    _nan_result = {
        "mean_seg_score": np.nan,
        "longest_seg_score": np.nan,
        "max_seg_score": np.nan,
        "seg_lengths": "",
        "seg_scores": "",
    }
    if arr is None or len(segments) == 0:
        return _nan_result

    seg_lens = []
    seg_means = []
    for start, end in segments:
        end_clipped = min(end, n_arr)
        if start >= n_arr:
            seg_lens.append(end - start)
            seg_means.append(np.nan)
            continue
        vals = arr[start:end_clipped].astype(float)
        vals = vals[~np.isnan(vals)]
        seg_lens.append(end - start)
        seg_means.append(float(np.mean(vals)) if len(vals) else np.nan)

    valid_means = [m for m in seg_means if not np.isnan(m)]
    mean_seg_score = float(np.mean(valid_means)) if valid_means else np.nan
    max_seg_score  = float(np.max(valid_means))  if valid_means else np.nan

    longest_idx = int(np.argmax(seg_lens)) if seg_lens else 0
    longest_seg_score = seg_means[longest_idx] if seg_lens else np.nan

    seg_lengths_str = ",".join(str(l) for l in seg_lens)
    seg_scores_str  = ",".join(f"{m:.4f}" if not np.isnan(m) else "nan" for m in seg_means)

    return {
        "mean_seg_score": mean_seg_score,
        "longest_seg_score": longest_seg_score,
        "max_seg_score": max_seg_score,
        "seg_lengths": seg_lengths_str,
        "seg_scores": seg_scores_str,
    }


rows = []
for uid in pos_ids:
    seq = seqs.get(uid)
    length = len(seq) if seq else None
    labels = build_region_labels(uid, length) if length else np.array([], dtype=object)
    n_annotated = int((labels != None).sum()) if length else 0  # noqa: E711

    row = {
        "UniProt_ID": uid,
        "length": length,
        "n_annotated": n_annotated,
        "frac_annotated": round(n_annotated / length, 4) if length else np.nan,
    }
    for region in REGIONS:
        row[f"n_{region}"] = int((labels == region).sum()) if length else 0

    # Segment counts are topology-only (not tool-specific) — compute once per region.
    region_segments = {r: (build_segments(labels, r) if length else []) for r in REGIONS}
    for region in REGIONS:
        row[f"n_segs_{region}"] = len(region_segments[region])

    for label, (tool, subkey, whole_key) in PHASEPRED_TOOLS.items():
        arr = parse_residue_array(uid, tool, subkey)
        n_arr = len(arr) if arr is not None else 0
        row[f"{label}_whole_official"] = whole_score(uid, tool, whole_key)
        row[f"{label}_whole"] = region_mean(arr, labels) if arr is not None and length else np.nan
        for region in REGIONS:
            row[f"{label}_{region}"] = region_mean(arr, labels, region) if arr is not None and length else np.nan
            stats = segment_stats_for_arr(arr if length else None, region_segments[region], n_arr)
            for k, v in stats.items():
                row[f"{label}_{region}_{k}"] = v

    if seq:
        p2_arr = parse2_residue_scores(seq)
        n_arr = len(p2_arr)
        row["ParSe2_whole_official"] = float(np.nansum(p2_arr)) if not np.all(np.isnan(p2_arr)) else np.nan
        row["ParSe2_whole"] = region_mean(p2_arr, labels)
        for region in REGIONS:
            row[f"ParSe2_{region}"] = region_mean(p2_arr, labels, region)
            stats = segment_stats_for_arr(p2_arr, region_segments[region], n_arr)
            for k, v in stats.items():
                row[f"ParSe2_{region}_{k}"] = v
    else:
        row["ParSe2_whole_official"] = np.nan
        row["ParSe2_whole"] = np.nan
        for region in REGIONS:
            row[f"ParSe2_{region}"] = np.nan
            stats = segment_stats_for_arr(None, [], 0)
            for k, v in stats.items():
                row[f"ParSe2_{region}_{k}"] = v

    rows.append(row)

topo_df = pd.DataFrame(rows)
topo_df.to_csv(OUT / "topology_domain_scores.csv", index=False)
print(f"\nSaved output/topology_domain_scores.csv  ({len(topo_df)} proteins)")

n_any_topo = (topo_df["n_annotated"] > 0).sum()
print(f"{n_any_topo}/{len(topo_df)} proteins have >=1 annotated topology residue")
for region in REGIONS:
    n_region = (topo_df[f"n_{region}"] > 0).sum()
    print(f"  {region:<24} annotated in {n_region}/{len(topo_df)} proteins")

# Cytoplasmic segment summary — shows whether the per-segment analysis is meaningful.
cyto_segs = topo_df["n_segs_Cytoplasmic"]
n_multi  = int((cyto_segs > 1).sum())
n_single = int((cyto_segs == 1).sum())
n_zero   = int((cyto_segs == 0).sum())
print(f"\nCytoplasmic segment breakdown (n_segs_Cytoplasmic):")
print(f"  >1 segment  : {n_multi} proteins  ← multi-loop, per-segment analysis meaningful")
print(f"   1 segment  : {n_single} proteins")
print(f"   0 segments : {n_zero} proteins  (no annotated cytoplasmic residues)")

# ── Paired comparison: region-mean vs whole-protein score ──────────────────

TOOLS = list(PHASEPRED_TOOLS.keys()) + ["ParSe2"]

summary_rows = []
for tool in TOOLS:
    for region in REGIONS:
        sub = topo_df[[f"{tool}_whole", f"{tool}_{region}"]].dropna()
        n = len(sub)
        if n < 3:
            summary_rows.append({"Tool": tool, "Region": region, "n": n,
                                  "mean_delta": np.nan, "median_delta": np.nan,
                                  "frac_region_higher": np.nan, "wilcoxon_p": np.nan})
            continue
        delta = sub[f"{tool}_{region}"] - sub[f"{tool}_whole"]
        try:
            _, p = wilcoxon(delta)
        except ValueError:
            p = np.nan
        summary_rows.append({
            "Tool": tool, "Region": region, "n": n,
            "mean_delta": round(float(delta.mean()), 4),
            "median_delta": round(float(delta.median()), 4),
            "frac_region_higher": round(float((delta > 0).mean()), 3),
            "wilcoxon_p": round(float(p), 4) if not np.isnan(p) else np.nan,
        })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(OUT / "topology_score_deltas.csv", index=False)
print("\nSaved output/topology_score_deltas.csv")
print(summary_df.to_string(index=False))

# ── Figure: mean delta (region - whole) per tool x region ──────────────────

REGION_COLOR = {"Cytoplasmic": "#0072B2", "Transmembrane": "#888888", "Extracellular/Lumenal": "#e74c3c"}

fig, ax = plt.subplots(figsize=(12, 5.5))
x = np.arange(len(TOOLS))
bar_w = 0.25
for i, region in enumerate(REGIONS):
    sub = summary_df[summary_df["Region"] == region].set_index("Tool").reindex(TOOLS)
    offset = (i - 1) * bar_w
    ax.bar(x + offset, sub["mean_delta"], bar_w, color=REGION_COLOR[region], alpha=0.85, label=region)
    for j, (tool, r) in enumerate(zip(TOOLS, sub.itertuples())):
        if pd.notna(r.wilcoxon_p) and r.wilcoxon_p < 0.05:
            ax.text(x[j] + offset, r.mean_delta + (0.01 if r.mean_delta >= 0 else -0.01),
                    "*", ha="center", va="bottom" if r.mean_delta >= 0 else "top", fontsize=11, color=REGION_COLOR[region])
ax.axhline(0, color="#333", linewidth=0.8)
ax.set_xticks(x)
ax.set_xticklabels(TOOLS, rotation=20, ha="right")
ax.set_ylabel("Mean Δ (region mean score − whole-protein score)")
ax.set_title("Does the region score higher than the whole protein?  (* = Wilcoxon p<0.05)")
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(FIG_DIR / "topology_score_deltas.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved output/figures/topology_score_deltas.png")
