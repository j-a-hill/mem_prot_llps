"""
rank_consensus.py — Rank-based consensus analysis of LLPS predictors on membrane proteins.

For each of 60 curated membrane LLPS proteins, predictor scores are converted to
normalised ranks (0–1, higher = more LLPS-prone), then summarised as an
AUROC-weighted mean rank and SD.  Outputs:
  - output/rank_consensus_table.csv   (per-protein summary)
  - output/figures/rank_consensus_spearman.png
  - output/figures/rank_consensus_heatmap.png
  - output/figures/rank_consensus_scatter.png
  - output/figures/rank_consensus_comparison.png
  - output/figures/rank_consensus_features.png
  - output/rank_consensus_report.md
"""

# %% ── Imports ────────────────────────────────────────────────────────────────
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

# %% ── Config ─────────────────────────────────────────────────────────────────
ROOT    = Path(__file__).parent
OUT_DIR = ROOT / "output"
FIG_DIR = OUT_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")
plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})

# 18 predictor columns and their human-readable names
PREDICTOR_COLS = {
    "PICNIC_score":      "PICNIC",
    "PICNIC_GO_score":   "PICNIC (GO)",
    "PSAP_score":        "PSAP",
    "PSPHunter_prob":    "PSPHunter",
    "PSPire_score":      "PSPire",
    "FuzDrop_pLLPS":     "FuzDrop",
    "catGRANULE_score":  "catGRANULE",
    "PLAAC_NLLR":        "PLAAC",
    "PScore_score":      "PScore",
    "ESpritz_score":     "ESpritz",
    "SEG_score":         "SEG",
    "SaPS_score":        "SaPS",
    "PdPS_score":        "PdPS",
    "DeepPhase_score":   "DeepPhase",
    "PDL_score":         "PDL",
    "RY_score":          "R+Y",
    "ParSe2_score":      "ParSe2",
    "LLPhyScore_score":  "LLPhyScore",
}
# Lower raw score = more LLPS-prone for this predictor; flip sign before ranking
INVERTED_PREDICTORS = {"LLPhyScore_score"}

# AUROC predictor name → background_scored column name
AUROC_TO_COL = {
    "PICNIC":       "PICNIC_score",
    "PICNIC (GO)":  "PICNIC_GO_score",
    "PSPire":       "PSPire_score",
    "PSPHunter":    "PSPHunter_prob",
    "SaPS":         "SaPS_score",
    "PdPS":         "PdPS_score",
    "PDL":          "PDL_score",
    "LLPhyScore":   "LLPhyScore_score",
    "PScore":       "PScore_score",
    "R+Y":          "RY_score",
    "ParSe2":       "ParSe2_score",
    "FuzDrop":      "FuzDrop_pLLPS",
    "PSAP":         "PSAP_score",
    "DeepPhase":    "DeepPhase_score",
    "catGRANULE":   "catGRANULE_score",
    "PLAAC":        "PLAAC_NLLR",
    "ESpritz":      "ESpritz_score",
    "SEG":          "SEG_score",
}

# Per-predictor leakage columns in leakage_map.csv
LEAKAGE_COL_MAP = {
    "PICNIC_score":     "pos_PICNIC",
    "PSAP_score":       "pos_PSAP",
    "PSPHunter_prob":   "pos_PSPHunter",
    "PSPire_score":     "pos_PSPire",
    "PDL_score":        "pos_PDL",
    "SaPS_score":       "pos_SaPS",
    "PdPS_score":       "pos_PdPS",
    "LLPhyScore_score": "pos_LLPhyScore",
    "DeepPhase_score":  "pos_DeePhase",
}


# %% ── Helper: Spearman CI via Fisher z-transform ─────────────────────────────
def _spearman_ci(r, n, alpha=0.05):
    """Asymptotic 95 % CI for Spearman ρ via Fisher z-transform."""
    z    = np.arctanh(np.clip(r, -0.9999, 0.9999))
    se   = 1.0 / np.sqrt(max(n - 3, 1))
    crit = stats.norm.ppf(1 - alpha / 2)
    return np.tanh(z - crit * se), np.tanh(z + crit * se)


# %% ── Step 1: Load data ──────────────────────────────────────────────────────
print("[step 1] Loading data ...")

# 1a. Metadata for the 60 membrane LLPS proteins
meta = pd.read_csv(OUT_DIR / "predictor_comparison.csv")
mem_ids = set(meta["Entry"])
print(f"  predictor_comparison: {len(meta)} proteins")

# 1b. Protein names from full dataset (best-effort; not critical)
try:
    full_df     = pd.read_csv(OUT_DIR / "full_dataset.csv", usecols=["Entry", "Protein names"])
    name_lookup = dict(zip(full_df["Entry"], full_df["Protein names"]))
except Exception:
    name_lookup = {}

# 1c. Score matrix — filter to the 60 membrane proteins
bg_all = pd.read_csv(OUT_DIR / "background_scored.csv")
scores = bg_all[bg_all["UniProt_ID"].isin(mem_ids)].copy()
scores = scores.set_index("UniProt_ID")
print(f"  background_scored filtered: {len(scores)} proteins")

# Flip LLPhyScore: lower raw score = more LLPS-prone (AUROC 0.376 < 0.5 inverted)
scores["LLPhyScore_score"] = -scores["LLPhyScore_score"]

# 1d. AUROC weights (Membrane vs membrane background scenario)
auroc_df = pd.read_csv(OUT_DIR / "clean_auroc_table.csv")
auroc_mem = auroc_df[auroc_df["Scenario"] == "Membrane vs membrane background"].copy()

# Fall back to AUROC_full when AUROC_clean is NaN
auroc_mem["weight"] = auroc_mem["AUROC_clean"].where(
    auroc_mem["AUROC_clean"].notna(), auroc_mem["AUROC_full"]
)

# Build col_name → AUROC weight map
weight_map = {}
for _, row in auroc_mem.iterrows():
    col = AUROC_TO_COL.get(row["Predictor"])
    if col:
        weight_map[col] = row["weight"]
print(f"  AUROC weights loaded for {len(weight_map)} predictors")

# 1e. Training-set leakage
leakage = pd.read_csv(OUT_DIR / "leakage_map.csv")
leakage = leakage.set_index("UniProt_ID")

# Normalise LLPhyScore column (has string 'unknown'/'True')
for c in leakage.columns:
    leakage[c] = leakage[c].map(lambda v: True if v is True or str(v) == "True" else False)

pos_cols  = [c for c in leakage.columns if c.startswith("pos_")]
leakage_any = leakage[pos_cols].any(axis=1).rename("leakage_any")
print(f"  Leakage: {leakage_any.sum()} proteins flagged in any predictor's positive training set")


# %% ── Step 2: Normalised ranks ───────────────────────────────────────────────
print("[step 2] Computing normalised ranks ...")

# Compute per-column ranks across all 60 proteins (rank over non-NaN values only)
rank_df = pd.DataFrame(index=scores.index, columns=list(PREDICTOR_COLS.values()), dtype=float)

for col, name in PREDICTOR_COLS.items():
    col_scores = scores[col].copy()
    valid_mask = col_scores.notna()
    valid_vals = col_scores[valid_mask].values
    n = len(valid_vals)
    if n == 0:
        continue
    # Higher rank = higher score = more LLPS-prone
    raw_ranks  = stats.rankdata(valid_vals, method="average")
    norm_ranks = 1.0 - (raw_ranks - 1) / (n - 1) if n > 1 else np.ones(n) * 0.5
    rank_df.loc[valid_mask, name] = norm_ranks

print(f"  Rank matrix: {rank_df.shape[0]} proteins × {rank_df.shape[1]} predictors")
nan_counts = rank_df.isna().sum()
print(f"  NaN per predictor: {nan_counts[nan_counts > 0].to_dict()}")


# %% ── Step 3: Per-protein summary ───────────────────────────────────────────
print("[step 3] Computing per-protein summary statistics ...")

summary_rows = []
for uid in rank_df.index:
    row = rank_df.loc[uid]
    valid = row.dropna()
    n_preds = len(valid)

    # Unweighted mean
    mean_rank = valid.mean() if n_preds > 0 else np.nan

    # AUROC-weighted mean — weights only for non-NaN predictors
    wt_vals, wt_weights = [], []
    for col, name in PREDICTOR_COLS.items():
        r = row[name]
        if not pd.isna(r) and col in weight_map:
            wt_vals.append(r)
            wt_weights.append(weight_map[col])
    if wt_weights:
        wt_arr = np.array(wt_weights)
        wt_arr = wt_arr / wt_arr.sum()
        weighted_mean_rank = float(np.dot(wt_vals, wt_arr))
    else:
        weighted_mean_rank = mean_rank  # fallback

    rank_sd = float(valid.std(ddof=1)) if n_preds >= 2 else np.nan

    summary_rows.append({
        "UniProt_ID":          uid,
        "mean_rank":           mean_rank,
        "weighted_mean_rank":  weighted_mean_rank,
        "rank_sd":             rank_sd,
        "n_predictors":        n_preds,
    })

summary = pd.DataFrame(summary_rows).set_index("UniProt_ID")

# Join metadata
meta_idx = meta.set_index("Entry")
summary["Entry_name"]    = meta_idx["Entry name"]
summary["p(LLPS)"]       = meta_idx["p(LLPS)"]
summary["pLLPS_Class"]   = meta_idx["pLLPS_Class"]
summary["TMD_count"]     = meta_idx["TMD_count"]
summary["Length"]        = meta_idx["Length"]
summary["Protein_Name"]  = summary.index.map(lambda x: name_lookup.get(x, ""))
summary["leakage_any"]   = leakage_any.reindex(summary.index).fillna(False)

summary = summary.sort_values("weighted_mean_rank", ascending=False)

# Save
out_csv = OUT_DIR / "rank_consensus_table.csv"
summary.reset_index().to_csv(out_csv, index=False)
print(f"  Saved {out_csv}  ({len(summary)} rows)")


# %% ── Step 4: Pairwise Spearman correlation heatmap ──────────────────────────
print("[step 4] Pairwise Spearman correlation clustermap ...")

corr_mat = rank_df.corr(method="spearman")

# Clustered heatmap
mask_tri = np.zeros_like(corr_mat.values, dtype=bool)  # show all cells

g = sns.clustermap(
    corr_mat,
    method="average",
    metric="euclidean",
    cmap="RdBu_r",
    vmin=-1, vmax=1,
    annot=True,
    fmt=".1f",
    annot_kws={"size": 7},
    figsize=(13, 12),
    linewidths=0.4,
    cbar_kws={"shrink": 0.6},
)
g.fig.suptitle(
    "Pairwise Spearman correlations between LLPS predictor normalised ranks",
    y=1.01, fontsize=13,
)
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_ylabel("")
g.fig.text(0.5, -0.01, "n = 60 membrane LLPS proteins", ha="center", fontsize=10, style="italic")

g.fig.savefig(FIG_DIR / "rank_consensus_spearman.png", dpi=150, bbox_inches="tight")
plt.close("all")
print(f"  Saved rank_consensus_spearman.png")

# Extract clustered column order for Step 5
clustered_col_order = rank_df.columns[g.dendrogram_col.reordered_ind].tolist()


# %% ── Step 5: Protein × predictor heatmap ────────────────────────────────────
print("[step 5] Protein × predictor heatmap ...")

# Row order: sorted by weighted_mean_rank (highest at top = first row)
row_order = summary.sort_values("weighted_mean_rank", ascending=True).index.tolist()
# (ascending=True so last protein in list is plotted at top in heatmap)

heat_data = rank_df.loc[row_order, clustered_col_order]

# TMD category for sidebar
def _tmd_cat(n):
    if n == 1:
        return "single-pass"
    elif 2 <= n <= 6:
        return "oligopass"
    else:
        return "multipass"

tmd_cats = summary.loc[row_order, "TMD_count"].apply(_tmd_cat)
tmd_colors = {"single-pass": "#0072B2", "oligopass": "#E69F00", "multipass": "#009E73"}
tmd_rgb = tmd_cats.map(tmd_colors)

# Build per-predictor leakage matrix (protein × predictor, human-readable names)
leak_mat = pd.DataFrame(False, index=row_order, columns=clustered_col_order)
for col, pos_col in LEAKAGE_COL_MAP.items():
    name = PREDICTOR_COLS.get(col)
    if name and name in leak_mat.columns and pos_col in leakage.columns:
        for uid in row_order:
            if uid in leakage.index:
                leak_mat.at[uid, name] = bool(leakage.at[uid, pos_col])

# Figure layout: TMD sidebar | main heatmap | mean_rank bar
n_rows = len(row_order)
n_cols = len(clustered_col_order)
fig_h  = max(14, n_rows * 0.30 + 3)
fig_w  = max(14, n_cols * 0.72 + 4)

fig = plt.figure(figsize=(fig_w, fig_h))
gs  = fig.add_gridspec(1, 3, width_ratios=[0.4, n_cols, 1.8], wspace=0.02)

ax_tmd  = fig.add_subplot(gs[0])
ax_heat = fig.add_subplot(gs[1])
ax_bar  = fig.add_subplot(gs[2])

# ── TMD sidebar ───────────────────────────────────────────────────────────────
for i, uid in enumerate(row_order):
    color = tmd_rgb[uid]
    ax_tmd.barh(i, 1, color=color, edgecolor="none")
ax_tmd.set_xlim(0, 1)
ax_tmd.set_ylim(-0.5, n_rows - 0.5)
ax_tmd.set_yticks(range(n_rows))
ax_tmd.set_yticklabels(summary.loc[row_order, "Entry_name"], fontsize=7)
ax_tmd.set_xticks([])
ax_tmd.set_xlabel("TMD", fontsize=9)
ax_tmd.spines[["top", "right", "bottom"]].set_visible(False)

legend_patches = [
    mpatches.Patch(color=tmd_colors["single-pass"], label="Single-pass (1)"),
    mpatches.Patch(color=tmd_colors["oligopass"],   label="Oligopass (2–6)"),
    mpatches.Patch(color=tmd_colors["multipass"],   label="Multipass (7+)"),
]
ax_tmd.legend(handles=legend_patches, loc="lower left", fontsize=7, framealpha=0.8,
              bbox_to_anchor=(0, -0.06))

# ── Main heatmap ──────────────────────────────────────────────────────────────
sns.heatmap(
    heat_data,
    ax=ax_heat,
    cmap="YlOrRd_r",
    vmin=0, vmax=1,
    linewidths=0.2,
    linecolor="#dddddd",
    annot=False,
    cbar=True,
    cbar_kws={"shrink": 0.4, "label": "Normalised rank\n(1 = highest score)"},
    yticklabels=False,
    xticklabels=True,
)

# Annotate leakage cells with '*'
for row_i, uid in enumerate(row_order):
    for col_i, pred_name in enumerate(clustered_col_order):
        if leak_mat.at[uid, pred_name]:
            ax_heat.text(
                col_i + 0.5, row_i + 0.5, "*",
                ha="center", va="center",
                fontsize=7, color="#222222", fontweight="bold",
            )

ax_heat.set_xlabel("")
ax_heat.set_ylabel("")
ax_heat.tick_params(axis="x", rotation=45, labelsize=8)
ax_heat.set_title(
    "Normalised LLPS predictor ranks: 60 membrane proteins\n"
    "(* = in positive training set for that predictor)",
    fontsize=11, pad=8,
)

# ── Weighted mean rank bar ────────────────────────────────────────────────────
wmr = summary.loc[row_order, "weighted_mean_rank"].values
bar_colors = plt.cm.YlOrRd_r(wmr)  # type: ignore[call-arg]

ax_bar.barh(range(n_rows), wmr, color=bar_colors[:, :3], edgecolor="none")
ax_bar.set_xlim(0, 1)
ax_bar.set_ylim(-0.5, n_rows - 0.5)
ax_bar.set_yticks([])
ax_bar.set_xlabel("Weighted\nmean rank", fontsize=9)
ax_bar.axvline(0.5, color="#555555", linewidth=0.8, linestyle="--", alpha=0.6)
ax_bar.spines[["top", "right"]].set_visible(False)
ax_bar.tick_params(axis="x", labelsize=8)

plt.savefig(FIG_DIR / "rank_consensus_heatmap.png", dpi=150, bbox_inches="tight")
plt.close("all")
print(f"  Saved rank_consensus_heatmap.png")


# %% ── Step 6: Mean vs SD scatter ────────────────────────────────────────────
print("[step 6] Mean-rank vs SD scatter ...")

scatter_df = summary.dropna(subset=["weighted_mean_rank", "rank_sd"]).copy()
scatter_df["TMD_class"] = scatter_df["TMD_count"].apply(
    lambda n: "Single-pass (1 TMD)" if n == 1 else "Multi-pass (≥2 TMDs)"
)

c_single = "#0072B2"
c_multi  = "#E69F00"
color_map_sc = {"Single-pass (1 TMD)": c_single, "Multi-pass (≥2 TMDs)": c_multi}
colors = scatter_df["TMD_class"].map(color_map_sc)

# Point size proportional to p(LLPS)
sizes = 40 + 180 * scatter_df["p(LLPS)"].fillna(0.5)

fig, ax = plt.subplots(figsize=(11, 8))

for cls, grp in scatter_df.groupby("TMD_class"):
    idx = grp.index
    ax.scatter(
        grp["weighted_mean_rank"], grp["rank_sd"],
        s=sizes[idx], color=color_map_sc[cls], alpha=0.80,
        edgecolors="white", linewidths=0.5, label=cls, zorder=3,
    )

# Reference lines and shaded regions
xmean = scatter_df["weighted_mean_rank"].mean()
xsd   = scatter_df["weighted_mean_rank"].std()
ax.axvline(xmean - xsd, color="#999999", linewidth=0.9, linestyle="--", alpha=0.7)
ax.axvline(xmean + xsd, color="#999999", linewidth=0.9, linestyle="--", alpha=0.7)

ax.axvspan(0.7, 1.0,  alpha=0.07, color="#2ecc71", zorder=1, label="Consistently high (>0.7)")
ax.axvspan(0.0, 0.3,  alpha=0.07, color="#e74c3c", zorder=1, label="Consistently low (<0.3)")

# Label: top-5 highest mean, bottom-5 lowest mean, top-5 highest SD
top5_mean = scatter_df.nlargest(5,  "weighted_mean_rank").index
bot5_mean = scatter_df.nsmallest(5, "weighted_mean_rank").index
top5_sd   = scatter_df.nlargest(5,  "rank_sd").index
to_label  = set(top5_mean) | set(bot5_mean) | set(top5_sd)

for uid in to_label:
    if uid not in scatter_df.index:
        continue
    x  = scatter_df.at[uid, "weighted_mean_rank"]
    y  = scatter_df.at[uid, "rank_sd"]
    lbl = scatter_df.at[uid, "Entry_name"]
    ax.annotate(
        lbl, (x, y),
        xytext=(5, 3), textcoords="offset points",
        fontsize=7.5, color="#333333",
        arrowprops={"arrowstyle": "-", "color": "#aaaaaa", "lw": 0.6},
    )

# Size legend
for pllps, sz in [(0.4, 40 + 180 * 0.4), (0.7, 40 + 180 * 0.7), (1.0, 40 + 180 * 1.0)]:
    ax.scatter([], [], s=sz, color="#888888", alpha=0.8, edgecolors="white",
               label=f"p(LLPS) = {pllps:.1f}")

ax.set_xlabel("Weighted mean normalised rank (↑ = consistently predicted)", fontsize=11)
ax.set_ylabel("SD of normalised ranks (↑ = more predictor disagreement)", fontsize=11)
ax.set_title("Predictor consensus vs disagreement for 60 membrane LLPS proteins", fontsize=12)
ax.legend(fontsize=8, loc="upper left", framealpha=0.85)
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(bottom=0)

fig.text(
    0.5, -0.02,
    "LLPhyScore sign-flipped before ranking.  "
    "Weights = AUROC (membrane vs membrane background).  "
    "Point size ∝ p(LLPS).",
    ha="center", fontsize=8, style="italic", color="#555555",
)

plt.tight_layout()
plt.savefig(FIG_DIR / "rank_consensus_scatter.png", dpi=150, bbox_inches="tight")
plt.close("all")
print(f"  Saved rank_consensus_scatter.png")


# %% ── Step 6b: Unweighted vs AUROC-weighted comparison ──────────────────────
print("[step 6b] Unweighted vs weighted rank comparison ...")

cmp_df = summary.dropna(subset=["mean_rank", "weighted_mean_rank", "rank_sd"]).copy()
cmp_df["TMD_class"] = cmp_df["TMD_count"].apply(
    lambda n: "Single-pass (1 TMD)" if n == 1 else "Multi-pass (≥2 TMDs)"
)
cmp_colors = cmp_df["TMD_class"].map({"Single-pass (1 TMD)": c_single, "Multi-pass (≥2 TMDs)": c_multi})

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# ── Panel A: unweighted mean rank vs SD ──────────────────────────────────────
ax = axes[0]
for cls, grp in cmp_df.groupby("TMD_class"):
    ax.scatter(
        grp["mean_rank"], grp["rank_sd"],
        s=60, color=color_map_sc[cls], alpha=0.80,
        edgecolors="white", linewidths=0.5, label=cls, zorder=3,
    )
ax.axvspan(0.7, 1.0, alpha=0.07, color="#2ecc71", zorder=1)
ax.axvspan(0.0, 0.3, alpha=0.07, color="#e74c3c", zorder=1)
ax.set_xlabel("Unweighted mean normalised rank", fontsize=10)
ax.set_ylabel("SD of normalised ranks", fontsize=10)
ax.set_title("A  Unweighted", fontsize=11, fontweight="bold")
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(bottom=0)
ax.legend(fontsize=8, loc="upper left", framealpha=0.85)
ax.spines[["top", "right"]].set_visible(False)

# ── Panel B: AUROC-weighted mean rank vs SD ───────────────────────────────────
ax = axes[1]
for cls, grp in cmp_df.groupby("TMD_class"):
    ax.scatter(
        grp["weighted_mean_rank"], grp["rank_sd"],
        s=60, color=color_map_sc[cls], alpha=0.80,
        edgecolors="white", linewidths=0.5, label=cls, zorder=3,
    )
ax.axvspan(0.7, 1.0, alpha=0.07, color="#2ecc71", zorder=1)
ax.axvspan(0.0, 0.3, alpha=0.07, color="#e74c3c", zorder=1)
ax.set_xlabel("AUROC-weighted mean normalised rank", fontsize=10)
ax.set_ylabel("SD of normalised ranks", fontsize=10)
ax.set_title("B  AUROC-weighted", fontsize=11, fontweight="bold")
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(bottom=0)
ax.spines[["top", "right"]].set_visible(False)

# ── Panel C: direct unweighted vs weighted comparison ─────────────────────────
ax = axes[2]
for cls, grp in cmp_df.groupby("TMD_class"):
    ax.scatter(
        grp["mean_rank"], grp["weighted_mean_rank"],
        s=60, color=color_map_sc[cls], alpha=0.80,
        edgecolors="white", linewidths=0.5, label=cls, zorder=3,
    )

# Identity line
ax.plot([0, 1], [0, 1], color="#999999", linewidth=1.0, linestyle="--", zorder=2)

# Label proteins where weighting shifts rank by > 0.1
cmp_df["rank_shift"] = cmp_df["weighted_mean_rank"] - cmp_df["mean_rank"]
shifted = cmp_df[cmp_df["rank_shift"].abs() > 0.1]
for uid, row in shifted.iterrows():
    ax.annotate(
        row["Entry_name"],
        (row["mean_rank"], row["weighted_mean_rank"]),
        xytext=(5, 3), textcoords="offset points",
        fontsize=7, color="#333333",
        arrowprops={"arrowstyle": "-", "color": "#aaaaaa", "lw": 0.6},
    )

ax.set_xlabel("Unweighted mean normalised rank", fontsize=10)
ax.set_ylabel("AUROC-weighted mean normalised rank", fontsize=10)
ax.set_title("C  Unweighted vs weighted", fontsize=11, fontweight="bold")
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
ax.spines[["top", "right"]].set_visible(False)

fig.suptitle(
    "Unweighted vs AUROC-weighted normalised ranks: 60 membrane LLPS proteins",
    fontsize=12,
)
fig.text(
    0.5, -0.03,
    "Weights = AUROC (membrane vs membrane background).  "
    "Panel C: dashed line = identity (no change from weighting); "
    "labelled points shifted >0.1 rank units.",
    ha="center", fontsize=8, style="italic", color="#555555",
)
plt.tight_layout()
plt.savefig(FIG_DIR / "rank_consensus_comparison.png", dpi=150, bbox_inches="tight")
plt.close("all")
print(f"  Saved rank_consensus_comparison.png")


# %% ── Step 7: Feature correlations ──────────────────────────────────────────
print("[step 7] Feature correlations ...")

fc_df = summary.copy()
fc_df["log10_Length"] = np.log10(fc_df["Length"].clip(lower=1))

features = {
    "log10_Length": "log₁₀(Length)",
    "TMD_count":    "TMD count",
}
targets = {
    "weighted_mean_rank": "Weighted mean rank",
    "rank_sd":            "Rank SD",
}

results = []
n_proteins = len(fc_df)
for feat_col, feat_label in features.items():
    for tgt_col, tgt_label in targets.items():
        valid = fc_df[[feat_col, tgt_col]].dropna()
        if len(valid) < 5:
            results.append({"feature": feat_label, "target": tgt_label,
                             "rho": np.nan, "p": np.nan})
            continue
        rho, p = stats.spearmanr(valid[feat_col], valid[tgt_col])
        results.append({"feature": feat_label, "target": tgt_label,
                         "rho": rho, "p": p, "n": len(valid)})

res_df = pd.DataFrame(results)

# BH FDR correction across all (feature × target) combinations
valid_mask = res_df["p"].notna()
p_vals = res_df.loc[valid_mask, "p"].values
_, p_adj, _, _ = multipletests(p_vals, method="fdr_bh")
res_df.loc[valid_mask, "p_adj"] = p_adj
res_df["significant"] = res_df["p_adj"] < 0.05

# CI via Fisher z-transform
n_obs = n_proteins  # approximate; exact n varies slightly per pair
res_df["ci_lo"] = res_df["rho"].apply(lambda r: _spearman_ci(r, n_obs)[0] if pd.notna(r) else np.nan)
res_df["ci_hi"] = res_df["rho"].apply(lambda r: _spearman_ci(r, n_obs)[1] if pd.notna(r) else np.nan)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

feat_labels = list(features.values())
tgt_cols    = list(targets.keys())
tgt_labels  = list(targets.values())

for ax_i, (tgt_col, tgt_label) in enumerate(zip(tgt_cols, tgt_labels)):
    ax   = axes[ax_i]
    sub  = res_df[res_df["target"] == tgt_label].copy()
    sub  = sub.set_index("feature").reindex(feat_labels).reset_index()

    rho_vals = sub["rho"].values
    ci_lo    = sub["ci_lo"].values
    ci_hi    = sub["ci_hi"].values
    err_lo   = rho_vals - ci_lo
    err_hi   = ci_hi - rho_vals
    sigs     = sub["significant"].fillna(False).values

    bar_cols = ["#2c7bb6" if s else "#a6cee3" for s in sigs]
    y_pos    = range(len(feat_labels))

    ax.barh(
        y_pos, rho_vals,
        xerr=[err_lo, err_hi],
        color=bar_cols, edgecolor="none", error_kw={"ecolor": "#555555", "capsize": 3},
    )
    ax.axvline(0, color="#333333", linewidth=0.9)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(feat_labels, fontsize=10)
    ax.set_xlabel("Spearman ρ", fontsize=11)
    ax.set_title(tgt_label, fontsize=11)
    ax.set_xlim(-1, 1)
    ax.spines[["top", "right"]].set_visible(False)

    # Significance annotations
    for y_i, (rho, sig) in enumerate(zip(rho_vals, sigs)):
        if sig and not np.isnan(rho):
            xpos = rho + (0.04 if rho >= 0 else -0.04)
            ax.text(xpos, y_i, "*", ha="center", va="center", fontsize=12, color="#e74c3c")

# Legend for significance
sig_patch   = mpatches.Patch(color="#2c7bb6",  label="BH-adjusted p < 0.05")
insig_patch = mpatches.Patch(color="#a6cee3", label="BH-adjusted p ≥ 0.05")
axes[1].legend(handles=[sig_patch, insig_patch], fontsize=8, loc="lower right")

fig.suptitle(
    "Feature correlations with consensus rank and predictor disagreement",
    fontsize=12,
)
plt.tight_layout()
plt.savefig(FIG_DIR / "rank_consensus_features.png", dpi=150, bbox_inches="tight")
plt.close("all")
print(f"  Saved rank_consensus_features.png")


# %% ── Step 8: Markdown report ────────────────────────────────────────────────
print("[step 8] Writing markdown report ...")

# Helper to build tables
def _md_table(df_sub, cols):
    header = "| " + " | ".join(cols) + " |"
    sep    = "| " + " | ".join(["---"] * len(cols)) + " |"
    rows   = []
    for _, r in df_sub.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                cells.append(f"{v:.3f}")
            else:
                cells.append(str(v))
        rows.append("| " + " | ".join(cells) + " |")
    return "\n".join([header, sep] + rows)

# Top-10 / bottom-10 by weighted_mean_rank
sum_reset = summary.reset_index()
display_cols = ["Entry_name", "weighted_mean_rank", "rank_sd", "leakage_any", "TMD_count"]

top10  = sum_reset.head(10)[display_cols]
bot10  = sum_reset.tail(10).sort_values("weighted_mean_rank")[display_cols]
disc10 = sum_reset.nlargest(10, "rank_sd").sort_values("rank_sd", ascending=False)[display_cols]

# Feature correlation prose
sig_results = res_df[res_df["significant"] == True].copy()

sig_prose_lines = []
for _, row in sig_results.iterrows():
    direction = "positively" if row["rho"] > 0 else "negatively"
    sig_prose_lines.append(
        f"  - **{row['feature']}** is {direction} correlated with **{row['target']}** "
        f"(ρ = {row['rho']:.2f}, BH-adjusted p = {row['p_adj']:.3f})."
    )
sig_prose = "\n".join(sig_prose_lines) if sig_prose_lines else (
    "  No features showed statistically significant correlations after BH correction "
    "(n = 60 limits power)."
)

# AUROC weight table
auroc_tbl_rows = []
for col, name in PREDICTOR_COLS.items():
    w = weight_map.get(col, np.nan)
    inv = " (sign-flipped)" if col in INVERTED_PREDICTORS else ""
    auroc_tbl_rows.append(f"| {name}{inv} | {w:.3f} |" if not np.isnan(w) else f"| {name}{inv} | N/A |")
auroc_tbl = "| Predictor | AUROC weight |\n| --- | --- |\n" + "\n".join(auroc_tbl_rows)

report = f"""# Rank-Based Consensus Analysis

## Methods

### Rationale
LLPS predictors output scores on heterogeneous scales (probabilities, log-likelihood
ratios, raw neural-network outputs, etc.) that cannot be meaningfully averaged
directly.  Converting each predictor's scores to normalised ranks removes
scale dependence and makes the consensus robust to outlier scores and
different score distributions.

### Score directionality
All 18 predictors are listed below.  For 17 of them, a higher raw score indicates
greater LLPS propensity.  The exception is **LLPhyScore**, for which a *lower* raw
score is more LLPS-prone (AUROC 0.376 < 0.5 when left as-is).  Its raw score is
sign-flipped before ranking so that higher rank still means more LLPS-prone.

{auroc_tbl}

### Normalisation formula
For each predictor, ranks are computed over the N valid (non-NaN) scores using
average-rank ties (scipy.stats.rankdata).  Normalised rank = 1 − (rank − 1) / (N − 1),
mapping the highest-scoring protein to 1.0 and the lowest to 0.0.  Proteins
without a score for a given predictor receive NaN for that predictor's normalised
rank and are excluded from that predictor's ranking.

### AUROC weighting
Weights for the weighted mean rank are taken from the **"Membrane vs membrane
background"** scenario in `output/clean_auroc_table.csv`, using `AUROC_clean`
(leakage-aware) when available, falling back to `AUROC_full` when `AUROC_clean`
is NaN (this affects LLPhyScore only).  Weights are normalised to sum to 1 over
non-NaN predictors per protein before taking the dot product.

### Training set leakage
`output/leakage_map.csv` records whether each protein appeared in the positive
training set of each predictor.  A protein is marked `leakage_any = True` if it
was in the positive training set of *any* predictor.  In the heatmap, individual
cells are asterisked where a protein was in the positive training set specifically
for that predictor.

---

## Key findings

### Top 10 consistently predicted proteins (highest weighted mean rank)

{_md_table(top10, display_cols)}

### Bottom 10 systematically underpredicted proteins (lowest weighted mean rank)

{_md_table(bot10, display_cols)}

### Most discordant proteins (highest rank SD)

{_md_table(disc10, display_cols)}

### Feature correlations
Spearman ρ was computed between each feature and (a) `weighted_mean_rank` and
(b) `rank_sd`, with Benjamini–Hochberg FDR correction across all {len(features) * len(targets)}
(feature × target) combinations.

{sig_prose}

See `output/figures/rank_consensus_features.png` for the full barplot with
95 % confidence intervals.

---

## Caveats

- **LLPhyScore** performs poorly on membrane proteins (AUROC 0.376 before sign
  flip, 0.624 after); it is included but interpret its contribution with caution.
- **Training-set leakage** may artificially inflate ranks for flagged proteins;
  use `AUROC_clean` weights (which use leakage-free test sets) to mitigate this.
- **n = 60** limits statistical power for feature correlations; treat non-significant
  results as inconclusive rather than evidence of no effect.
- The consensus reflects **inter-predictor agreement**, not ground truth.  High
  agreement among predictors that share training data may reflect correlated
  biases rather than true LLPS propensity.

---

## Output files

| File | Description |
| --- | --- |
| `output/rank_consensus_table.csv` | Per-protein summary (mean rank, weighted mean rank, SD, leakage flag, metadata) |
| `output/figures/rank_consensus_spearman.png` | Clustered Spearman correlation heatmap between predictor normalised ranks |
| `output/figures/rank_consensus_heatmap.png` | Protein × predictor normalised rank heatmap |
| `output/figures/rank_consensus_scatter.png` | AUROC-weighted mean rank vs rank SD scatter plot |
| `output/figures/rank_consensus_comparison.png` | Side-by-side unweighted vs weighted scatters; direct comparison panel |
| `output/figures/rank_consensus_features.png` | Feature correlations with consensus rank and disagreement |
| `output/rank_consensus_report.md` | This report |
"""

report_path = OUT_DIR / "rank_consensus_report.md"
report_path.write_text(report)
print(f"  Saved {report_path}")

print("\n[done] All outputs written.")
print(f"  CSV:     {out_csv}")
print(f"  Figures: {FIG_DIR}/rank_consensus_*.png")
print(f"  Report:  {report_path}")
