"""
plot_roc.py

ROC/precision-recall curves, MCC analysis, training-set contamination masking,
and matched-background sensitivity analysis for LLPS predictors.

Scenarios
---------
  All experimental vs proteome        All experimental LLPS (882) vs proteome background
  Membrane vs proteome background     Membrane experimental (60)  vs proteome background (LLPS-free negatives)
  Membrane vs membrane background     Membrane experimental (60)  vs membrane background (TMD>0)
  Membrane vs single-pass background  Membrane experimental (60)  vs single-pass membrane background
  Membrane vs multi-pass background   Membrane experimental (60)  vs multi-pass membrane background

Outputs — output/figures/v{N}/ (auto-incrementing — every run gets a new
version directory, e.g. v1/, v2/, ...; nothing is overwritten in place)
  roc_curves.png                  ROC grid, all predictors × all scenarios
  auroc_summary.png               AUROC bar chart with 95% bootstrap CI
  auprc_summary.png               AUPRC bar chart
  mcc_curves.png                  MCC-vs-threshold grid
  mcc_summary.png                 Max-MCC bar chart with CI
  auroc_clean_delta.png           Full vs per-tool-clean AUROC, membrane vs membrane background
  auroc_clean_all.png             Per-tool-clean AUROC + Max-MCC, all scenarios
  roc_curves_clean.png            ROC: full vs per-tool-clean, membrane vs membrane background
  roc_curves_clean_scenarios.png  ROC: per-tool-clean, all 4 background scenarios × all predictors
  roc_curves_clean_vs_global_*.png  ROC: per-tool-clean vs globally-clean, one file per scenario
  auroc_matched_proteome.png      Matched vs full — proteome background
  auroc_matched_membrane.png      Matched vs full — membrane background

Outputs — output/
  clean_auroc_table.csv       AUROC/MaxMCC results (full / clean / global)
  matched_auroc_table.csv     Matched-background mean AUROC ± SD (100 draws)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc, average_precision_score

ROOT     = Path(__file__).parent
FIG_BASE = ROOT / "output" / "figures"


def _next_version_dir(base):
    """Auto-incrementing output/figures/v{N}/ — every run gets a fresh version."""
    base.mkdir(parents=True, exist_ok=True)
    existing = [int(p.name[1:]) for p in base.iterdir()
                if p.is_dir() and p.name.startswith("v") and p.name[1:].isdigit()]
    next_v = max(existing, default=0) + 1
    d = base / f"v{next_v}"
    d.mkdir(parents=True, exist_ok=True)
    return d


FIG_DIR = _next_version_dir(FIG_BASE)
print(f"Writing figures to {FIG_DIR.relative_to(ROOT)}")

# ── Data ──────────────────────────────────────────────────────────────────────

df    = pd.read_csv(ROOT / "output" / "background_scored.csv")
masks = pd.read_csv(ROOT / "output" / "clean_masks.csv")
df    = df.merge(masks, on="UniProt_ID", how="left")

# Ordered: LLPS-specific predictors first (grouped by mechanism), then the
# IDR/disorder predictors that are general-purpose tools repurposed as LLPS
# proxies (PLAAC, ESpritz, SEG) — these were never trained on LLPS labels.
PREDICTORS = {
    "PICNIC":      "PICNIC_score",
    "PICNIC (GO)": "PICNIC_GO_score",
    "PSPire":      "PSPire_score",
    "PSPHunter":   "PSPHunter_prob",
    "SaPS":        "SaPS_score",
    "PdPS":        "PdPS_score",
    "PDL":         "PDL_score",
    "LLPhyScore":  "LLPhyScore_score",
    "PScore":      "PScore_score",
    "R+Y":         "RY_score",
    "ParSe2":      "ParSe2_score",
    "FuzDrop":     "FuzDrop_pLLPS",
    "PSAP":        "PSAP_score",
    "DeepPhase":   "DeepPhase_score",
    "catGRANULE":  "catGRANULE_score",
    "PLAAC":       "PLAAC_NLLR",
    "ESpritz":     "ESpritz_score",
    "SEG":         "SEG_score",
}

for c in PREDICTORS.values():
    df[c] = pd.to_numeric(df[c], errors="coerce")

# ── Predictor categorisation ──────────────────────────────────────────────────
# Target type: what the tool was actually designed/trained to predict.
IDR_PROXY_TOOLS = {"PLAAC", "ESpritz", "SEG"}
TARGET_TYPE = {p: ("IDR proxy" if p in IDR_PROXY_TOOLS else "LLPS") for p in PREDICTORS}

# Mechanism: the dominant feature/theory driving the method. Per
# predictor_classification.md — only PICNIC/PICNIC(GO)/PSPire actually use
# AlphaFold2 structure; PSPHunter uses annotations (not structure), and
# SaPS/PdPS/PDL/LLPhyScore are non-structural multi-feature ML ensembles.
# NOTE: SaPS and PdPS are the same PhaSePred "10-feature gradient boost"
# model trained on two different target labels (scaffold vs partner-
# dependent/client), and both take catGRANULE/PLAAC/PScore/ESpritz/SEG/
# DeepPhase scores as *input features* — they are not independent of those
# six tools also in this comparison. PDL is unrelated (separate paper/
# architecture, ProtT5+KmerConv) and only shares a distribution file with
# SaPS/PdPS for benchmarking convenience.
MECHANISM = {
    "PICNIC": "AF2 structure + sequence", "PICNIC (GO)": "AF2 structure + sequence",
    "PSPire": "AF2 structure + sequence",
    "PSPHunter": "Multi-feature ML ensemble", "SaPS": "Multi-feature ML ensemble",
    "PdPS": "Multi-feature ML ensemble", "PDL": "Multi-feature ML ensemble",
    "LLPhyScore": "Multi-feature ML ensemble",
    "PScore": "Pi-pi / cation-pi", "R+Y": "Pi-pi / cation-pi",
    "ParSe2": "Sticker-spacer / polymer theory", "FuzDrop": "Sticker-spacer / polymer theory",
    "PSAP": "IDR / disorder", "DeepPhase": "IDR / disorder", "catGRANULE": "IDR / disorder",
    "PLAAC": "IDR / disorder", "ESpritz": "IDR / disorder", "SEG": "IDR / disorder",
}

MECHANISM_COLOR = {
    "AF2 structure + sequence":         "#3366CC",
    "Multi-feature ML ensemble":        "#7755AA",
    "Pi-pi / cation-pi":                "#CC3333",
    "Sticker-spacer / polymer theory":  "#CC8800",
    "IDR / disorder":                   "#228855",
}

# Short forms for box labels — some mechanism blocks are only 2 predictors
# wide, too narrow for the full name without bleeding into the next box.
CATEGORY_LABEL_SHORT = {
    "AF2 structure + sequence":         "AF2 structure",
    "Multi-feature ML ensemble":        "ML ensemble",
    "Pi-pi / cation-pi":                "Pi-pi",
    "Sticker-spacer / polymer theory":  "Sticker-spacer",
    "IDR / disorder (LLPS-trained)":    "IDR/disorder (LLPS)",
    "IDR / disorder (IDR proxy)":       "IDR/disorder (proxy)",
}

# Category blocks shown as a labelled box on every figure — one per mechanism,
# split further by target type where a mechanism spans both (IDR / disorder
# covers both LLPS-trained tools and the repurposed IDR-proxy tools).
CATEGORY_BLOCKS = [
    (mech, [p for p in PREDICTORS if MECHANISM[p] == mech and TARGET_TYPE[p] == "LLPS"])
    for mech in dict.fromkeys(MECHANISM.values())
] + [("IDR / disorder (IDR proxy)",
      [p for p in PREDICTORS if MECHANISM[p] == "IDR / disorder" and TARGET_TYPE[p] == "IDR proxy"])]
# Drop the now-empty plain "IDR / disorder" entry's IDR-proxy members were
# pulled out above; relabel the LLPS-trained one for clarity.
CATEGORY_BLOCKS = [
    ("IDR / disorder (LLPS-trained)" if name == "IDR / disorder" else name, members)
    for name, members in CATEGORY_BLOCKS if members
]

# label_col, filter_col, colour, linestyle
SCENARIOS = {
    "All experimental vs proteome":       ("label_all",               None,                       "#888888", "-"),
    "Membrane vs proteome background":    ("label_mem_fixed",         "label_mem_fixed",          "#0072B2", "-"),
    "Membrane vs membrane background":    ("label_mem_bg",            "label_mem_bg",             "#e74c3c", "-"),
    "Membrane vs single-pass background": ("label_mem_singlepass_bg", "label_mem_singlepass_bg",  "#009E73", "--"),
    "Membrane vs multi-pass background":  ("label_mem_multipass_bg",  "label_mem_multipass_bg",   "#CC79A7", "--"),
}

# Compact tags for in-panel legends (full names already given by the
# figure-level scenario legend, so these just need to be short).
SCENARIO_TAG = {
    "All experimental vs proteome":       "All-exp",
    "Membrane vs proteome background":    "Proteome bg",
    "Membrane vs membrane background":    "Membrane bg",
    "Membrane vs single-pass background": "Single-pass bg",
    "Membrane vs multi-pass background":  "Multi-pass bg",
}

# Filename-safe slugs for the per-scenario clean-vs-global ROC figures (§3c).
SCENARIO_SLUG = {
    "Membrane vs proteome background":    "proteome_bg",
    "Membrane vs membrane background":    "membrane_bg",
    "Membrane vs single-pass background": "singlepass_bg",
    "Membrane vs multi-pass background":  "multipass_bg",
}

PRIMARY_SC  = "Membrane vs membrane background"
pred_labels = list(PREDICTORS.keys())
sc_names    = list(SCENARIOS.keys())

ncols = 4


def grid_layout(ncols):
    """(row, col) per predictor; each category block starts on a new row."""
    positions = {}
    row = 0
    block_rows = []  # (start_row, end_row) per block, in CATEGORY_BLOCKS order
    for _, members in CATEGORY_BLOCKS:
        start_row = row
        for j, p in enumerate(members):
            positions[p] = (row + j // ncols, j % ncols)
        row += -(-len(members) // ncols)  # ceil division
        block_rows.append((start_row, row - 1))
    return positions, row, block_rows


GRID_POS, nrows, BLOCK_ROWS = grid_layout(ncols)


def new_predictor_grid(figsize_cell, suptitle):
    """Figure/axes grid with unused cells hidden and the category row-break applied."""
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * figsize_cell[0], nrows * figsize_cell[1]))
    axes = np.atleast_2d(axes)
    used = set(GRID_POS.values())
    for r in range(nrows):
        for c in range(ncols):
            if (r, c) not in used:
                axes[r, c].set_visible(False)
    fig.suptitle(suptitle, fontsize=12, y=1.01)
    return fig, axes


def style_predictor_title(ax, pred_label, fontsize=10):
    """IDR-proxy tools flagged in italics; mechanism is shown by the category box instead."""
    is_proxy = TARGET_TYPE[pred_label] == "IDR proxy"
    title = pred_label + ("  (IDR proxy)" if is_proxy else "")
    ax.set_title(title, fontsize=fontsize, fontstyle="italic" if is_proxy else "normal",
                 fontweight="bold")


def style_xticklabels(ax, labels=None):
    """IDR-proxy tools in italics; mechanism is shown by the category box instead."""
    labels = labels or pred_labels
    for tick, p in zip(ax.get_xticklabels(), labels):
        if TARGET_TYPE[p] == "IDR proxy":
            tick.set_fontstyle("italic")


def style_yticklabels(ax, labels):
    for tick, p in zip(ax.get_yticklabels(), labels):
        if TARGET_TYPE[p] == "IDR proxy":
            tick.set_fontstyle("italic")


def add_category_boxes(fig, axes):
    """Draw a labelled box around each mechanism category's block of subplots."""
    for (label, members), (r0, r1) in zip(CATEGORY_BLOCKS, BLOCK_ROWS):
        color = MECHANISM_COLOR[MECHANISM[members[0]]]
        p0, p1 = axes[r0, 0].get_position(), axes[r1, axes.shape[1] - 1].get_position()
        pad = 0.006
        x0, y1 = p0.x0 - pad, p0.y1 + pad
        x1, y0 = p1.x1 + pad, p1.y0 - pad
        fig.add_artist(mpatches.FancyBboxPatch(
            (x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0.002",
            linewidth=1.4, edgecolor=color, facecolor="none",
            transform=fig.transFigure, zorder=10, clip_on=False))
        fig.text(x0 + 0.004, y1 - 0.003, label, color=color, fontsize=8,
                 fontweight="bold", va="top", ha="left", zorder=11)


def add_category_boxes_bar(ax):
    """Draw a labelled box around each mechanism category's bars/ticks on the x-axis."""
    x = 0
    for label, members in CATEGORY_BLOCKS:
        n = len(members)
        color = MECHANISM_COLOR[MECHANISM[members[0]]]
        x0, x1 = x - 0.5, x + n - 0.5
        ax.add_patch(mpatches.Rectangle(
            (x0, 0), x1 - x0, 1, transform=ax.get_xaxis_transform(),
            linewidth=1.3, edgecolor=color, facecolor="none", zorder=5, clip_on=False))
        ax.text((x0 + x1) / 2, 0.99, CATEGORY_LABEL_SHORT[label], transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=7, color=color, fontweight="bold", zorder=6,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.8))
        x += n


# ── Shared helpers ────────────────────────────────────────────────────────────

def get_xy(df, label_col, filter_col, score_col, pos_mask_col=None):
    """Return (y_true, scores); optionally restrict positives to clean subset."""
    sub = df if filter_col is None else df[df[filter_col].notna()]
    if pos_mask_col is not None and pos_mask_col in sub.columns:
        is_neg   = sub[label_col] == 0
        is_clean = (sub[label_col] == 1) & sub[pos_mask_col].fillna(True).astype(bool)
        sub = sub[is_neg | is_clean]
    valid = sub[score_col].notna()
    return sub.loc[valid, label_col].astype(int).values, sub.loc[valid, score_col].values


def compute_auroc(y, scores):
    if len(np.unique(y)) < 2 or len(y) == 0:
        return np.nan
    fpr, tpr, _ = roc_curve(y, scores)
    return float(auc(fpr, tpr))


def mcc_from_rates(tpr, fpr, n_pos, n_neg):
    tp, fp = tpr * n_pos, fpr * n_neg
    tn, fn = n_neg - fp, n_pos - tp
    num   = tp * tn - fp * fn
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return num / np.where(denom == 0, 1e-12, denom)


def compute_max_mcc(y, scores):
    if len(np.unique(y)) < 2 or len(y) == 0:
        return np.nan
    n_pos, n_neg = int(y.sum()), int((y == 0).sum())
    fpr, tpr, _ = roc_curve(y, scores)
    return float(np.max(mcc_from_rates(tpr, fpr, n_pos, n_neg)))


def bootstrap_auroc_ci(y, scores, n=1000, seed=42):
    rng = np.random.default_rng(seed)
    aucs = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        yb, sb = y[idx], scores[idx]
        if len(np.unique(yb)) < 2:
            continue
        fpr, tpr, _ = roc_curve(yb, sb)
        aucs.append(auc(fpr, tpr))
    aucs = np.array(aucs)
    return (np.percentile(aucs, 2.5), np.percentile(aucs, 97.5)) if len(aucs) else (np.nan, np.nan)


def bootstrap_mcc_ci(y, scores, n=1000, seed=42):
    rng = np.random.default_rng(seed)
    mccs = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        yb, sb = y[idx], scores[idx]
        if len(np.unique(yb)) < 2:
            continue
        mccs.append(compute_max_mcc(yb, sb))
    mccs = np.array(mccs)
    return (np.percentile(mccs, 2.5), np.percentile(mccs, 97.5)) if len(mccs) else (np.nan, np.nan)


def scenario_legend():
    return [
        mlines.Line2D([], [], color=v[2], linewidth=2, linestyle=v[3], label=sc)
        for sc, v in SCENARIOS.items()
    ]


def bar_chart(ax, summary_df, value_col, err_lo_col, err_hi_col, ylabel, title, ref_line,
              scenario_subset=None):
    scs    = scenario_subset or sc_names
    sc_map = {k: v for k, v in SCENARIOS.items() if k in scs}
    x      = np.arange(len(pred_labels))
    bar_w  = 0.72 / len(scs)
    for i, sc in enumerate(scs):
        sub    = summary_df[summary_df["Scenario"] == sc].set_index("Predictor")
        vals   = [sub.loc[p, value_col]  if p in sub.index else np.nan for p in pred_labels]
        lo     = [sub.loc[p, err_lo_col] if p in sub.index else np.nan for p in pred_labels]
        hi     = [sub.loc[p, err_hi_col] if p in sub.index else np.nan for p in pred_labels]
        color  = sc_map[sc][2]
        offset = (i - len(scs) / 2 + 0.5) * bar_w
        ax.bar(x + offset, vals, bar_w, color=color, alpha=0.85,
               label=sc)
        ax.errorbar(x + offset, vals,
                    yerr=[[v - l if not np.isnan(v) else 0 for v, l in zip(vals, lo)],
                          [h - v if not np.isnan(v) else 0 for v, h in zip(vals, hi)]],
                    fmt="none", color="#333", linewidth=0.9, capsize=2)
    ax.axhline(ref_line, color="#333", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(pred_labels, rotation=35, ha="right", fontsize=9)
    style_xticklabels(ax)
    add_category_boxes_bar(ax)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8, ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0))


# ═════════════════════════════════════════════════════════════════════════════
# § 1  ROC curves + AUROC/AUPRC
# ═════════════════════════════════════════════════════════════════════════════

fig, axes = new_predictor_grid((3.5, 3.2), "ROC curves by scenario — LLPS predictors")

roc_records = []

for pred_label, score_col in PREDICTORS.items():
    r, c = GRID_POS[pred_label]
    ax = axes[r, c]
    ax.plot([0, 1], [0, 1], "k--", linewidth=0.7, alpha=0.35)

    for sc_name, (lc, fc, color, ls) in SCENARIOS.items():
        y, scores = get_xy(df, lc, fc, score_col)
        if len(np.unique(y)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y, scores)
        roc_auc      = auc(fpr, tpr)
        ap            = average_precision_score(y, scores)
        ci_lo, ci_hi  = bootstrap_auroc_ci(y, scores)
        n_pos, n_neg  = int(y.sum()), int((y == 0).sum())
        ax.plot(fpr, tpr, color=color, linewidth=1.4, linestyle=ls,
                label=f"{SCENARIO_TAG[sc_name]}  {roc_auc:.2f} [{ci_lo:.2f}–{ci_hi:.2f}]")
        roc_records.append({"Predictor": pred_label, "Scenario": sc_name,
                            "AUROC": roc_auc, "CI_lo": ci_lo, "CI_hi": ci_hi,
                            "AUPRC": ap, "n_pos": n_pos, "n_neg": n_neg})

    ax.set(xlim=(0, 1), ylim=(0, 1))
    style_predictor_title(ax, pred_label, fontsize=9)
    ax.set_xlabel("FPR", fontsize=8)
    ax.set_ylabel("TPR" if c == 0 else "", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=5.5, loc="lower right")

fig.legend(handles=scenario_legend(), loc="lower right", bbox_to_anchor=(0.99, 0.01),
           fontsize=7.5, title="Scenario", framealpha=0.9)
fig.tight_layout()
add_category_boxes(fig, axes)
fig.savefig(FIG_DIR / "roc_curves.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved roc_curves.png")

roc_summary = pd.DataFrame(roc_records)
x     = np.arange(len(pred_labels))
bar_w = 0.16

fig, ax = plt.subplots(figsize=(16, 5.5))
bar_chart(ax, roc_summary, "AUROC", "CI_lo", "CI_hi",
          "AUROC  (error bars = 95% bootstrap CI)",
          "AUROC by predictor and scenario", ref_line=0.5)
ax.set_ylim(0.3, 1.05)
fig.tight_layout()
fig.savefig(FIG_DIR / "auroc_summary.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved auroc_summary.png")

fig, ax = plt.subplots(figsize=(16, 5))
for i, sc in enumerate(sc_names):
    sub    = roc_summary[roc_summary["Scenario"] == sc].set_index("Predictor")
    vals   = [sub.loc[p, "AUPRC"] if p in sub.index else np.nan for p in pred_labels]
    color  = list(SCENARIOS.values())[i][2]
    offset = (i - len(sc_names) / 2 + 0.5) * bar_w
    ax.bar(x + offset, vals, bar_w, color=color, alpha=0.85,
           label=sc)

baselines = {
    sc: roc_summary[roc_summary["Scenario"] == sc].iloc[0]["n_pos"] /
        (roc_summary[roc_summary["Scenario"] == sc].iloc[0]["n_pos"] +
         roc_summary[roc_summary["Scenario"] == sc].iloc[0]["n_neg"])
    for sc in sc_names if not roc_summary[roc_summary["Scenario"] == sc].empty
}
for sc, bl in baselines.items():
    ax.axhline(bl, color=list(SCENARIOS.values())[sc_names.index(sc)][2],
               linewidth=0.8, linestyle=":", alpha=0.7)

ax.set_xticks(x)
ax.set_xticklabels(pred_labels, rotation=35, ha="right", fontsize=9)
style_xticklabels(ax)
add_category_boxes_bar(ax)
ax.set_ylabel("AUPRC")
ax.set_title("AUPRC by predictor and scenario (dotted lines = baseline positive rate)")
ax.legend(fontsize=8, ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0))
fig.tight_layout()
fig.savefig(FIG_DIR / "auprc_summary.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved auprc_summary.png")


# ═════════════════════════════════════════════════════════════════════════════
# § 2  MCC-vs-threshold curves + Max-MCC
# ═════════════════════════════════════════════════════════════════════════════

fig, axes = new_predictor_grid((3.5, 3.2), "MCC vs score threshold by scenario — LLPS predictors")

mcc_records = []

for pred_label, score_col in PREDICTORS.items():
    r, c = GRID_POS[pred_label]
    ax = axes[r, c]
    ax.axhline(0, color="#333", linestyle="--", linewidth=0.7, alpha=0.35)

    for sc_name, (lc, fc, color, ls) in SCENARIOS.items():
        y, scores = get_xy(df, lc, fc, score_col)
        if len(np.unique(y)) < 2:
            continue
        n_pos, n_neg     = int(y.sum()), int((y == 0).sum())
        fpr, tpr, thresholds = roc_curve(y, scores)
        mcc_vals         = mcc_from_rates(tpr, fpr, n_pos, n_neg)
        best_mcc         = float(np.max(mcc_vals))
        ci_lo, ci_hi     = bootstrap_mcc_ci(y, scores)

        s_min, s_max = scores.min(), scores.max()
        t_norm = (np.clip((thresholds - s_min) / (s_max - s_min), 0, 1)
                  if s_max > s_min else np.zeros_like(thresholds))

        ax.plot(t_norm, mcc_vals, color=color, linewidth=1.4, linestyle=ls,
                label=f"{SCENARIO_TAG[sc_name]}  {best_mcc:.2f}")
        mcc_records.append({"Predictor": pred_label, "Scenario": sc_name,
                            "MaxMCC": best_mcc, "CI_lo": ci_lo, "CI_hi": ci_hi,
                            "n_pos": n_pos, "n_neg": n_neg})

    ax.set(xlim=(0, 1), ylim=(-0.05, 1.0))
    style_predictor_title(ax, pred_label, fontsize=9)
    ax.set_xlabel("Norm. threshold", fontsize=7)
    ax.set_ylabel("MCC" if c == 0 else "", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=5.5, loc="upper left")

fig.legend(handles=scenario_legend(), loc="lower right", bbox_to_anchor=(0.99, 0.01),
           fontsize=7.5, title="Scenario", framealpha=0.9)
fig.tight_layout()
add_category_boxes(fig, axes)
fig.savefig(FIG_DIR / "mcc_curves.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved mcc_curves.png")

mcc_summary = pd.DataFrame(mcc_records)
fig, ax = plt.subplots(figsize=(16, 5.5))
bar_chart(ax, mcc_summary, "MaxMCC", "CI_lo", "CI_hi",
          "Max-MCC  (error bars = 95% bootstrap CI)",
          "Max-MCC by predictor and scenario", ref_line=0.0)
ax.set_ylim(-0.05, 0.8)
fig.tight_layout()
fig.savefig(FIG_DIR / "mcc_summary.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved mcc_summary.png")


# ═════════════════════════════════════════════════════════════════════════════
# § 3  Training-set contamination — per-tool-clean AUROC
# ═════════════════════════════════════════════════════════════════════════════

CLEAN_SCENARIOS = {k: v for k, v in SCENARIOS.items() if k != "All experimental vs proteome"}

clean_records = []

for pred_label, score_col in PREDICTORS.items():
    mask_col = f"clean_{score_col}"
    glob_col = "clean_global"

    for sc_name, (lc, fc, color, ls) in CLEAN_SCENARIOS.items():
        y_full,  s_full  = get_xy(df, lc, fc, score_col)
        y_clean, s_clean = get_xy(df, lc, fc, score_col, pos_mask_col=mask_col)
        y_glob,  s_glob  = get_xy(df, lc, fc, score_col, pos_mask_col=glob_col)

        auroc_full  = compute_auroc(y_full,  s_full)
        auroc_clean = compute_auroc(y_clean, s_clean)
        auroc_glob  = compute_auroc(y_glob,  s_glob)
        mcc_full    = compute_max_mcc(y_full,  s_full)
        mcc_clean   = compute_max_mcc(y_clean, s_clean)
        mcc_glob    = compute_max_mcc(y_glob,  s_glob)
        ci_lo, ci_hi = bootstrap_auroc_ci(y_clean, s_clean)
        delta = (auroc_clean - auroc_full) if not np.isnan(auroc_clean) else np.nan

        clean_records.append({
            "Predictor":     pred_label,
            "Scenario":      sc_name,
            "AUROC_full":    round(auroc_full,  3),
            "AUROC_clean":   round(auroc_clean, 3),
            "AUROC_global":  round(auroc_glob,  3),
            "MaxMCC_full":   round(mcc_full,    3),
            "MaxMCC_clean":  round(mcc_clean,   3),
            "MaxMCC_global": round(mcc_glob,    3),
            "Delta_AUROC":   round(delta,       3),
            "CI_lo":         round(ci_lo, 3),
            "CI_hi":         round(ci_hi, 3),
            "N_pos_full":    int(y_full.sum()),
            "N_pos_clean":   int(y_clean.sum()),
            "N_pos_global":  int(y_glob.sum()),
        })

clean_results = pd.DataFrame(clean_records)
clean_results.to_csv(ROOT / "output" / "clean_auroc_table.csv", index=False)
print("Saved output/clean_auroc_table.csv")

# Delta plot — scenario C, sorted by AUROC_full
prim     = clean_results[clean_results["Scenario"] == PRIMARY_SC].set_index("Predictor")
prim_ord = prim.reindex(pred_labels).sort_values("AUROC_full", ascending=True)
y_pos    = np.arange(len(prim_ord))
bar_h    = 0.35

fig, axes = plt.subplots(1, 2, figsize=(14, 7), gridspec_kw={"width_ratios": [3, 1]})
ax = axes[0]
ax.barh(y_pos + bar_h/2, prim_ord["AUROC_full"],  bar_h, color="#888888", alpha=0.75, label="Full (all 60)")
ax.barh(y_pos - bar_h/2, prim_ord["AUROC_clean"], bar_h, color="#e74c3c", alpha=0.85, label="Per-tool clean")
ax.errorbar(prim_ord["AUROC_clean"], y_pos - bar_h/2,
            xerr=[prim_ord["AUROC_clean"] - prim_ord["CI_lo"],
                  prim_ord["CI_hi"] - prim_ord["AUROC_clean"]],
            fmt="none", color="#333", linewidth=0.9, capsize=3)
ax.axvline(0.5, color="#333", linestyle="--", linewidth=0.8, alpha=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels([f"{p}  (n={int(n)})" for p, n in
                    zip(prim_ord.index, prim_ord["N_pos_clean"])], fontsize=9)
style_yticklabels(ax, prim_ord.index)
ax.set_xlabel("AUROC")
ax.set_title(f"Full vs per-tool clean AUROC — {PRIMARY_SC}")
ax.set_xlim(0.3, 1.05)
ax.legend(fontsize=9)

ax2 = axes[1]
deltas  = prim_ord["Delta_AUROC"].values
colours = ["#2ecc71" if d >= 0 else "#e74c3c" for d in deltas]
ax2.barh(y_pos, deltas, bar_h * 1.2, color=colours, alpha=0.8)
ax2.axvline(0, color="#333", linewidth=0.8)
ax2.set_yticks(y_pos)
ax2.set_yticklabels([])
ax2.set_xlabel("Δ AUROC\n(clean − full)")
ax2.set_title("Delta\n(green=rises, red=drops)")
for i, d in enumerate(deltas):
    ax2.text(d + (0.002 if d >= 0 else -0.002), i, f"{d:+.3f}", va="center",
             ha="left" if d >= 0 else "right", fontsize=7.5)
fig.tight_layout()
fig.savefig(FIG_DIR / "auroc_clean_delta.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved auroc_clean_delta.png")

# Per-tool-clean AUROC + Max-MCC across all scenarios
clean_sc_names = list(CLEAN_SCENARIOS.keys())
fig, axes = plt.subplots(1, 2, figsize=(32, 6))
for ax_idx, (value_col, ylabel) in enumerate([
    ("AUROC_clean",  "AUROC (per-tool clean)"),
    ("MaxMCC_clean", "Max-MCC (per-tool clean)"),
]):
    ax    = axes[ax_idx]
    x     = np.arange(len(pred_labels))
    bar_w = 0.72 / len(clean_sc_names)
    for i, sc in enumerate(clean_sc_names):
        sub    = clean_results[clean_results["Scenario"] == sc].set_index("Predictor")
        vals   = [sub.loc[p, value_col] if p in sub.index else np.nan for p in pred_labels]
        color  = CLEAN_SCENARIOS[sc][2]
        offset = (i - len(clean_sc_names) / 2 + 0.5) * bar_w
        ax.bar(x + offset, vals, bar_w, color=color, alpha=0.85,
               label=sc)
    ref = 0.5 if "AUROC" in value_col else 0.0
    ax.axhline(ref, color="#333", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(pred_labels, rotation=35, ha="right", fontsize=8.5)
    style_xticklabels(ax)
    add_category_boxes_bar(ax)
    ax.set_ylabel(ylabel)
    ax.set_ylim((0.3 if "AUROC" in value_col else -0.05), 1.05)
    ax.legend(fontsize=8, ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0))
fig.tight_layout()
fig.savefig(FIG_DIR / "auroc_clean_all.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved auroc_clean_all.png")

# ROC curves: full vs per-tool-clean, scenario C
sc_lc, sc_fc, sc_color, _ = CLEAN_SCENARIOS[PRIMARY_SC]
fig, axes = new_predictor_grid((3.4, 3.0), f"ROC curves: full vs per-tool clean — {PRIMARY_SC}")

for pred_label, score_col in PREDICTORS.items():
    r, c     = GRID_POS[pred_label]
    ax       = axes[r, c]
    mask_col = f"clean_{score_col}"
    ax.plot([0, 1], [0, 1], "k--", linewidth=0.6, alpha=0.3)

    y_full, s_full = get_xy(df, sc_lc, sc_fc, score_col)
    if len(np.unique(y_full)) >= 2:
        fpr_f, tpr_f, _ = roc_curve(y_full, s_full)
        ax.plot(fpr_f, tpr_f, color="#888888", linewidth=1.2, linestyle="--",
                label=f"Full  {auc(fpr_f, tpr_f):.2f}")

    y_c, s_c = get_xy(df, sc_lc, sc_fc, score_col, pos_mask_col=mask_col)
    if len(np.unique(y_c)) >= 2:
        fpr_c, tpr_c, _ = roc_curve(y_c, s_c)
        ax.plot(fpr_c, tpr_c, color=sc_color, linewidth=1.5,
                label=f"Clean {auc(fpr_c, tpr_c):.2f}  (n={int(y_c.sum())})")

    ax.set(xlim=(0, 1), ylim=(0, 1))
    style_predictor_title(ax, pred_label, fontsize=9)
    ax.set_xlabel("FPR", fontsize=7.5)
    ax.set_ylabel("TPR" if c == 0 else "", fontsize=7.5)
    ax.tick_params(labelsize=6.5)
    ax.legend(fontsize=6, loc="lower right")

fig.tight_layout()
add_category_boxes(fig, axes)
fig.savefig(FIG_DIR / "roc_curves_clean.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved roc_curves_clean.png")


# ═════════════════════════════════════════════════════════════════════════════
# § 3b  ROC curves, per-tool-clean, all background scenarios overlaid
# ═════════════════════════════════════════════════════════════════════════════

fig, axes = new_predictor_grid((3.5, 3.2), "ROC curves by scenario — per-tool training-leakage removed")

for pred_label, score_col in PREDICTORS.items():
    r, c     = GRID_POS[pred_label]
    ax       = axes[r, c]
    mask_col = f"clean_{score_col}"
    ax.plot([0, 1], [0, 1], "k--", linewidth=0.7, alpha=0.35)

    for sc_name, (lc, fc, color, ls) in CLEAN_SCENARIOS.items():
        y, scores = get_xy(df, lc, fc, score_col, pos_mask_col=mask_col)
        if len(np.unique(y)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y, scores)
        roc_auc = auc(fpr, tpr)
        ax.plot(fpr, tpr, color=color, linewidth=1.4, linestyle=ls,
                label=f"{SCENARIO_TAG[sc_name]}  {roc_auc:.2f}  (n={int(y.sum())})")

    ax.set(xlim=(0, 1), ylim=(0, 1))
    style_predictor_title(ax, pred_label, fontsize=9)
    ax.set_xlabel("FPR", fontsize=8)
    ax.set_ylabel("TPR" if c == 0 else "", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=5.5, loc="lower right")

clean_scenario_legend = [
    mlines.Line2D([], [], color=v[2], linewidth=2, linestyle=v[3], label=sc)
    for sc, v in CLEAN_SCENARIOS.items()
]
fig.legend(handles=clean_scenario_legend, loc="lower right", bbox_to_anchor=(0.99, 0.01),
           fontsize=7.5, title="Scenario (per-tool clean)", framealpha=0.9)
fig.tight_layout()
add_category_boxes(fig, axes)
fig.savefig(FIG_DIR / "roc_curves_clean_scenarios.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved roc_curves_clean_scenarios.png")


# ═════════════════════════════════════════════════════════════════════════════
# § 3c  ROC curves, per-tool-clean vs globally-clean — one figure per scenario
# ═════════════════════════════════════════════════════════════════════════════

GLOBAL_COLOR = "#222222"

for sc_name, (lc, fc, color, ls) in CLEAN_SCENARIOS.items():
    fig, axes = new_predictor_grid((3.4, 3.0), f"ROC curves: per-tool clean vs globally clean — {sc_name}")

    for pred_label, score_col in PREDICTORS.items():
        r, c     = GRID_POS[pred_label]
        ax       = axes[r, c]
        mask_col = f"clean_{score_col}"
        ax.plot([0, 1], [0, 1], "k--", linewidth=0.6, alpha=0.3)

        y_c, s_c = get_xy(df, lc, fc, score_col, pos_mask_col=mask_col)
        if len(np.unique(y_c)) >= 2:
            fpr_c, tpr_c, _ = roc_curve(y_c, s_c)
            ax.plot(fpr_c, tpr_c, color=color, linewidth=1.5,
                    label=f"Per-tool clean {auc(fpr_c, tpr_c):.2f}  (n={int(y_c.sum())})")

        y_g, s_g = get_xy(df, lc, fc, score_col, pos_mask_col="clean_global")
        if len(np.unique(y_g)) >= 2:
            fpr_g, tpr_g, _ = roc_curve(y_g, s_g)
            ax.plot(fpr_g, tpr_g, color=GLOBAL_COLOR, linewidth=1.2, linestyle="--",
                    label=f"Globally clean {auc(fpr_g, tpr_g):.2f}  (n={int(y_g.sum())})")

        ax.set(xlim=(0, 1), ylim=(0, 1))
        style_predictor_title(ax, pred_label, fontsize=9)
        ax.set_xlabel("FPR", fontsize=7.5)
        ax.set_ylabel("TPR" if c == 0 else "", fontsize=7.5)
        ax.tick_params(labelsize=6.5)
        ax.legend(fontsize=6, loc="lower right")

    fig.tight_layout()
    add_category_boxes(fig, axes)
    fname = f"roc_curves_clean_vs_global_{SCENARIO_SLUG[sc_name]}.png"
    fig.savefig(FIG_DIR / fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {fname}")


# ═════════════════════════════════════════════════════════════════════════════
# § 4  Matched-background sensitivity (100 bootstrap subsamples)
# ═════════════════════════════════════════════════════════════════════════════

N_BOOT = 100

pos_ids       = set(df.loc[df["label_mem"] == 1, "UniProt_ID"])
pool_proteome = df[(~df["is_membrane_protein"]) & (df["label_all"] == 0)]["UniProt_ID"].values
pool_membrane = df[df["is_membrane_protein"] & (df["label_all"] == 0) &
                   (~df["UniProt_ID"].isin(pos_ids))]["UniProt_ID"].values

print(f"Pools: proteome={len(pool_proteome)}, membrane={len(pool_membrane)}")


def matched_auroc_bootstrap(score_col, neg_pool, n_neg, n_boot=N_BOOT):
    pos_df  = df[df["UniProt_ID"].isin(pos_ids) & df[score_col].notna()]
    neg_df  = df[df["UniProt_ID"].isin(neg_pool) & df[score_col].notna()]
    neg_ids = neg_df["UniProt_ID"].values
    replace = len(neg_ids) < n_neg
    aucs, rng = [], np.random.default_rng(seed=42)
    for _ in range(n_boot):
        sampled = rng.choice(neg_ids, size=n_neg, replace=replace)
        sub     = pd.concat([pos_df, neg_df[neg_df["UniProt_ID"].isin(sampled)]])
        y       = sub["UniProt_ID"].isin(pos_ids).astype(int).values
        aucs.append(compute_auroc(y, sub[score_col].values))
    return np.array([a for a in aucs if not np.isnan(a)])


def full_auroc_sc(score_col, label_col, filter_col):
    sub = df[df[filter_col].notna()] if filter_col else df
    sub = sub[sub[score_col].notna()]
    return compute_auroc(sub[label_col].astype(int).values, sub[score_col].values)


MATCHED = {
    "Proteome background (1:1 matched)":  (pool_proteome,  60),
    "Proteome background (1:10 matched)": (pool_proteome, 600),
    "Membrane background (1:1 matched)":  (pool_membrane,  60),
    "Membrane background (1:10 matched)": (pool_membrane, 600),
}
REFERENCE = {
    "Proteome background (full)": ("label_mem_fixed", "label_mem_fixed"),
    "Membrane background (full)": ("label_mem_bg",    "label_mem_bg"),
}

matched_records = []
for pred, score_col in PREDICTORS.items():
    for ref_name, (lc, fc) in REFERENCE.items():
        auroc_ref = full_auroc_sc(score_col, lc, fc)
        n_neg_ref = int((df[fc] == 0).sum())
        matched_records.append({"Predictor": pred, "Scenario": ref_name,
                                 "AUROC_mean": auroc_ref, "AUROC_sd": 0.0,
                                 "AUROC_lo": auroc_ref, "AUROC_hi": auroc_ref,
                                 "n_neg": n_neg_ref})
    for sc_name, (pool, n_neg) in MATCHED.items():
        aucs = matched_auroc_bootstrap(score_col, pool, n_neg)
        if len(aucs) == 0:
            matched_records.append({"Predictor": pred, "Scenario": sc_name,
                                     "AUROC_mean": np.nan, "AUROC_sd": np.nan,
                                     "AUROC_lo": np.nan, "AUROC_hi": np.nan,
                                     "n_neg": n_neg})
        else:
            matched_records.append({"Predictor": pred, "Scenario": sc_name,
                                     "AUROC_mean": float(np.mean(aucs)),
                                     "AUROC_sd":   float(np.std(aucs)),
                                     "AUROC_lo":   float(np.percentile(aucs, 2.5)),
                                     "AUROC_hi":   float(np.percentile(aucs, 97.5)),
                                     "n_neg": n_neg})
    print(f"  {pred:<14} done")

matched_results = pd.DataFrame(matched_records).round(3)
matched_results.to_csv(ROOT / "output" / "matched_auroc_table.csv", index=False)
print("Saved output/matched_auroc_table.csv")

MATCHED_COLORS = {
    "Proteome background (full)":         "#aaaaaa",
    "Proteome background (1:1 matched)":  "#0072B2",
    "Proteome background (1:10 matched)": "#56B4E9",
    "Membrane background (full)":         "#aaaaaa",
    "Membrane background (1:1 matched)":  "#e74c3c",
    "Membrane background (1:10 matched)": "#f39c12",
}


def make_matched_fig(scenarios, title, fname):
    x     = np.arange(len(pred_labels))
    bar_w = 0.72 / len(scenarios)
    fig, ax = plt.subplots(figsize=(16, 5.5))
    for i, sc in enumerate(scenarios):
        sub    = matched_results[matched_results["Scenario"] == sc].set_index("Predictor")
        means  = [sub.loc[p, "AUROC_mean"] if p in sub.index else np.nan for p in pred_labels]
        lo     = [sub.loc[p, "AUROC_lo"]   if p in sub.index else np.nan for p in pred_labels]
        hi     = [sub.loc[p, "AUROC_hi"]   if p in sub.index else np.nan for p in pred_labels]
        offset = (i - len(scenarios) / 2 + 0.5) * bar_w
        hatch  = "//" if "full" in sc else ""
        ax.bar(x + offset, means, bar_w, color=MATCHED_COLORS[sc], alpha=0.85,
               label=sc, hatch=hatch,
               edgecolor="white" if not hatch else MATCHED_COLORS[sc])
        if "full" not in sc:
            ax.errorbar(x + offset, means,
                        yerr=[[m - l if not np.isnan(m) else 0 for m, l in zip(means, lo)],
                              [h - m if not np.isnan(m) else 0 for m, h in zip(means, hi)]],
                        fmt="none", color="#333", linewidth=0.9, capsize=3)
    ax.axhline(0.5, color="#333", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(pred_labels, rotation=35, ha="right", fontsize=9)
    style_xticklabels(ax)
    add_category_boxes_bar(ax)
    ax.set_ylabel("AUROC  (matched: mean ± 95% CI, 100 draws)")
    ax.set_ylim(0.3, 1.05)
    ax.set_title(title)
    ax.legend(fontsize=8.5, ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(FIG_DIR / fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {fname}")


make_matched_fig(
    ["Proteome background (full)", "Proteome background (1:1 matched)", "Proteome background (1:10 matched)"],
    "Matched vs full — proteome background (hatched = full)",
    "auroc_matched_proteome.png",
)
make_matched_fig(
    ["Membrane background (full)", "Membrane background (1:1 matched)", "Membrane background (1:10 matched)"],
    "Matched vs full — membrane background (hatched = full)",
    "auroc_matched_membrane.png",
)
