"""
plot_distributions.py

Score distributions and summary comparisons for all LLPS predictors:
  all experimental proteins (n=882, PhasePDB + LLPSDB) vs membrane subset (n=60).

Figures — output/figures/
  distributions_kde.png               Two-row KDE: main predictors / PSPspredict sub-scores
  distributions_all_vs_membrane.png   KDE grid (4×4, all predictors)
  median_comparison_all_vs_membrane.png  Median + IQR bar chart, Mann-Whitney p
  membrane_score_vs_fuzdrop.png       Predictor score vs FuzDrop p(LLPS) scatter
"""

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde, mannwhitneyu, spearmanr

sns.set_style("whitegrid")
plt.rcParams.update({"font.size": 10, "axes.labelsize": 11,
                     "xtick.labelsize": 9, "ytick.labelsize": 9})

ROOT    = Path(__file__).parent
FIG_DIR = ROOT / "output" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

df_all = pd.read_csv(ROOT / "output" / "predictor_comparison_all.csv")
df_mem = pd.read_csv(ROOT / "output" / "predictor_comparison_mem.csv")

MAIN = {
    "PICNIC":      "PICNIC_score",
    "PICNIC (GO)": "PICNIC_GO_score",
    "PSAP":        "PSAP_score",
    "PSPHunter":   "PSPHunter_prob",
    "PSPire":      "PSPire_score",
}
SUBS = {
    "PDL":        "PSPspredict_PDL",
    "catGRANULE": "PSPspredict_catGRANULE",
    "PLAAC":      "PSPspredict_PLAAC",
    "PScore":     "PSPspredict_PScore",
    "ESpritz":    "PSPspredict_ESpritz_DisProt",
    "SEG":        "PSPspredict_SEG",
    "DeepPhase":  "PSPspredict_DeepPhase",
    "SaPS":       "PSPspredict_SaPS_10fea",
    "PdPS":       "PSPspredict_PdPS_10fea",
    "PSPHunter*": "PSPspredict_Psphunter",
    "PSPire*":    "PSPspredict_PSPire",
}
PREDICTORS = {**MAIN, **SUBS}

for _df in (df_all, df_mem):
    for _c in PREDICTORS.values():
        if _c in _df.columns:
            _df[_c] = pd.to_numeric(_df[_c], errors="coerce")

COL_ALL = "#888888"
COL_MEM = "#0072B2"

legend_handles = [
    mpatches.Patch(color=COL_ALL, alpha=0.7, label=f"All experimental (n={len(df_all)})"),
    mpatches.Patch(color=COL_MEM, alpha=0.7, label=f"Membrane subset (n={len(df_mem)})"),
]

n_pred = len(PREDICTORS)
ncols  = 4
nrows  = int(np.ceil(n_pred / ncols))


# ── Shared helper ─────────────────────────────────────────────────────────────

def kde_plot(ax, vals_all, vals_mem, label, shaded=False):
    for vals, color in [(vals_all, COL_ALL), (vals_mem, COL_MEM)]:
        v = vals.dropna().to_numpy(dtype=float)
        if len(v) < 4:
            continue
        kde = gaussian_kde(v, bw_method="scott")
        pad = (v.max() - v.min()) * 0.05
        x   = np.linspace(max(0, v.min() - pad), min(1, v.max() + pad), 400)
        ax.fill_between(x, kde(x), alpha=0.35, color=color)
        ax.plot(x, kde(x), color=color, linewidth=1.6)
        ax.axvline(float(np.median(v)), color=color, linestyle="--", linewidth=1.1, alpha=0.85)
    a = vals_all.dropna().to_numpy(dtype=float)
    m = vals_mem.dropna().to_numpy(dtype=float)
    if len(a) >= 4 and len(m) >= 4:
        _, p = mannwhitneyu(a, m, alternative="two-sided")
        sig  = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
        ax.text(0.97, 0.97, sig, transform=ax.transAxes, ha="right", va="top", fontsize=10,
                color="#c0392b" if sig != "ns" else "#888888", fontweight="bold")
    ax.set_title(label, fontsize=10)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Score", fontsize=8)
    ax.tick_params(labelsize=8)
    ax.set_yticks([])
    if shaded:
        ax.set_facecolor("#f7f7f7")


# ═════════════════════════════════════════════════════════════════════════════
# § 1  Two-row KDE: main predictors / PSPspredict sub-scores
# ═════════════════════════════════════════════════════════════════════════════

ncols_k = max(len(MAIN), len(SUBS))
fig = plt.figure(figsize=(ncols_k * 3.0, 7.5))

axes_main = [fig.add_subplot(2, ncols_k, i + 1)          for i in range(len(MAIN))]
axes_subs = [fig.add_subplot(2, ncols_k, ncols_k + i + 1) for i in range(len(SUBS))]

for ax, (label, col) in zip(axes_main, MAIN.items()):
    kde_plot(ax, df_all[col], df_mem[col], label)
axes_main[0].set_ylabel("Density", fontsize=9)

for ax, (label, col) in zip(axes_subs, SUBS.items()):
    kde_plot(ax, df_all[col], df_mem[col], label, shaded=True)
axes_subs[0].set_ylabel("Density", fontsize=9)

fig.text(0.01, 0.97, "Main predictors",        fontsize=10, fontweight="bold", va="top")
fig.text(0.01, 0.50, "PSPspredict sub-scores", fontsize=10, fontweight="bold", va="top")
fig.legend(handles=legend_handles + [plt.Line2D([0], [0], color="#555", linestyle="--",
           linewidth=1.2, label="Median")],
           loc="upper right", bbox_to_anchor=(0.99, 0.99), fontsize=9, framealpha=0.9)
fig.suptitle("Score distributions: all experimental LLPS vs membrane subset",
             fontsize=12, y=1.02)
plt.tight_layout()
fig.savefig(FIG_DIR / "distributions_kde.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved distributions_kde.png")


# ═════════════════════════════════════════════════════════════════════════════
# § 2  KDE grid — all predictors (4×4)
# ═════════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.4, nrows * 2.8))
axes_flat = axes.flatten()
fig.suptitle(
    "Score distributions: all experimental LLPS (grey, n=882) vs membrane subset (blue, n=60)",
    fontsize=11, y=1.01,
)

for idx, (label, col) in enumerate(PREDICTORS.items()):
    ax       = axes_flat[idx]
    vals_all = df_all[col].dropna()
    vals_mem = df_mem[col].dropna()
    for vals, color, zo in [(vals_all, COL_ALL, 1), (vals_mem, COL_MEM, 2)]:
        if len(vals) < 3:
            continue
        v   = vals.to_numpy(dtype=float)
        kde = gaussian_kde(v, bw_method="scott")
        x   = np.linspace(v.min(), v.max(), 300)
        ax.fill_between(x, kde(x), alpha=0.42, color=color, zorder=zo)
        ax.plot(x, kde(x), color=color, linewidth=1.4, zorder=zo + 1)
    ax.axvline(vals_all.median(), color=COL_ALL, linestyle="--", linewidth=1.0, alpha=0.9)
    ax.axvline(vals_mem.median(), color=COL_MEM, linestyle="--", linewidth=1.0, alpha=0.9)
    if len(vals_all) >= 3 and len(vals_mem) >= 3:
        _, p = mannwhitneyu(vals_all, vals_mem, alternative="two-sided")
        pstr = f"p={p:.3f}" if p >= 0.001 else "p<0.001"
        ax.text(0.97, 0.96, pstr, transform=ax.transAxes,
                ha="right", va="top", fontsize=7.5, color="#333333")
    ax.set_title(label)
    ax.set_xlabel("Score")
    ax.set_ylabel("Density" if idx % ncols == 0 else "")
    ax.set_xlim(left=0)

for idx in range(n_pred, len(axes_flat)):
    axes_flat[idx].set_visible(False)

fig.legend(handles=legend_handles, loc="lower right", bbox_to_anchor=(0.99, 0.01), fontsize=9)
fig.tight_layout()
fig.savefig(FIG_DIR / "distributions_all_vs_membrane.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved distributions_all_vs_membrane.png")


# ═════════════════════════════════════════════════════════════════════════════
# § 3  Median + IQR summary, Mann-Whitney significance
# ═════════════════════════════════════════════════════════════════════════════

labels = list(PREDICTORS.keys())
cols   = list(PREDICTORS.values())
x      = np.arange(len(labels))
width  = 0.35

fig, ax = plt.subplots(figsize=(14, 5))
for offset, df, color, group_label in [
    (-width / 2, df_all, COL_ALL, f"All experimental (n={len(df_all)})"),
    ( width / 2, df_mem, COL_MEM, f"Membrane (n={len(df_mem)})"),
]:
    medians = [df[c].median()                          for c in cols]
    yerr_lo = [df[c].median() - df[c].quantile(0.25)  for c in cols]
    yerr_hi = [df[c].quantile(0.75) - df[c].median()  for c in cols]
    ax.bar(x + offset, medians, width, color=color, alpha=0.8, label=group_label, zorder=2)
    ax.errorbar(x + offset, medians, yerr=[yerr_lo, yerr_hi],
                fmt="none", color="#333333", linewidth=1.0, capsize=3, zorder=3)

for i, col in enumerate(cols):
    a = df_all[col].dropna()
    m = df_mem[col].dropna()
    if len(a) >= 3 and len(m) >= 3:
        _, p = mannwhitneyu(a, m, alternative="two-sided")
        sig  = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else ""))
        if sig:
            y_top = min(max(df_all[col].quantile(0.75), df_mem[col].quantile(0.75)) + 0.04, 1.02)
            ax.text(i, y_top, sig, ha="center", va="bottom", fontsize=9, color="#c0392b")

ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=35, ha="right")
ax.set_ylabel("Median score (bars = IQR)")
ax.set_ylim(0, 1.08)
ax.set_title("Median predictor scores: all experimental vs membrane subset\n"
             "(* p<0.05  ** p<0.01  *** p<0.001, Mann-Whitney)")
ax.legend(fontsize=9)
ax.axvline(4.5, color="#cccccc", linewidth=1, linestyle=":")
ax.text(4.6, 1.05, "PSPspredict sub-scores →", fontsize=7.5, color="#888888", va="top")
fig.tight_layout()
fig.savefig(FIG_DIR / "median_comparison_all_vs_membrane.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved median_comparison_all_vs_membrane.png")


# ═════════════════════════════════════════════════════════════════════════════
# § 4  Membrane only: predictor score vs FuzDrop p(LLPS)
# ═════════════════════════════════════════════════════════════════════════════

CLASS_ORDER  = ["High", "Medium", "Low"]
CLASS_COLORS = {"High": "#e74c3c", "Medium": "#f39c12", "Low": "#3498db"}

fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.4, nrows * 3.0))
axes_flat = axes.flatten()
fig.suptitle(
    "Membrane proteins: predictor score vs FuzDrop p(LLPS) (coloured by pLLPS class)",
    fontsize=11, y=1.01,
)

for idx, (label, col) in enumerate(PREDICTORS.items()):
    ax  = axes_flat[idx]
    sub = df_mem[["p(LLPS)", col, "pLLPS_Class"]].dropna()
    for cls in CLASS_ORDER:
        s = sub[sub["pLLPS_Class"] == cls]
        ax.scatter(s["p(LLPS)"], s[col], c=CLASS_COLORS[cls],
                   label=cls, s=28, alpha=0.8, edgecolors="none", zorder=2)
    if len(sub) >= 5:
        r, p = spearmanr(sub["p(LLPS)"], sub[col])
        pstr = "p<0.001" if p < 0.001 else f"p={p:.3f}"
        ax.text(0.03, 0.97, f"ρ={r:.2f} {pstr}", transform=ax.transAxes,
                va="top", fontsize=7.5, color="#333333")
    ax.set_title(label)
    ax.set_xlabel("FuzDrop p(LLPS)", fontsize=8)
    ax.set_ylabel("Score", fontsize=8)

for idx in range(n_pred, len(axes_flat)):
    axes_flat[idx].set_visible(False)

handles = [plt.Line2D([0], [0], marker="o", color="w",
                       markerfacecolor=CLASS_COLORS[c], markersize=7, label=c)
           for c in CLASS_ORDER]
fig.legend(handles=handles, title="pLLPS class", loc="lower right",
           bbox_to_anchor=(0.99, 0.01), fontsize=9)
fig.tight_layout()
fig.savefig(FIG_DIR / "membrane_score_vs_fuzdrop.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved membrane_score_vs_fuzdrop.png")


# ── Console summary ───────────────────────────────────────────────────────────

print("\n── Median scores + Mann-Whitney p ────────────────────────────────────")
print(f"{'Predictor':<14} {'All median':>11} {'Mem median':>11} {'p-value':>10} {'sig':>5}")
print("-" * 57)
for label, col in PREDICTORS.items():
    a = df_all[col].dropna()
    m = df_mem[col].dropna()
    _, p = (mannwhitneyu(a, m, alternative="two-sided")
            if len(a) >= 3 and len(m) >= 3 else (None, float("nan")))
    sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
    print(f"{label:<14} {a.median():>11.3f} {m.median():>11.3f} {p:>10.4f} {sig:>5}")
