"""
STEP 8 -- figures for the background and threshold results.

These are the paper's opening figures under the current arc: the offset comes first,
then the background flip, then the cutoff disagreement.

Needs 06_background.py and 07_thresholds.py to have run.

Run:  python 08_figures_background.py
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "legend.frameon": False,
})

bgres = pd.read_csv(C.TABLES / "background_offset.csv")
m = pd.read_csv(C.BUILD / "master.csv")

C_BIAS = "#762a83"     # penalises membrane proteins
C_NEUT = "#b0b0b0"     # no membrane bias


def save(fig, name):
    out = C.FIGS / name
    fig.savefig(out)
    plt.close(fig)
    print(f"  wrote {out.name}")


# =============================================================== Figure 5
def fig_offset():
    """
    Two panels. Left: how far each tool's score alone separates membrane from soluble
    proteins, among proteins none of the databases call phase-separating. 0.5 is no
    bias. Right: that offset against how much the tool's apparent skill moves when
    you swap the background under it.
    """
    d = bgres.sort_values("offset_auroc")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.2),
                                   gridspec_kw={"width_ratios": [1.15, 1]})

    # ---- left: the offset. Drawn as a point per tool with a stem back to the 0.5
    # chance line, so the axis reads in real AUROC units and the mark's distance from
    # 0.5 is the bias. (A bar from the axis origin would imply a ratio to zero, which
    # is meaningless for an AUROC.)
    y = np.arange(len(d))
    cols = [C_BIAS if v < 0.5 else C_NEUT for v in d.offset_auroc]
    ax1.hlines(y, 0.5, d.offset_auroc, color=cols, lw=1.6)
    ax1.scatter(d.offset_auroc, y, s=26, color=cols, zorder=3)
    ax1.axvline(0.5, color="black", lw=0.8, zorder=2)
    ax1.set_yticks(y, d.tool)
    ax1.set_xlabel("membrane vs soluble AUROC\n(0.5 = no membrane bias)")
    ax1.invert_yaxis()
    ax1.margins(x=0.10)

    # ---- right: offset vs the shift it buys
    ax2.scatter(d.offset_auroc, d.auroc_shift, s=24, color="#333333", zorder=3)
    fit = np.polyfit(d.offset_auroc, d.auroc_shift, 1)
    xs = np.linspace(d.offset_auroc.min(), d.offset_auroc.max(), 50)
    ax2.plot(xs, np.polyval(fit, xs), color="#762a83", lw=1.2, zorder=2)

    rho = stats.spearmanr(d.offset_auroc, d.auroc_shift)
    ax2.set_xlabel("membrane vs soluble AUROC\n(0.5 = no membrane bias)")
    ax2.set_ylabel("AUROC gain from using a\nmembrane background")
    ax2.axhline(0, color="#999999", ls="--", lw=0.7, zorder=1)

    # Label the two ends and the one tool that runs the other way -- the extremes are
    # the informative points and a reader will look for them.
    for _, r in d.iterrows():
        if r.tool in ("PICNIC (GO)", "LLPhyScore", "SEG"):
            ax2.annotate(r.tool, (r.offset_auroc, r.auroc_shift),
                         textcoords="offset points", xytext=(6, 5), fontsize=7)
    ax2.margins(0.12)
    ax2.text(0.97, 0.95, f"rho = {rho.statistic:.2f}", transform=ax2.transAxes,
             ha="right", va="top", fontsize=7)

    fig.subplots_adjust(wspace=0.42)
    save(fig, "fig5_membrane_offset.png")


# =============================================================== Figure 6
def fig_background_flip():
    """
    The same 473 positives, the same scores, two different backgrounds. Each tool is
    a line between its two AUROC values. The dashed line is chance -- a tool whose
    line crosses it changes verdict on nothing more than the choice of comparison.
    """
    d = bgres.sort_values("auroc_shift")
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    y = np.arange(len(d))

    crosses = (d.auroc_vs_soluble < 0.5) & (d.auroc_vs_membrane > 0.5)

    for i, (_, r) in enumerate(d.iterrows()):
        col = C_BIAS if crosses.iloc[i] else "#888888"
        ax.plot([r.auroc_vs_soluble, r.auroc_vs_membrane], [i, i],
                color=col, lw=1.2, zorder=1)
    ax.scatter(d.auroc_vs_soluble, y, s=22, facecolor="white",
               edgecolor="#333333", linewidths=0.9, zorder=3)
    ax.scatter(d.auroc_vs_membrane, y, s=22, color="#333333", zorder=3)

    ax.axvline(0.5, color="#999999", ls="--", lw=0.8, zorder=0)
    ax.set_yticks(y, d.tool)
    ax.set_xlabel("AUROC for the same 473 membrane LLPS positives")
    ax.margins(y=0.02, x=0.06)

    ax.scatter([], [], s=22, facecolor="white", edgecolor="#333333",
               linewidths=0.9, label="vs soluble background")
    ax.scatter([], [], s=22, color="#333333", label="vs membrane background")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=2, fontsize=7)
    save(fig, "fig6_background_flip.png")


# =============================================================== Figure 7
def fig_pass_rates():
    """
    Pass rate at each tool's own published cutoff, with the cutoff printed against
    each bar. Sixteen-fold spread on identical proteins.
    """
    d = pd.read_csv(C.TABLES / "per_tool_pass_rates.csv").sort_values("pass_pct")
    fig, ax = plt.subplots(figsize=(4.4, 2.9))
    y = np.arange(len(d))
    ax.barh(y, d.pass_pct, color="#4a7fb5", height=0.66)
    for yi, (pct, cut, n) in enumerate(zip(d.pass_pct, d.cutoff, d.n_pass)):
        ax.text(pct + 1.2, yi, f"{n} pass  (cutoff {cut:g})",
                va="center", fontsize=7, color="#333333")
    ax.set_yticks(y, d.tool)
    ax.set_xlabel("percent of the membrane LLPS set called positive")
    ax.margins(x=0.30)
    save(fig, "fig7_pass_rates.png")


# =============================================================== Figure 8
def fig_agreement():
    """
    Left: how many of its available cutoffs each protein clears. Right: the same
    proteins on the consensus-rank axes, with the >=75% group picked out -- it
    separates on rank but not on cross-tool disagreement.
    """
    cutoffs = {c: t for c, t in C.PUBLISHED_CUTOFF.items() if c in m.columns}
    passes = pd.DataFrame({c: m[c] >= t for c, t in cutoffs.items()}).sum(axis=1)
    n_eval = pd.DataFrame({c: m[c].notna() for c in cutoffs}).sum(axis=1)
    frac = passes / n_eval

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0))

    strict = frac >= 0.75
    ax1.hist(passes[~strict], bins=np.arange(-0.5, 11.5), color="#c8c8c8")
    ax1.hist(passes[strict], bins=np.arange(-0.5, 11.5), color="#4a7fb5")
    ax1.set_xlabel("number of published cutoffs cleared")
    ax1.set_ylabel("proteins")
    ax1.set_xticks(range(0, 11, 2))

    ax2.scatter(m.mean_rank[~strict], m.rank_sd[~strict], s=7, color="#c8c8c8",
                linewidths=0, label="rest")
    ax2.scatter(m.mean_rank[strict], m.rank_sd[strict], s=16, color="#4a7fb5",
                linewidths=0, label=f"clears >=75% (n = {int(strict.sum())})")
    ax2.set_xlabel("consensus rank across 18 predictors")
    ax2.set_ylabel("cross-predictor disagreement (rank SD)")
    ax2.legend(loc="upper left", fontsize=7, markerscale=1.4)
    ax2.margins(0.04)

    fig.subplots_adjust(wspace=0.32)
    save(fig, "fig8_cutoff_agreement.png")


if __name__ == "__main__":
    fig_offset()
    fig_background_flip()
    fig_pass_rates()
    fig_agreement()
    print(f"\nfigures in {C.FIGS}")
