"""
STEP 5 -- the figures.

Four figures, one per idea, each built by its own function so you can run one
without the others. They go to minimal/figures/.

House rules followed here (they are this project's conventions):
  - no titles on the canvas. The claim goes in the caption you write, not the plot.
  - no commentary about how the figure was built.
  - where a number comes from a particular database, the database name goes in the
    axis label.
  - a colour key only where the colour is not obvious from the mark.

Run:  python 05_figures.py
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")            # write files, do not try to open a window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C


def style():
    """Set the look once. Three font sizes, thin spines, no top/right box."""
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "legend.frameon": False,
    })


style()

m = pd.read_csv(C.BUILD / "master.csv")
for col in C.INVERTED:                       # same sign flip as 04_stats.py
    if col in m.columns:
        m[col] = -m[col]

# Proteins that span no membrane are excluded from the topology comparisons.
mt = m[m.topology != "intramembrane-only"].copy()
COLOUR = {"single-pass": C.C_SINGLE, "multi-pass": C.C_MULTI}


def save(fig, name):
    out = C.FIGS / name
    fig.savefig(out)
    plt.close(fig)
    print(f"  wrote {out.name}")
    return out


# =============================================================== Figure 1
def fig_database_coverage():
    """
    How many proteins each database contributes, split by whether any other
    database also has it. The point is source overlap: a database whose bar is
    almost all 'shared' is not adding independent evidence to this set.
    """
    dbs = ["CD-CODE", "PhaSepDB", "DrLLPS", "PhaSePro", "LLPSDB"]
    uniq = [int((m[f"in_{d}"] & (m.n_databases == 1)).sum()) for d in dbs]
    shared = [int((m[f"in_{d}"] & (m.n_databases > 1)).sum()) for d in dbs]

    fig, ax = plt.subplots(figsize=(4.6, 2.5))
    y = np.arange(len(dbs))

    ax.barh(y, shared, color="#b0b0b0", height=0.62, label="also in another database")
    ax.barh(y, uniq, left=shared, color="#762a83", height=0.62,
            label="only source for this protein")

    # Label every segment with its count, including the dominant one, so the bar
    # can be read without consulting the key. A segment too short to hold its own
    # text gets the label outside instead, otherwise it is unreadable.
    span = max(np.array(shared) + np.array(uniq))
    for yi, s, u in zip(y, shared, uniq):
        if s and s > 0.08 * span:
            ax.text(s / 2, yi, f"{s} shared", ha="center", va="center",
                    color="white", fontsize=7)
            if u:
                ax.text(s + u + 0.012 * span, yi, f"{u} unique", va="center",
                        fontsize=7, color="#762a83")
        else:
            # Short bar: put the whole breakdown to the right of it.
            txt = f"{s} shared, {u} unique" if s else f"{u} unique"
            ax.text(s + u + 0.012 * span, yi, txt, va="center", fontsize=7,
                    color="#333333")

    ax.set_yticks(y, dbs)
    ax.invert_yaxis()
    ax.set_xlabel(f"membrane proteins in the set (n = {len(m)})")
    ax.legend(loc="lower right", fontsize=7)
    ax.margins(x=0.14)
    return save(fig, "fig1_database_coverage.png")


# =============================================================== Figure 2
def fig_topology():
    """
    Consensus predictor rank, single-pass vs multi-pass. Box for the quartiles,
    points for every protein so the reader sees the spread rather than a summary.
    """
    fig, ax = plt.subplots(figsize=(3.0, 2.8))
    groups = ["single-pass", "multi-pass"]
    data = [mt.loc[mt.topology == g, "mean_rank"].dropna() for g in groups]

    bp = ax.boxplot(data, positions=[0, 1], widths=0.5, showfliers=False,
                    patch_artist=True, medianprops=dict(color="black", lw=1.2))
    for patch, g in zip(bp["boxes"], groups):
        patch.set(facecolor=COLOUR[g], alpha=0.30, edgecolor=COLOUR[g], lw=0.8)

    rng = np.random.default_rng(0)
    for i, (g, d) in enumerate(zip(groups, data)):
        ax.scatter(i + rng.uniform(-0.16, 0.16, len(d)), d, s=5, alpha=0.45,
                   color=COLOUR[g], linewidths=0)

    ax.set_xticks([0, 1], [f"single-pass\nn = {len(data[0])}",
                           f"multi-pass\nn = {len(data[1])}"])
    ax.set_ylabel("consensus rank across 18 predictors")
    ax.margins(y=0.10)

    # The test result, as an effect size with its interval. Placed in headroom
    # opened above the data so it cannot land on a point or a whisker.
    s = pd.read_csv(C.TABLES / "headline_stats.csv").iloc[0]
    lo_y, hi_y = ax.get_ylim()
    ax.set_ylim(lo_y, hi_y + 0.20 * (hi_y - lo_y))
    ax.text(0.5, 0.99, f"rank-biserial {s.rank_biserial:+.2f} "
                       f"[{s.ci_lo:+.2f}, {s.ci_hi:+.2f}]\n"
                       f"P = {s.mannwhitney_p:.1e}",
            transform=ax.transAxes, ha="center", va="top", fontsize=7)
    return save(fig, "fig2_topology_consensus.png")


# =============================================================== Figure 3
def fig_hydropathy():
    """
    The confound, in one panel. Hydropathy on x, consensus rank on y, the two
    topology classes coloured. Both classes sit on ONE downward relationship --
    they differ in where along it they sit, not in having separate relationships.
    The heavy lines are per-class medians within hydropathy quintiles.
    """
    d = mt.dropna(subset=["mean_rank", "whole_hydropathy"])

    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    for g in ["single-pass", "multi-pass"]:
        sub = d[d.topology == g]
        ax.scatter(sub.whole_hydropathy, sub.mean_rank, s=7, alpha=0.40,
                   color=COLOUR[g], linewidths=0, label=g)

    edges = np.percentile(d.whole_hydropathy, [0, 20, 40, 60, 80, 100])
    centres = (edges[:-1] + edges[1:]) / 2
    for g in ["single-pass", "multi-pass"]:
        sub = d[d.topology == g]
        med = [sub.loc[(sub.whole_hydropathy >= lo) & (sub.whole_hydropathy <= hi),
                       "mean_rank"].median()
               for lo, hi in zip(edges[:-1], edges[1:])]
        ax.plot(centres, med, "-o", color=COLOUR[g], lw=1.6, ms=4,
                markeredgecolor="white", markeredgewidth=0.5)

    ax.set_xlabel("whole-protein hydropathy")
    ax.set_ylabel("consensus rank across 18 predictors")
    ax.legend(loc="upper right", fontsize=7, markerscale=1.6)
    ax.margins(0.04)
    return save(fig, "fig3_hydropathy_confound.png")


# =============================================================== Figure 4
def fig_per_predictor():
    """
    One row per predictor: the topology effect it shows on its own, with its 95%
    confidence interval. Filled where it survives multiple-testing correction,
    open where it does not. The dashed line is no effect.
    """
    per = pd.read_csv(C.TABLES / "per_predictor_topology.csv")
    per = per.sort_values("rank_biserial")

    fig, ax = plt.subplots(figsize=(4.0, 4.2))
    y = np.arange(len(per))
    sig = per.p_adj_BH < 0.05

    ax.axvline(0, color="#999999", ls="--", lw=0.7, zorder=0)
    ax.hlines(y, per.ci_lo, per.ci_hi, color="#444444", lw=1.0, zorder=1)
    ax.scatter(per.rank_biserial[sig], y[sig], s=26, color="#333333", zorder=2)
    ax.scatter(per.rank_biserial[~sig], y[~sig], s=26, facecolor="white",
               edgecolor="#333333", linewidths=0.9, zorder=2)

    ax.set_yticks(y, per.predictor)
    ax.set_xlabel("rank-biserial effect size")
    ax.margins(y=0.02)

    # Which end of the axis means what. Two short cues beat one long wrapped label
    # that runs off the figure.
    ax.annotate("ranks multi-pass higher", xy=(0.02, -0.105), xycoords="axes fraction",
                ha="left", va="top", fontsize=7, color="#555555")
    ax.annotate("ranks single-pass higher", xy=(0.98, -0.105), xycoords="axes fraction",
                ha="right", va="top", fontsize=7, color="#555555")

    # Filled vs open is the only encoding a reader cannot guess, so it is the only
    # thing keyed. Placed above the plot to stay clear of the SaPS/LLPhyScore rows.
    ax.scatter([], [], s=26, color="#333333", label="P < 0.05 after BH correction")
    ax.scatter([], [], s=26, facecolor="white", edgecolor="#333333",
               linewidths=0.9, label="not significant")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=2, fontsize=7)
    return save(fig, "fig4_per_predictor_topology.png")


if __name__ == "__main__":
    fig_database_coverage()
    fig_topology()
    fig_hydropathy()
    fig_per_predictor()
    print(f"\nfigures in {C.FIGS}")
