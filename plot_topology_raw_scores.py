"""
Raw score strip and spaghetti — unified layout for all 14 tools.

X-axis is the same for every tool:
  x=0:   Whole protein
  x=1.5: Cytoplasmic — approach A (concatenated or pdb_region)
  x=2.0: Cytoplasmic — approach B (segment; absent for PDB tools)
  x=3.5: Transmembrane — approach A
  x=4.0: Transmembrane — approach B
  x=5.5: Extracellular/Lumenal — approach A
  x=6.0: Extracellular/Lumenal — approach B

Approach A = concatenated (sequence tools) or pdb_region (PICNIC, PSPire)
Approach B = max-segment per protein (sequence tools only)

Two output figures:
  output/figures/topology_raw_strip.png     — strip + box at each position
  output/figures/topology_raw_spaghetti.png — per-protein lines, colour = p(LLPS)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

REGIONS   = ["Cytoplasmic", "Transmembrane", "Extracellular/Lumenal"]

WHOLE_APPROACH = {"ESpritz": "whole_db", "PLAAC": "whole_db", "PScore": "whole_db",
                  "SEG": "whole_db", "catGRANULE": "whole_db", "PSPire": "whole_db"}

# Unified x positions
X_WHOLE = 0.0
XA = {"Cytoplasmic": 1.5, "Transmembrane": 3.5, "Extracellular/Lumenal": 5.5}
XB = {"Cytoplasmic": 2.0, "Transmembrane": 4.0, "Extracellular/Lumenal": 6.0}

XTICKS  = [0.0, 1.75, 3.75, 5.75]
XLABELS = ["Whole", "Cytoplasmic", "Transmembrane", "Extracellular /\nLumenal"]
XSEPS   = [2.75, 4.75]
XLIM    = (-0.5, 6.6)

COL = {
    "whole":       "#333333",
    "concatenated":"#3182bd",
    "pdb_region":  "#e6550d",
    "segment":     "#756bb1",
}

pllps_norm = mcolors.Normalize(vmin=0.4, vmax=1.0)
cmap_p     = cm.plasma

# ── load ──────────────────────────────────────────────────────────────────────
master = pd.read_csv("output/topology_scores_master.csv")
full   = pd.read_csv("output/full_dataset.csv")[["Entry","p(LLPS)"]].rename(
             columns={"Entry":"UniProt_ID"})

# Best whole-protein score per tool per protein
whole_rows = []
for tool in master.tool.unique():
    df_t = master[master.tool == tool]
    wa   = WHOLE_APPROACH.get(tool, "whole")
    w    = df_t[df_t.approach == wa][["UniProt_ID","score"]].drop_duplicates("UniProt_ID")
    if w.empty:
        alt = "whole" if wa == "whole_db" else "whole_db"
        w   = df_t[df_t.approach == alt][["UniProt_ID","score"]].drop_duplicates("UniProt_ID")
    w["tool"] = tool
    whole_rows.append(w)
whole_df = pd.concat(whole_rows).rename(columns={"score":"whole_score"})

# Max-segment per protein per tool per region
seg_max = (master[master.approach == "segment"]
           .groupby(["UniProt_ID","tool","region"])["score"].max()
           .reset_index().rename(columns={"score":"seg_max"}))

PDB_TOOLS = ["PICNIC","PSPire"]
SEQ_TOOLS = sorted(t for t in master.tool.unique() if t not in PDB_TOOLS)
ALL_TOOLS = SEQ_TOOLS + PDB_TOOLS


def _approach_a(tool):
    return "pdb_region" if tool in PDB_TOOLS else "concatenated"


def _get_region_a(df_t, tool, region):
    """Scores for approach A (concat or pdb_region) for a region."""
    return df_t[(df_t.approach == _approach_a(tool)) &
                (df_t.region == region)][["UniProt_ID","score"]].drop_duplicates("UniProt_ID")


def _get_region_b(tool, region):
    """Max-segment scores per protein for approach B."""
    sub = seg_max[(seg_max.tool == tool) & (seg_max.region == region)]
    return sub.rename(columns={"seg_max":"score"})[["UniProt_ID","score"]]


def _strip_panel(ax, tool, df_t, w_scores):
    """Draw strip+box on ax for a single tool."""
    rng = np.random.default_rng(hash(tool) % 99999)

    # whole
    vals = w_scores["whole_score"].dropna()
    jit  = rng.uniform(-0.10, 0.10, len(vals))
    ax.scatter(X_WHOLE + jit, vals, color=COL["whole"], s=16,
               alpha=0.55, linewidths=0, zorder=3)
    if len(vals) > 2:
        ax.boxplot(vals, positions=[X_WHOLE], widths=0.22, patch_artist=True,
                   manage_ticks=False,
                   boxprops=dict(facecolor=COL["whole"], alpha=0.25, lw=0.8),
                   medianprops=dict(color="white", lw=1.5),
                   whiskerprops=dict(lw=0.8, color=COL["whole"]),
                   capprops=dict(lw=0.8, color=COL["whole"]),
                   flierprops=dict(marker="", ms=0))

    for region in REGIONS:
        for xpos, get_fn, app_key in [
            (XA[region], lambda r=region: _get_region_a(df_t, tool, r), _approach_a(tool)),
            (XB[region], lambda r=region: _get_region_b(tool, r),       "segment"),
        ]:
            sub  = get_fn()
            vals2 = sub["score"].dropna()
            if len(vals2) == 0:
                continue
            col = COL[app_key]
            jit2 = rng.uniform(-0.09, 0.09, len(vals2))
            ax.scatter(xpos + jit2, vals2, color=col, s=14,
                       alpha=0.45, linewidths=0, zorder=3)
            if len(vals2) > 2:
                ax.boxplot(vals2, positions=[xpos], widths=0.20, patch_artist=True,
                           manage_ticks=False,
                           boxprops=dict(facecolor=col, alpha=0.25, lw=0.8),
                           medianprops=dict(color="black", lw=1.5),
                           whiskerprops=dict(lw=0.8, color=col),
                           capprops=dict(lw=0.8, color=col),
                           flierprops=dict(marker="", ms=0))

    for xsep in XSEPS:
        ax.axvline(xsep, color="#cccccc", lw=0.7, ls="--", zorder=1)
    ax.set_xticks(XTICKS)
    ax.set_xticklabels(XLABELS, fontsize=7.5)
    ax.set_xlim(*XLIM)
    ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_ylabel("Raw score", fontsize=7)
    app_a_label = "pdb_region / pdb_seg" if tool in PDB_TOOLS else "concat / seg"
    ax.set_title(f"{tool}  [{app_a_label}]", fontsize=9, fontweight="bold")


def _spag_panel(ax, tool, df_t, w_scores):
    """Draw spaghetti on ax for a single tool."""
    merged = w_scores.merge(full, on="UniProt_ID", how="left")

    for _, row in merged.iterrows():
        uid   = row["UniProt_ID"]
        ws    = row["whole_score"]
        pllps = row["p(LLPS)"]
        if pd.isna(ws):
            continue
        col = cmap_p(pllps_norm(pllps)) if pd.notna(pllps) else "#aaaaaa"

        # approach A path (solid)
        a_xs, a_ys = [X_WHOLE], [ws]
        for region in REGIONS:
            sub = _get_region_a(df_t, tool, region)
            v   = sub[sub.UniProt_ID == uid]["score"]
            if len(v) and pd.notna(v.iloc[0]):
                a_xs.append(XA[region])
                a_ys.append(float(v.iloc[0]))
        if len(a_xs) > 1:
            ax.plot(a_xs, a_ys, color=col, alpha=0.30, lw=0.8, ls="-", zorder=2)
            ax.scatter(a_xs, a_ys, color=col, s=13, alpha=0.55, linewidths=0, zorder=3)

        # approach B path (dashed)
        b_xs, b_ys = [X_WHOLE], [ws]
        for region in REGIONS:
            sub = _get_region_b(tool, region)
            v   = sub[sub.UniProt_ID == uid]["score"]
            if len(v) and pd.notna(v.iloc[0]):
                b_xs.append(XB[region])
                b_ys.append(float(v.iloc[0]))
        if len(b_xs) > 1:
            ax.plot(b_xs, b_ys, color=col, alpha=0.20, lw=0.7, ls="--", zorder=2)

    # medians — approach A (solid black)
    m_xs, m_ys = [X_WHOLE], [w_scores["whole_score"].median()]
    for region in REGIONS:
        vals = _get_region_a(df_t, tool, region)["score"].dropna()
        if len(vals):
            m_xs.append(XA[region])
            m_ys.append(vals.median())
    ax.plot(m_xs, m_ys, color="black", lw=2.0, ls="-", zorder=4, alpha=0.9)
    ax.scatter(m_xs, m_ys, color="black", s=55, marker="D", zorder=5, linewidths=0)

    # medians — approach B (dashed black)
    m_xs2, m_ys2 = [X_WHOLE], [w_scores["whole_score"].median()]
    for region in REGIONS:
        vals = _get_region_b(tool, region)["score"].dropna()
        if len(vals):
            m_xs2.append(XB[region])
            m_ys2.append(vals.median())
    if len(m_xs2) > 1:
        ax.plot(m_xs2, m_ys2, color="black", lw=1.5, ls="--", zorder=4, alpha=0.85)
        ax.scatter(m_xs2, m_ys2, color="black", s=40, marker="D",
                   zorder=5, linewidths=0)

    for xsep in XSEPS:
        ax.axvline(xsep, color="#cccccc", lw=0.7, ls=":", zorder=1)
    ax.set_xticks(XTICKS)
    ax.set_xticklabels(XLABELS, fontsize=7.5)
    ax.set_xlim(*XLIM)
    ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_ylabel("Raw score", fontsize=7)
    app_a_label = "pdb_region / pdb_seg" if tool in PDB_TOOLS else "concat / seg"
    ax.set_title(f"{tool}  [{app_a_label}]", fontsize=9, fontweight="bold")


# ── build figures ─────────────────────────────────────────────────────────────
ncols = 4
nrows = int(np.ceil(len(ALL_TOOLS) / ncols))

for fig_type in ["strip", "spaghetti"]:
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(22, nrows * (4.0 if fig_type=="strip" else 4.5)),
                             constrained_layout=True)
    flat = axes.flatten()

    for ax_i, tool in enumerate(ALL_TOOLS):
        ax   = flat[ax_i]
        df_t = master[master.tool == tool]
        w    = whole_df[whole_df.tool == tool]

        if fig_type == "strip":
            _strip_panel(ax, tool, df_t, w)
        else:
            _spag_panel(ax, tool, df_t, w)

    for ax in flat[len(ALL_TOOLS):]:
        ax.set_visible(False)

    legend_handles = [
        mpatches.Patch(color=COL["whole"],        label="Whole protein"),
        mpatches.Patch(color=COL["concatenated"],  label="Concatenated (isolated region, sequence tools)"),
        mpatches.Patch(color=COL["pdb_region"],    label="PDB region (AF2 structural slice; PICNIC, PSPire)"),
        mpatches.Patch(color=COL["segment"],       label="Max-segment per protein (sequence tools)"),
    ]
    fig.legend(handles=legend_handles, ncol=2, fontsize=8.5, frameon=False,
               loc="upper center", bbox_to_anchor=(0.5, 1.03))

    if fig_type == "strip":
        fig.suptitle(
            "Raw LLPS scores: whole protein vs topology region approaches\n"
            "Region group: left = concatenated/pdb_region; right = max-segment\n"
            "Segment points = one per individual topological span",
            fontsize=11, y=1.06
        )
    else:
        sm = cm.ScalarMappable(cmap=cmap_p, norm=pllps_norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=flat[:len(ALL_TOOLS)], shrink=0.3,
                            location="right", aspect=30, pad=0.01)
        cbar.set_label("Whole-protein p(LLPS)", fontsize=9)
        fig.legend(
            handles=legend_handles +
                    [Line2D([0],[0], color="black", lw=1.5, ls="-",  label="Median (concat/pdb)  ◆"),
                     Line2D([0],[0], color="black", lw=1.5, ls="--", label="Median (max-segment) ◆")],
            ncol=3, fontsize=8, frameon=False,
            loc="upper center", bbox_to_anchor=(0.5, 1.04)
        )
        fig.suptitle(
            "Raw LLPS scores per protein — solid = concatenated/pdb_region path; dashed = max-segment path",
            fontsize=11, y=1.07
        )

    out = f"output/figures/topology_raw_{fig_type}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")
    plt.close(fig)
