"""
Spaghetti / slope plot: per-predictor scores across topology regions.

Each protein is one set of connected points. X positions are
Whole → Cytoplasmic → Transmembrane → Extracellular/Lumenal.
Lines are coloured by whole-protein pLLPS score from the main dataset.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# ── build wide table from master CSV ──────────────────────────────────────────
REGION_APPROACH = {"PICNIC": "pdb_region", "PSPire": "pdb_region"}
WHOLE_APPROACH  = {"ESpritz": "whole_db", "PLAAC": "whole_db", "PScore": "whole_db",
                   "SEG": "whole_db", "catGRANULE": "whole_db", "PSPire": "whole_db"}
REGIONS = ["Cytoplasmic", "Transmembrane", "Extracellular/Lumenal"]

master = pd.read_csv("output/topology_scores_master.csv")
proteins = sorted(master.UniProt_ID.unique())
topo = pd.DataFrame({"UniProt_ID": proteins})

for tool in sorted(master.tool.unique()):
    df_t = master[master.tool == tool]
    # whole
    wa = WHOLE_APPROACH.get(tool, "whole")
    whole = df_t[df_t.approach == wa][["UniProt_ID", "score"]].drop_duplicates("UniProt_ID")
    if whole.empty:
        whole = df_t[df_t.approach == ("whole" if wa == "whole_db" else "whole_db")][
            ["UniProt_ID", "score"]].drop_duplicates("UniProt_ID")
    topo = topo.merge(whole.rename(columns={"score": f"{tool}_whole"}),
                      on="UniProt_ID", how="left")
    # regions
    ra = REGION_APPROACH.get(tool, "concatenated")
    for region in REGIONS:
        reg = df_t[(df_t.approach == ra) & (df_t.region == region)][
            ["UniProt_ID", "score"]].drop_duplicates("UniProt_ID")
        topo = topo.merge(reg.rename(columns={"score": f"{tool}_{region}"}),
                          on="UniProt_ID", how="left")

full = pd.read_csv("output/full_dataset.csv")[["Entry", "p(LLPS)"]]
topo = topo.merge(full, left_on="UniProt_ID", right_on="Entry", how="left")

TOOLS = sorted({
    c.rsplit("_", 1)[0]
    for c in topo.columns
    if any(c.endswith(f"_{r}") for r in REGIONS)
    and f"{c.rsplit('_',1)[0]}_whole" in topo.columns
    and c.rsplit("_", 1)[0] not in ("n",)
})

X_LABELS = ["Whole", "Cytoplasmic", "Transmembrane", "Extracellular /\nLumenal"]
X_POS    = [0, 1, 2, 3]

pllps_norm = mcolors.Normalize(vmin=0.4, vmax=1.0)
cmap       = cm.plasma

# ── layout ────────────────────────────────────────────────────────────────────
ncols = 4
nrows = int(np.ceil(len(TOOLS) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(20, nrows * 4.5))
axes_flat  = axes.flatten()

for ax_idx, tool in enumerate(TOOLS):
    ax   = axes_flat[ax_idx]
    wcol = f"{tool}_whole"
    rcols = [wcol] + [
        f"{tool}_{r}" if f"{tool}_{r}" in topo.columns else None
        for r in REGIONS
    ]

    for _, row in topo.iterrows():
        ys = []
        xs = []
        for xp, col in zip(X_POS, rcols):
            if col is None:
                continue
            v = row.get(col, np.nan)
            if pd.notna(v):
                ys.append(float(v))
                xs.append(xp)

        if len(xs) < 2:
            continue

        colour = cmap(pllps_norm(row["p(LLPS)"])) if pd.notna(row["p(LLPS)"]) else "#aaaaaa"
        ax.plot(xs, ys, color=colour, alpha=0.35, lw=0.9, zorder=2)
        ax.scatter(xs, ys, color=colour, s=14, alpha=0.6, zorder=3, linewidths=0)

    for xp, col in zip(X_POS, rcols):
        if col is None:
            continue
        med = topo[col].median()
        if pd.notna(med):
            ax.scatter([xp], [med], color="black", s=60, zorder=5,
                       marker="D", linewidths=0)

    med_xs, med_ys = [], []
    for xp, col in zip(X_POS, rcols):
        if col is None:
            continue
        med = topo[col].median()
        if pd.notna(med):
            med_xs.append(xp)
            med_ys.append(med)
    if len(med_xs) > 1:
        ax.plot(med_xs, med_ys, color="black", lw=2, zorder=4, ls="--", alpha=0.8)

    ra = REGION_APPROACH.get(tool, "concatenated")
    ax.set_title(f"{tool}\n({ra})", fontsize=9, fontweight="bold")
    ax.set_xticks([xp for xp, c in zip(X_POS, rcols) if c is not None])
    ax.set_xticklabels(
        [X_LABELS[i] for i, c in enumerate(rcols) if c is not None],
        fontsize=7.5
    )
    ax.tick_params(axis="y", labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.axhline(0, color="#cccccc", lw=0.7, ls=":", zorder=1)

for ax in axes_flat[len(TOOLS):]:
    ax.set_visible(False)

sm = cm.ScalarMappable(cmap=cmap, norm=pllps_norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes_flat[:len(TOOLS)], shrink=0.4, pad=0.02,
                    location="right", aspect=30)
cbar.set_label("Whole-protein p(LLPS)", fontsize=10)

from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([0], [0], color="black", lw=2, ls="--", marker="D",
                            markersize=6, label="Median")],
           loc="lower right", fontsize=9, frameon=False, bbox_to_anchor=(0.98, 0.01))

fig.suptitle(
    "LLPS predictor scores across topology regions\n"
    "Each line = one protein; colour = whole-protein p(LLPS); ◆ = median\n"
    "Region score: concatenated sequence (most tools) or AF2 PDB slice (PICNIC, PSPire)",
    fontsize=12, y=1.002
)

plt.tight_layout()
plt.savefig("output/figures/topology_spaghetti.png", dpi=150, bbox_inches="tight")
print("Saved output/figures/topology_spaghetti.png")
