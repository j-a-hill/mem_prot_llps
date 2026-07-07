"""
Consensus vote heatmap for topology region scores.

For each tool × protein × region, vote +1 if region score > whole-protein score,
-1 if lower, grey if missing. Three panels (Cytoplasmic / Transmembrane /
Extracellular+Lumenal) share the same protein row order.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec

ROOT_OUT = "output/figures"

REGIONS       = ["Cytoplasmic", "Transmembrane", "Extracellular/Lumenal"]
REGION_LABELS = ["Cytoplasmic", "Transmembrane", "Extracellular /\nLumenal"]

# ── build wide table from master CSV ──────────────────────────────────────────
REGION_APPROACH = {"PICNIC": "pdb_region", "PSPire": "pdb_region"}
WHOLE_APPROACH  = {"ESpritz": "whole_db", "PLAAC": "whole_db", "PScore": "whole_db",
                   "SEG": "whole_db", "catGRANULE": "whole_db", "PSPire": "whole_db"}

master = pd.read_csv("output/topology_scores_master.csv")
proteins = sorted(master.UniProt_ID.unique())
topo = pd.DataFrame({"UniProt_ID": proteins})

for tool in sorted(master.tool.unique()):
    df_t = master[master.tool == tool]
    wa = WHOLE_APPROACH.get(tool, "whole")
    whole = df_t[df_t.approach == wa][["UniProt_ID", "score"]].drop_duplicates("UniProt_ID")
    if whole.empty:
        whole = df_t[df_t.approach == ("whole" if wa == "whole_db" else "whole_db")][
            ["UniProt_ID", "score"]].drop_duplicates("UniProt_ID")
    topo = topo.merge(whole.rename(columns={"score": f"{tool}_whole"}),
                      on="UniProt_ID", how="left")
    ra = REGION_APPROACH.get(tool, "concatenated")
    for region in REGIONS:
        reg = df_t[(df_t.approach == ra) & (df_t.region == region)][
            ["UniProt_ID", "score"]].drop_duplicates("UniProt_ID")
        topo = topo.merge(reg.rename(columns={"score": f"{tool}_{region}"}),
                          on="UniProt_ID", how="left")

meta = pd.read_csv("output/master_table.csv")[["UniProt_ID", "Gene_Names"]]
meta["gene"] = meta["Gene_Names"].str.split().str[0]
topo = topo.merge(meta[["UniProt_ID", "gene"]], on="UniProt_ID", how="left")

# ── identify tools with both whole and region columns ─────────────────────────
TOOLS = sorted({
    c.rsplit("_", 1)[0]
    for c in topo.columns
    if any(c.endswith(f"_{r}") for r in REGIONS)
    and f"{c.rsplit('_', 1)[0]}_whole" in topo.columns
    and c.rsplit("_", 1)[0] not in ("n",)
})

# ── build vote matrices: proteins × tools, one per region ─────────────────────
def vote_matrix(region):
    mat = pd.DataFrame(index=topo["UniProt_ID"], columns=TOOLS, dtype=float)
    for tool in TOOLS:
        wcol = f"{tool}_whole"
        rcol = f"{tool}_{region}"
        if rcol not in topo.columns:
            continue
        delta = topo[rcol].values - topo[wcol].values
        mat[tool] = np.where(np.isnan(delta), np.nan, np.sign(delta))
    return mat

votes = {r: vote_matrix(r) for r in REGIONS}

# ── sort proteins by cytoplasmic positive-vote count (descending) ─────────────
cyto_pos = (votes["Cytoplasmic"] == 1).sum(axis=1)
protein_order = cyto_pos.sort_values(ascending=False).index.tolist()
gene_labels = (
    topo.set_index("UniProt_ID").loc[protein_order, "gene"]
    .fillna(pd.Series(protein_order, index=protein_order))
    .tolist()
)

# ── sort tools by cytoplasmic positive-vote fraction (descending) ─────────────
cyto_frac = (votes["Cytoplasmic"] == 1).mean()
tool_order = cyto_frac.sort_values(ascending=False).index.tolist()

# ── colormap: -1=blue, 0=lightgrey (NaN), +1=red ─────────────────────────────
cmap = ListedColormap(["#3182bd", "#d4d4d4", "#de2d26"])

def to_plot_array(mat, prot_order, tool_order):
    m = mat.loc[prot_order, tool_order].values.astype(float)
    out = np.full(m.shape, 1.0)    # grey = NaN
    out[m == -1] = 0.0             # blue
    out[m ==  1] = 2.0             # red
    return out

def frac_pos(mat, prot_order):
    sub = mat.loc[prot_order]
    return (sub == 1).sum(axis=1) / sub.notna().sum(axis=1)

# ── figure ────────────────────────────────────────────────────────────────────
n_prot  = len(protein_order)
n_tools = len(tool_order)

fig = plt.figure(figsize=(max(22, n_tools * 1.3), 14))
fig.suptitle("Consensus vote: does each tool score the region higher than the whole protein?",
             fontsize=13, y=0.98)

gs = GridSpec(1, 4, figure=fig,
              width_ratios=[n_tools, n_tools, n_tools, 1.5],
              wspace=0.08)

axes_heat = [fig.add_subplot(gs[0, i]) for i in range(3)]
ax_bar    = fig.add_subplot(gs[0, 3])

for ax, region, region_label in zip(axes_heat, REGIONS, REGION_LABELS):
    arr = to_plot_array(votes[region], protein_order, tool_order)
    ax.imshow(arr, aspect="auto", cmap=cmap, vmin=0, vmax=2,
              interpolation="none")

    # add region approach annotation to x-axis label
    ra = REGION_APPROACH.get(tool_order[0], "concatenated") if tool_order else "concatenated"
    ax.set_xticks(range(n_tools))
    # annotate per-tool approach in label
    xlabels = [
        f"{t}\n(pdb)" if REGION_APPROACH.get(t) == "pdb_region" else t
        for t in tool_order
    ]
    ax.set_xticklabels(xlabels, rotation=40, ha="right", fontsize=7)
    ax.set_title(region_label, fontsize=11, pad=6)

    if ax is axes_heat[0]:
        ax.set_yticks(range(n_prot))
        ax.set_yticklabels(gene_labels, fontsize=6.5)
        ax.set_ylabel("Protein (gene name)", fontsize=9)
    else:
        ax.set_yticks([])

    frac = frac_pos(votes[region], protein_order).values
    for i, f in enumerate(frac):
        if not np.isnan(f):
            ax.text(n_tools - 0.45, i, f"{f:.0%}", va="center", ha="left",
                    fontsize=5, color="black")

# ── right bar: cytoplasmic positive-vote fraction per protein ─────────────────
frac_cyto = frac_pos(votes["Cytoplasmic"], protein_order).values
colors_bar = ["#de2d26" if f >= 0.5 else "#3182bd" for f in frac_cyto]
ax_bar.barh(range(n_prot), frac_cyto, color=colors_bar, height=0.8)
ax_bar.axvline(0.5, color="black", lw=0.8, ls="--")
ax_bar.set_xlim(0, 1)
ax_bar.set_ylim(-0.5, n_prot - 0.5)
ax_bar.invert_yaxis()
ax_bar.set_yticks([])
ax_bar.set_xticks([0, 0.5, 1])
ax_bar.set_xticklabels(["0%", "50%", "100%"], fontsize=7)
ax_bar.set_title("Cyto\nvote %", fontsize=9, pad=6)
ax_bar.spines[["top", "right", "left"]].set_visible(False)

legend_patches = [
    mpatches.Patch(color="#de2d26", label="Region > Whole (+)"),
    mpatches.Patch(color="#3182bd", label="Region < Whole (−)"),
    mpatches.Patch(color="#d4d4d4", label="Missing data"),
]
fig.legend(handles=legend_patches, loc="lower center", ncol=3,
           fontsize=9, frameon=False, bbox_to_anchor=(0.45, -0.01))

plt.savefig(f"{ROOT_OUT}/topology_consensus_vote.png",
            dpi=150, bbox_inches="tight")
print(f"Saved {ROOT_OUT}/topology_consensus_vote.png")
