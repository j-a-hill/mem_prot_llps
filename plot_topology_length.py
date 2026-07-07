"""
Does region length drive LLPS predictions?

Panel A: scatter of region length vs region score, faceted by tool × region
Panel B: heatmap of Spearman r between region length and region score, per tool
Panel C: non-TMD region length vs whole-protein score, per tool
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from scipy.stats import spearmanr

# ── build wide score table from master CSV ────────────────────────────────────
REGION_APPROACH = {"PICNIC": "pdb_region", "PSPire": "pdb_region"}
WHOLE_APPROACH  = {"ESpritz": "whole_db", "PLAAC": "whole_db", "PScore": "whole_db",
                   "SEG": "whole_db", "catGRANULE": "whole_db", "PSPire": "whole_db"}
REGIONS    = ["Cytoplasmic", "Transmembrane", "Extracellular/Lumenal"]

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

# ── add region length columns from annotation-based wide file ─────────────────
lengths = pd.read_csv("output/topology_domain_scores.csv")[
    ["UniProt_ID", "n_Cytoplasmic", "n_Transmembrane", "n_Extracellular/Lumenal"]
]
topo = topo.merge(lengths, on="UniProt_ID", how="left")

meta = pd.read_csv("output/master_table.csv")[["UniProt_ID", "Gene_Names"]]
meta["gene"] = meta["Gene_Names"].str.split().str[0]
topo = topo.merge(meta[["UniProt_ID", "gene"]], on="UniProt_ID", how="left")

LEN_COLS   = ["n_Cytoplasmic", "n_Transmembrane", "n_Extracellular/Lumenal"]
REG_LABELS = ["Cytoplasmic", "Transmembrane", "Extracellular /\nLumenal"]

TOOLS = sorted({
    c.rsplit("_", 1)[0]
    for c in topo.columns
    if any(c.endswith(f"_{r}") for r in REGIONS)
    and f"{c.rsplit('_', 1)[0]}_whole" in topo.columns
    and c.rsplit("_", 1)[0] not in ("n",)
})

# ── per-tool Spearman r: region length vs region score ────────────────────────
corr_rows = []
for tool in TOOLS:
    row = {"Tool": tool}
    for region, lcol in zip(REGIONS, LEN_COLS):
        rcol = f"{tool}_{region}"
        if rcol not in topo.columns:
            row[region] = np.nan
            continue
        sub = topo[[lcol, rcol]].dropna()
        if len(sub) < 5:
            row[region] = np.nan
        else:
            row[region] = spearmanr(sub[lcol], sub[rcol]).statistic
    corr_rows.append(row)
corr_df = pd.DataFrame(corr_rows).set_index("Tool")

corr_df = corr_df.sort_values("Cytoplasmic", ascending=False)
TOOLS_SORTED = list(corr_df.index)

# ── figure ────────────────────────────────────────────────────────────────────
n_tools = len(TOOLS_SORTED)
n_reg   = len(REGIONS)

colours = ["#e6550d", "#756bb1", "#31a354"]

panel_a_height = n_tools * 1.6
panel_b_height = 5.0
panel_c_height = n_tools * 1.1
fig = plt.figure(figsize=(14, panel_a_height + panel_b_height + panel_c_height + 2.0))
gs_outer = gridspec.GridSpec(3, 1, figure=fig,
                              height_ratios=[panel_a_height, panel_b_height, panel_c_height],
                              hspace=0.35)

# ── Panel A: length vs score, faceted tool × region ──────────────────────────
gs_top = gs_outer[0].subgridspec(n_tools, n_reg, wspace=0.3, hspace=0.15)

for row_i, tool in enumerate(TOOLS_SORTED):
    ra = REGION_APPROACH.get(tool, "concatenated")
    for col_i, (region, lcol, label, col) in enumerate(
            zip(REGIONS, LEN_COLS, REG_LABELS, colours)):
        ax = fig.add_subplot(gs_top[row_i, col_i])
        rcol = f"{tool}_{region}"
        if rcol not in topo.columns:
            ax.set_visible(False)
            continue

        sub = topo[[lcol, rcol]].dropna()
        x, y = sub[lcol], sub[rcol]

        ax.scatter(x, y, color=col, alpha=0.5, s=12, edgecolors="none")

        if len(x) >= 2:
            m, b = np.polyfit(x, y, 1)
            xr = np.linspace(x.min(), x.max(), 100)
            ax.plot(xr, m * xr + b, color=col, lw=1.2, ls="--")

        r, p = spearmanr(x, y)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        ax.text(0.97, 0.05, f"ρ={r:.2f}{sig}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=6,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7))

        ax.tick_params(labelsize=6)
        if row_i == n_tools - 1:
            ax.set_xlabel(label.replace("\n", " ") + "\nlength (res)", fontsize=6)
        else:
            ax.set_xticklabels([])
        if col_i == 0:
            tool_label = f"{tool}\n({ra})" if ra == "pdb_region" else tool
            ax.set_ylabel(tool_label, fontsize=7, rotation=0, ha="right", va="center",
                          labelpad=50)
        if row_i == 0:
            ax.set_title(label, fontsize=8, color=col)

fig.text(0.5, 0.995, "A   Region length vs region score, per tool",
         ha="center", va="top", fontsize=11, fontweight="bold",
         transform=fig.transFigure)

# ── Panel B: heatmap of Spearman r (tool × region) ───────────────────────────
ax_heat = fig.add_subplot(gs_outer[1])

cmap_b = plt.cm.RdBu_r
norm   = mcolors.TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1)

mat = corr_df[REGIONS].values

im = ax_heat.imshow(mat.T, aspect="auto", cmap=cmap_b, norm=norm,
                    interpolation="none")

ax_heat.set_xticks(range(n_tools))
ax_heat.set_xticklabels(
    [f"{t}\n(pdb)" if REGION_APPROACH.get(t) == "pdb_region" else t
     for t in corr_df.index],
    rotation=35, ha="right", fontsize=7
)
ax_heat.set_yticks(range(n_reg))
ax_heat.set_yticklabels(REG_LABELS, fontsize=9)

for i in range(n_reg):
    for j in range(n_tools):
        v = mat[j, i]
        if not np.isnan(v):
            ax_heat.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=6.5,
                        color="white" if abs(v) > 0.5 else "black")

plt.colorbar(im, ax=ax_heat, label="Spearman ρ (length vs score)",
             shrink=0.6, pad=0.01)

ax_heat.set_title(
    "B   Length bias per tool: Spearman ρ between region length and region score\n"
    "    (positive = longer regions score higher, independent of sequence content)",
    fontsize=10, loc="left")

# ── Panel C: non-TMD region length vs whole-protein score ────────────────────
NON_TMD_REGIONS = ["Cytoplasmic", "Extracellular/Lumenal"]
NON_TMD_LCOLS   = ["n_Cytoplasmic", "n_Extracellular/Lumenal"]
NON_TMD_LABELS  = ["Cytoplasmic", "Extracellular /\nLumenal"]
NON_TMD_COLOURS = ["#e6550d", "#31a354"]

gs_bot = gs_outer[2].subgridspec(n_tools, 2, wspace=0.3, hspace=0.15)

for row_i, tool in enumerate(TOOLS_SORTED):
    wcol = f"{tool}_whole"
    if wcol not in topo.columns:
        for col_i in range(2):
            fig.add_subplot(gs_bot[row_i, col_i]).set_visible(False)
        continue
    for col_i, (lcol, label, col) in enumerate(
            zip(NON_TMD_LCOLS, NON_TMD_LABELS, NON_TMD_COLOURS)):
        ax = fig.add_subplot(gs_bot[row_i, col_i])
        sub = topo[[lcol, wcol]].dropna()
        x, y = sub[lcol], sub[wcol]

        ax.scatter(x, y, color=col, alpha=0.5, s=12, edgecolors="none")

        if len(x) >= 2:
            m, b = np.polyfit(x, y, 1)
            xr = np.linspace(x.min(), x.max(), 100)
            ax.plot(xr, m * xr + b, color=col, lw=1.2, ls="--")

        r, p = spearmanr(x, y)
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        ax.text(0.97, 0.05, f"ρ={r:.2f}{sig}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=6,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7))

        ax.tick_params(labelsize=6)
        if row_i == n_tools - 1:
            ax.set_xlabel(label.replace("\n", " ") + "\nlength (res)", fontsize=6)
        else:
            ax.set_xticklabels([])
        if col_i == 0:
            ax.set_ylabel(tool, fontsize=7, rotation=0, ha="right", va="center",
                          labelpad=40)
        if row_i == 0:
            ax.set_title(label, fontsize=8, color=col)

fig.text(0.5, (panel_b_height + panel_c_height) / (panel_a_height + panel_b_height + panel_c_height + 2.0) - 0.005,
         "C   Non-TMD region length vs whole-protein score, per tool\n"
         "    (does longer cytoplasmic / extracellular sequence drive the overall prediction?)",
         ha="center", va="top", fontsize=11, fontweight="bold",
         transform=fig.transFigure)

plt.savefig("output/figures/topology_length_bias.png", dpi=150, bbox_inches="tight")
print("Saved output/figures/topology_length_bias.png")
