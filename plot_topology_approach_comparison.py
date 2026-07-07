"""
Approach comparison: concatenated sequence vs AF2 PDB region slice.

For each topology region, shows the distribution of delta scores
(region - whole) across all proteins, split by approach family and tool.

Layout:
  Row = region (Cytoplasmic / Transmembrane / Extracellular/Lumenal)
  Within each row: strip + box plot per tool, coloured by approach
    - blue: concatenated (isolated sequence scored as one)
    - orange: pdb_region (structural AF2 slice; PICNIC, PSPire)
    - grey: region_mean / segment (derived from full-protein per-residue scores)

A separate panel (row 4) directly compares the two structural tools' deltas
across regions to confirm cytoplasmic > extracellular > TM.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

REGIONS       = ["Cytoplasmic", "Transmembrane", "Extracellular/Lumenal"]
REG_LABELS    = ["Cytoplasmic", "Transmembrane", "Extracellular / Lumenal"]

WHOLE_APPROACH = {"ESpritz": "whole_db", "PLAAC": "whole_db", "PScore": "whole_db",
                  "SEG": "whole_db", "catGRANULE": "whole_db", "PSPire": "whole_db"}

APPROACH_COLOUR = {
    "concatenated": "#3182bd",
    "pdb_region":   "#e6550d",
    "region_mean":  "#74c476",
    "segment":      "#9e9e9e",
}
APPROACH_LABEL = {
    "concatenated": "concatenated (isolated sequence)",
    "pdb_region":   "pdb_region (AF2 structural slice)",
    "region_mean":  "region_mean (full-protein per-residue, avg)",
    "segment":      "segment (per-span, max score used)",
}

master = pd.read_csv("output/topology_scores_master.csv")

# ── build delta for each (protein, tool, approach, region) ────────────────────
# Whole-protein score per tool per protein
whole_rows = []
for tool in master.tool.unique():
    df_t = master[master.tool == tool]
    wa = WHOLE_APPROACH.get(tool, "whole")
    whole = df_t[df_t.approach == wa][["UniProt_ID", "score"]].drop_duplicates("UniProt_ID")
    if whole.empty:
        alt = "whole" if wa == "whole_db" else "whole_db"
        whole = df_t[df_t.approach == alt][["UniProt_ID", "score"]].drop_duplicates("UniProt_ID")
    whole["tool"] = tool
    whole_rows.append(whole)

whole_df = pd.concat(whole_rows, ignore_index=True).rename(columns={"score": "whole_score"})

# Region rows: exclude 'whole', 'whole_db', 'segment' for the main comparison
# Use max segment score as an aggregate for the segment approach
region_master = master[master.region.isin(REGIONS)].copy()

# For segment approach: take max score per (protein, tool, region)
seg = (region_master[region_master.approach == "segment"]
       .groupby(["UniProt_ID", "tool", "region"])["score"].max()
       .reset_index())
seg["approach"] = "segment"

non_seg = region_master[region_master.approach != "segment"].copy()

all_region = pd.concat([non_seg[["UniProt_ID","tool","approach","region","score"]], seg],
                       ignore_index=True)

# Merge whole scores and compute delta
delta_df = all_region.merge(whole_df, on=["UniProt_ID", "tool"], how="inner")
delta_df["delta"] = delta_df["score"] - delta_df["whole_score"]

# Drop region_mean — identical to concatenated for 5/6 tools (verified)
# Keep it only where it meaningfully differs (ParSe2 has slight difference)
delta_df = delta_df[delta_df.approach != "region_mean"]

# ── figure layout: 3 region rows + 1 structural-tools summary ─────────────────
fig, axes = plt.subplots(4, 1, figsize=(18, 22), constrained_layout=True)

tools_ordered = sorted(master.tool.unique())

for row_i, region in enumerate(REGIONS):
    ax = axes[row_i]
    sub = delta_df[delta_df.region == region]

    approaches_in_data = sorted(sub.approach.unique())
    offsets = {a: (i - (len(approaches_in_data) - 1) / 2) * 0.28
               for i, a in enumerate(approaches_in_data)}

    for tool_i, tool in enumerate(tools_ordered):
        for approach in approaches_in_data:
            vals = sub[(sub.tool == tool) & (sub.approach == approach)]["delta"].dropna()
            if len(vals) == 0:
                continue

            xc  = tool_i + offsets[approach]
            col = APPROACH_COLOUR.get(approach, "#aaaaaa")

            # jitter strip
            jitter = np.random.default_rng(tool_i + hash(approach) % 1000).uniform(
                -0.08, 0.08, len(vals))
            ax.scatter(xc + jitter, vals, color=col, alpha=0.45, s=16,
                       linewidths=0, zorder=3)

            # box
            bp = ax.boxplot(vals, positions=[xc], widths=0.18, patch_artist=True,
                            manage_ticks=False,
                            boxprops=dict(facecolor=col, alpha=0.3, linewidth=0.8),
                            medianprops=dict(color="black", linewidth=1.5),
                            whiskerprops=dict(linewidth=0.8, color=col),
                            capprops=dict(linewidth=0.8, color=col),
                            flierprops=dict(marker=".", ms=3, alpha=0.4, color=col))

    ax.axhline(0, color="#444444", lw=1.0, ls="--", alpha=0.6, zorder=1)
    ax.set_xticks(range(len(tools_ordered)))
    ax.set_xticklabels(
        [f"{t}\n(pdb)" if t in ("PICNIC", "PSPire") else t
         for t in tools_ordered],
        fontsize=7.5, rotation=30, ha="right"
    )
    ax.set_ylabel("Δ score (region − whole)", fontsize=9)
    ax.set_title(f"{REG_LABELS[row_i]}", fontsize=11, fontweight="bold", pad=4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(-0.7, len(tools_ordered) - 0.3)

# ── Panel 4: structural tools (PICNIC, PSPire) across regions ─────────────────
ax4 = axes[3]
struct_tools = ["PICNIC", "PSPire"]
colours_struct = {"PICNIC": "#1f77b4", "PSPire": "#ff7f0e"}

x_positions = {r: i for i, r in enumerate(REGIONS)}
width = 0.35

for t_i, tool in enumerate(struct_tools):
    sub = delta_df[(delta_df.tool == tool) & (delta_df.approach == "pdb_region")]
    for region in REGIONS:
        vals = sub[sub.region == region]["delta"].dropna()
        if len(vals) == 0:
            continue
        xc = x_positions[region] + (t_i - 0.5) * width
        col = colours_struct[tool]

        jitter = np.random.default_rng(t_i * 10 + x_positions[region]).uniform(
            -0.06, 0.06, len(vals))
        ax4.scatter(xc + jitter, vals, color=col, alpha=0.55, s=20,
                    linewidths=0, zorder=3)
        ax4.boxplot(vals, positions=[xc], widths=0.25, patch_artist=True,
                    manage_ticks=False,
                    boxprops=dict(facecolor=col, alpha=0.3, linewidth=0.8),
                    medianprops=dict(color="black", linewidth=1.8),
                    whiskerprops=dict(linewidth=0.8, color=col),
                    capprops=dict(linewidth=0.8, color=col),
                    flierprops=dict(marker=".", ms=3, alpha=0.4, color=col))

ax4.axhline(0, color="#444444", lw=1.0, ls="--", alpha=0.6, zorder=1)
ax4.set_xticks(range(len(REGIONS)))
ax4.set_xticklabels(REG_LABELS, fontsize=10)
ax4.set_ylabel("Δ score (pdb_region − whole)", fontsize=9)
ax4.set_title("Structure-based tools (pdb_region approach): Δ score by region",
              fontsize=11, fontweight="bold", pad=4)
ax4.spines[["top", "right"]].set_visible(False)
ax4.set_xlim(-0.7, len(REGIONS) - 0.3)
ax4.legend(handles=[mpatches.Patch(color=colours_struct[t], label=t, alpha=0.7)
                    for t in struct_tools],
           fontsize=9, frameon=False, loc="upper right")

# ── shared legend ─────────────────────────────────────────────────────────────
legend_patches = [
    mpatches.Patch(color=APPROACH_COLOUR[a], label=APPROACH_LABEL[a], alpha=0.7)
    for a in ["concatenated", "pdb_region", "segment"]
]
fig.legend(handles=legend_patches, loc="upper center", ncol=3, fontsize=9,
           frameon=False, bbox_to_anchor=(0.5, 1.01))

fig.suptitle(
    "Approach comparison: Δ score (region − whole) by tool and scoring method\n"
    "Row 1–3: all tools, all approaches; Row 4: structure-based tools only (PICNIC & PSPire)",
    fontsize=13, y=1.03
)

plt.savefig("output/figures/topology_approach_comparison.png",
            dpi=150, bbox_inches="tight")
print("Saved output/figures/topology_approach_comparison.png")
