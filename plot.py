# %% Imports and config
from ast import literal_eval
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

sns.set_style("whitegrid")
plt.rcParams.update({"font.size": 11, "axes.labelsize": 12, "xtick.labelsize": 9, "ytick.labelsize": 9})

ROOT    = Path(__file__).parent
DATA    = ROOT / "output/full_dataset.csv"
FIG_DIR = ROOT / "output/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

CUTOFF = 0.6   # pLLPS threshold for enrichment plots
TOP_N  = 10    # categories shown in violin/box plots

C = {
    "Membrane":     "#0072B2",
    "Non-Membrane": "#E69F00",
    "Cytosolic":    "#009E73",
    "High":         "#e74c3c",
    "Medium":       "#f39c12",
    "Low":          "#3498db",
}


# %% Load data
df = pd.read_csv(DATA)

df["TMD_count"]   = pd.to_numeric(df["TMD_count"],  errors="coerce").fillna(0).astype(int)
df["Is_Membrane"] = df["TMD_count"] > 0
df["Length"]      = pd.to_numeric(df["Length"], errors="coerce")


def _as_list(value):
    if isinstance(value, list):
        return value
    if isinstance(value, str):
        text = value.strip()
        if text.startswith("["):
            try:
                parsed = literal_eval(text)
                return parsed if isinstance(parsed, list) else []
            except (ValueError, SyntaxError):
                return []
    return []


df["Function_Slim"]      = df["Function_Slim"].apply(_as_list)
df["Location_Slim"]      = df["Location_Slim"].apply(_as_list)
df["GO_Slim_Categories"] = df["GO_Slim_Categories"].apply(_as_list)

df["Function_Slim_Plot"] = df["Function_Slim"].apply(lambda x: x if x else ["Unannotated"])
df["Location_Slim_Plot"] = df["Location_Slim"].apply(lambda x: x if x else ["Unannotated"])

df["GO_BP_list"] = df["GO_BP"].apply(_as_list)
df["GO_MF_list"] = df["GO_MF"].apply(_as_list)
df["GO_CC_list"] = df["GO_CC"].apply(_as_list)

print(f"{len(df):,} proteins  |  {df['pLLPS_Class'].value_counts().to_dict()}")
print(f"Membrane (TMD>0): {df['Is_Membrane'].sum():,}")
print(
    f"GO slim function terms: {df['Function_Slim'].explode().nunique():,} unique  |  "
    f"GO slim location terms: {df['Location_Slim'].explode().nunique():,} unique"
)
df.head(2)


# %% Fig 1 -- pLLPS distribution: count histogram + KDE density by class
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for cls in ["High", "Medium", "Low"]:
    s = df[df["pLLPS_Class"] == cls]["p(LLPS)"]
    axes[0].hist(s, bins=40, alpha=0.72, color=C[cls], label=f"{cls}  (n={len(s):,})")
axes[0].axvline(CUTOFF, color="#333333", linestyle="--", linewidth=1, alpha=0.6, label=f"cutoff {CUTOFF}")
axes[0].set(xlabel="p(LLPS)", ylabel="Count", title="pLLPS distribution by class")
axes[0].legend(title="Class")

for cls in ["High", "Medium", "Low"]:
    s = df[df["pLLPS_Class"] == cls]["p(LLPS)"]
    axes[1].hist(s, bins=40, density=True, alpha=0.28, color=C[cls])
    sns.kdeplot(s, ax=axes[1], color=C[cls], linewidth=2, label=cls)
axes[1].axvline(CUTOFF, color="#333333", linestyle="--", linewidth=1, alpha=0.6)
axes[1].set(xlabel="p(LLPS)", ylabel="Density", title="pLLPS density by class")
axes[1].legend(title="Class")

plt.tight_layout()
plt.savefig(FIG_DIR / "fig1_pllps_distribution.png", dpi=300)
plt.show()


# %% Fig 2 -- Membrane vs non-membrane
groups = {
    "Membrane":     df[df["Is_Membrane"]]["p(LLPS)"].dropna(),
    "Non-Membrane": df[~df["Is_Membrane"]]["p(LLPS)"].dropna(),
}
long   = pd.concat([s.rename("pLLPS").to_frame().assign(Group=k) for k, s in groups.items()])
_, p_mw = stats.mannwhitneyu(*groups.values(), alternative="two-sided")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for label, s in groups.items():
    axes[0].hist(s, bins=40, density=True, alpha=0.45, color=C[label],
                 label=f"{label}  (n={len(s):,})")
    sns.kdeplot(s, ax=axes[0], color=C[label], linewidth=2)
axes[0].axvline(CUTOFF, color="#333333", linestyle="--", linewidth=1, alpha=0.6)
axes[0].set(xlabel="p(LLPS)", ylabel="Density", title="Histogram + KDE overlay")
axes[0].legend()

palette = {g: C[g] for g in groups}
sns.violinplot(data=long, x="Group", y="pLLPS", inner=None, cut=0,
               hue="Group", palette=palette, legend=False, ax=axes[1], alpha=0.6)
sns.boxplot(data=long, x="Group", y="pLLPS", width=0.18, showfliers=False,
            hue="Group", palette=palette, legend=False, ax=axes[1],
            boxprops={"facecolor": "white", "edgecolor": "#333333", "linewidth": 0.9},
            whiskerprops={"color": "#555555"}, capprops={"color": "#555555"},
            medianprops={"color": "#111111", "linewidth": 1.5})
axes[1].set(ylim=(0, 1.1), xlabel="", ylabel="p(LLPS)",
            title=f"Violin + box  (Mann-Whitney p = {p_mw:.2e})")

plt.suptitle("Membrane vs Non-Membrane pLLPS", fontsize=13)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig2_membrane_vs_nonmembrane.png", dpi=300)
plt.show()


# %% Fig 3 -- pLLPS histograms by top GO slim function categories
func_long = (
    df[["p(LLPS)", "Function_Slim_Plot"]]
    .explode("Function_Slim_Plot")
    .rename(columns={"Function_Slim_Plot": "Function"})
    .query("Function not in ['Unannotated', 'Other', '']")
    .dropna(subset=["p(LLPS)"])
)
top_funcs = func_long["Function"].value_counts().nlargest(6).index.tolist()
s_all     = df["p(LLPS)"].dropna()

fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True)
axes = axes.flatten()
for i, func in enumerate(top_funcs):
    ax = axes[i]
    s  = func_long.loc[func_long["Function"] == func, "p(LLPS)"]
    ax.hist(s_all, bins=40, density=True, alpha=0.15, color="#999999", label="All proteins")
    ax.hist(s,     bins=40, density=True, alpha=0.50, color=C["Membrane"], label=func)
    sns.kdeplot(s, ax=ax, color=C["Membrane"], linewidth=2)
    ax.axvline(CUTOFF, color="#e74c3c", linestyle="--", linewidth=0.9, alpha=0.7)
    ax.set(title=f"{func}  (n={len(s):,})", xlabel="p(LLPS)", ylabel="Density")
    ax.legend(fontsize=8)

plt.suptitle("pLLPS by top GO slim function  (grey = all proteins)", fontsize=13)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig3_pllps_by_go_slim_function.png", dpi=300, bbox_inches="tight")
plt.show()


# %% Fig 4 -- pLLPS by GO slim location (all proteins, top 10 by count)
loc_long = (
    df[["p(LLPS)", "Location_Slim_Plot"]]
    .explode("Location_Slim_Plot")
    .dropna(subset=["Location_Slim_Plot"])
    .query("Location_Slim_Plot not in ['Other', '', 'Unannotated']")
    .rename(columns={"Location_Slim_Plot": "Location"})
    .dropna(subset=["p(LLPS)"])
)
top_locs  = loc_long["Location"].value_counts().nlargest(TOP_N).index
loc_long  = loc_long[loc_long["Location"].isin(top_locs)]
order_locs = (
    loc_long.groupby("Location")["p(LLPS)"]
    .median().sort_values(ascending=False).index.tolist()
)

fig, ax = plt.subplots(figsize=(12, max(6, 0.42 * len(order_locs) + 1.5)))
sns.violinplot(data=loc_long, y="Location", x="p(LLPS)", order=order_locs,
               inner=None, cut=0, color=C["Cytosolic"], alpha=0.55, ax=ax)
sns.boxplot(data=loc_long, y="Location", x="p(LLPS)", order=order_locs,
            width=0.15, showfliers=False, ax=ax,
            boxprops={"facecolor": "white", "edgecolor": "#333333", "linewidth": 0.9},
            whiskerprops={"color": "#555555"}, capprops={"color": "#555555"},
            medianprops={"color": "#111111", "linewidth": 1.2})
ax.axvline(CUTOFF, color="#e74c3c", linestyle="--", linewidth=0.9, alpha=0.6)
ax.set(xlim=(0, 1.1), xlabel="p(LLPS)", ylabel="",
       title=f"pLLPS by GO slim location -- all proteins  (top {TOP_N} by count)")
ax.tick_params(axis="y", labelsize=10)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig4_by_location.png", dpi=300, bbox_inches="tight")
plt.show()


# %% Fig 5 -- 2D KDE: pLLPS vs length, DPR count, TMD count
rng = np.random.default_rng(42)

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

s_len   = df.dropna(subset=["p(LLPS)", "Length"])
log_len = np.log10(s_len["Length"].clip(lower=1))
sns.kdeplot(x=log_len, y=s_len["p(LLPS)"], ax=axes[0],
            fill=True, cmap="Blues", thresh=0.04, levels=14)
axes[0].scatter(log_len, s_len["p(LLPS)"],
                alpha=0.015, s=4, color="#333333", rasterized=True)
axes[0].axhline(CUTOFF, color="#e74c3c", linestyle="--", linewidth=1, alpha=0.7)
axes[0].set(xlabel="log₁₀(Length, aa)", ylabel="p(LLPS)",
            title=f"pLLPS vs protein length  (n={len(s_len):,})")

s_dpr = df.dropna(subset=["p(LLPS)", "n(DPR=> 25)"])
dpr   = s_dpr["n(DPR=> 25)"].clip(upper=40)
sns.kdeplot(x=dpr, y=s_dpr["p(LLPS)"], ax=axes[1],
            fill=True, cmap="Greens", thresh=0.04, levels=14)
axes[1].scatter(dpr, s_dpr["p(LLPS)"],
                alpha=0.015, s=4, color="#333333", rasterized=True)
axes[1].axhline(CUTOFF, color="#e74c3c", linestyle="--", linewidth=1, alpha=0.7)
axes[1].set(xlabel="n(DPR >= 25)  [clipped at 40]", ylabel="p(LLPS)",
            title=f"pLLPS vs disordered region count  (n={len(s_dpr):,})")

tmd_grp   = df["TMD_count"].clip(upper=8).astype(str)
tmd_order = [str(i) for i in sorted(df["TMD_count"].clip(upper=8).unique())]
jitter    = rng.uniform(-0.25, 0.25, len(df))
axes[2].scatter(df["TMD_count"].clip(upper=8) + jitter, df["p(LLPS)"],
                alpha=0.03, s=4, color="#888888", rasterized=True)
sns.violinplot(data=df.assign(_tmd=tmd_grp), x="_tmd", y="p(LLPS)",
               order=tmd_order, inner=None, cut=0,
               color=C["Cytosolic"], alpha=0.5, ax=axes[2])
sns.boxplot(data=df.assign(_tmd=tmd_grp), x="_tmd", y="p(LLPS)",
            order=tmd_order, width=0.18, showfliers=False, ax=axes[2],
            boxprops={"facecolor": "white", "edgecolor": "#333333", "linewidth": 0.9},
            whiskerprops={"color": "#555555"}, capprops={"color": "#555555"},
            medianprops={"color": "#111111", "linewidth": 1.5})
axes[2].axhline(CUTOFF, color="#e74c3c", linestyle="--", linewidth=1, alpha=0.7)
axes[2].set(xlabel="TMD count  (8 = 8+)", ylabel="p(LLPS)",
            title="pLLPS by TMD count  (violin + box)")

plt.suptitle("pLLPS vs structural features", fontsize=13)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig5_2d_kde.png", dpi=300, bbox_inches="tight")
plt.show()


# %% Fig 6 -- Membrane proteins: pLLPS by GO slim function (top 10)
mem       = df[df["Is_Membrane"]].copy()
func_long = (
    mem[["p(LLPS)", "Function_Slim_Plot"]]
    .explode("Function_Slim_Plot")
    .dropna(subset=["Function_Slim_Plot"])
    .query("Function_Slim_Plot not in ['Other', '', 'Unannotated']")
    .rename(columns={"Function_Slim_Plot": "Function"})
    .dropna(subset=["p(LLPS)"])
)
func_counts = func_long["Function"].value_counts()
top_funcs   = func_counts.nlargest(TOP_N).index
func_long   = func_long[func_long["Function"].isin(top_funcs)].copy()
order_func  = func_counts.loc[top_funcs].sort_values(ascending=False).index.tolist()

fig, ax = plt.subplots(figsize=(12, max(7, 0.35 * len(order_func) + 1.5)))
sns.violinplot(data=func_long, y="Function", x="p(LLPS)", order=order_func,
               inner=None, cut=0, color=C["Membrane"], alpha=0.55, ax=ax)
sns.boxplot(data=func_long, y="Function", x="p(LLPS)", order=order_func,
            width=0.15, showfliers=False, ax=ax,
            boxprops={"facecolor": "white", "edgecolor": "#333333", "linewidth": 0.9},
            whiskerprops={"color": "#555555"}, capprops={"color": "#555555"},
            medianprops={"color": "#111111", "linewidth": 1.2})
ax.axvline(CUTOFF, color="#e74c3c", linestyle="--", linewidth=0.9, alpha=0.6)
ax.set(xlim=(0, 1.1), xlabel="p(LLPS)", ylabel="",
       title=f"Membrane proteins -- pLLPS by GO slim function  (top {TOP_N},  n={len(func_long):,})")
ax.tick_params(axis="y", labelsize=9)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig6_membrane_by_function.png", dpi=300, bbox_inches="tight")
plt.show()


# %% Fig 7 -- Membrane proteins: pLLPS by GO slim location (top 10, excl. cytosol)
loc_long = (
    mem[["p(LLPS)", "Location_Slim_Plot"]]
    .explode("Location_Slim_Plot")
    .dropna(subset=["Location_Slim_Plot"])
    .query("Location_Slim_Plot not in ['Other', '', 'Unannotated', 'cytosol']")
    .rename(columns={"Location_Slim_Plot": "Location"})
    .dropna(subset=["p(LLPS)"])
)
loc_counts = loc_long["Location"].value_counts()
top_locs   = loc_counts.nlargest(TOP_N).index
loc_long   = loc_long[loc_long["Location"].isin(top_locs)].copy()
order_locs = loc_counts.loc[top_locs].sort_values(ascending=False).index.tolist()

fig, ax = plt.subplots(figsize=(12, max(7, 0.42 * len(order_locs) + 1.5)))
sns.violinplot(data=loc_long, y="Location", x="p(LLPS)", order=order_locs,
               inner=None, cut=0, color=C["Cytosolic"], alpha=0.55, ax=ax)
sns.boxplot(data=loc_long, y="Location", x="p(LLPS)", order=order_locs,
            width=0.15, showfliers=False, ax=ax,
            boxprops={"facecolor": "white", "edgecolor": "#333333", "linewidth": 0.9},
            whiskerprops={"color": "#555555"}, capprops={"color": "#555555"},
            medianprops={"color": "#111111", "linewidth": 1.2})
ax.axvline(CUTOFF, color="#e74c3c", linestyle="--", linewidth=0.9, alpha=0.6)
ax.set(xlim=(0, 1.1), xlabel="p(LLPS)", ylabel="",
       title=f"Membrane proteins -- pLLPS by GO slim location  (top {TOP_N}, excl. cytosol,  n={len(loc_long):,})")
ax.tick_params(axis="y", labelsize=10)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig7_membrane_by_location.png", dpi=300, bbox_inches="tight")
plt.show()


# %% Fig 8 -- % above cutoff by GO slim function (membrane, n > 10)
overall_pct = 100 * (mem["p(LLPS)"] > CUTOFF).sum() / len(mem)

func_stats = (
    mem[["Function_Slim_Plot", "p(LLPS)"]]
    .explode("Function_Slim_Plot")
    .dropna(subset=["Function_Slim_Plot"])
    .query("Function_Slim_Plot not in ['Other', '', 'Unannotated']")
    .groupby("Function_Slim_Plot")["p(LLPS)"]
    .agg(n="count", n_above=lambda s: (s > CUTOFF).sum())
    .reset_index()
)
func_stats["pct"]   = 100 * func_stats["n_above"] / func_stats["n"]
func_stats_filt     = func_stats[func_stats["n"] > 10].sort_values("pct", ascending=True)

fig, ax = plt.subplots(figsize=(12, max(7, 0.3 * len(func_stats_filt) + 1.5)))
ax.barh(func_stats_filt["Function_Slim_Plot"], func_stats_filt["pct"],
        color=C["Membrane"], alpha=0.85, edgecolor="white")
ax.axvline(overall_pct, color="#e74c3c", linestyle="--", linewidth=1.5,
           label=f"Membrane overall  {overall_pct:.1f}%")
for y, (_, row) in enumerate(func_stats_filt.iterrows()):
    ax.text(row["pct"] + 0.8, y, f"n={int(row['n'])}", va="center", fontsize=8, color="#333333")
ax.set(
    xlim=(0, max(100, func_stats_filt["pct"].max() * 1.12)),
    xlabel=f"% with p(LLPS) > {CUTOFF}", ylabel="",
    title=f"% above p(LLPS) {CUTOFF} by GO slim function  (membrane, n > 10)",
)
ax.tick_params(axis="y", labelsize=9)
ax.legend()
plt.tight_layout()
plt.savefig(FIG_DIR / "fig8_pct_above_cutoff_function.png", dpi=300, bbox_inches="tight")
plt.show()


# %% Fig 9 -- % above cutoff by GO slim location (membrane, n > 10, excl. cytosol)
loc_stats = (
    mem[["Location_Slim_Plot", "p(LLPS)"]]
    .explode("Location_Slim_Plot")
    .dropna(subset=["Location_Slim_Plot"])
    .query("Location_Slim_Plot not in ['Other', '', 'Unannotated', 'cytosol']")
    .groupby("Location_Slim_Plot")["p(LLPS)"]
    .agg(n="count", n_above=lambda s: (s > CUTOFF).sum())
    .reset_index()
)
loc_stats["pct"]  = 100 * loc_stats["n_above"] / loc_stats["n"]
loc_stats_filt    = loc_stats[loc_stats["n"] > 10].sort_values("pct", ascending=True)

fig, ax = plt.subplots(figsize=(12, max(7, 0.42 * len(loc_stats_filt) + 1.5)))
ax.barh(loc_stats_filt["Location_Slim_Plot"], loc_stats_filt["pct"],
        color=C["Cytosolic"], alpha=0.85, edgecolor="white")
ax.axvline(overall_pct, color="#e74c3c", linestyle="--", linewidth=1.5,
           label=f"Membrane overall  {overall_pct:.1f}%")
for y, (_, row) in enumerate(loc_stats_filt.iterrows()):
    ax.text(row["pct"] + 0.8, y, f"n={int(row['n'])}", va="center", fontsize=8, color="#333333")
ax.set(
    xlim=(0, max(100, loc_stats_filt["pct"].max() * 1.12)),
    xlabel=f"% with p(LLPS) > {CUTOFF}", ylabel="",
    title=f"% above p(LLPS) {CUTOFF} by GO slim location  (membrane, n > 10, excl. cytosol)",
)
ax.tick_params(axis="y", labelsize=9)
ax.legend()
plt.tight_layout()
plt.savefig(FIG_DIR / "fig9_pct_above_cutoff_location.png", dpi=300, bbox_inches="tight")
plt.show()


# %% GO term filter -- ad-hoc exploration
# Edit FILTER_GO_TERMS and re-run this cell to compare any subset vs the rest.
FILTER_GO_TERMS: list[str] = [
    # "GO:0005886",  # plasma membrane (CC)
    # "GO:0003723",  # RNA binding (MF)
    # "GO:0000785",  # chromatin (CC)
]
FILTER_LABEL = "GO term subset"

if not FILTER_GO_TERMS:
    print("No GO terms specified — edit FILTER_GO_TERMS above and re-run.\n")
    print("Top GO IDs in dataset (by frequency):")
    for col, label in [("GO_BP_list", "BP"), ("GO_MF_list", "MF"), ("GO_CC_list", "CC")]:
        sample = df[col].explode().dropna().value_counts().head(5)
        print(f"  {label}: {sample.index.tolist()}")
else:
    def _has_go_term(row, terms):
        all_ids = row["GO_BP_list"] + row["GO_MF_list"] + row["GO_CC_list"]
        return any(t in all_ids for t in terms)

    mask   = df.apply(_has_go_term, axis=1, args=(FILTER_GO_TERMS,))
    subset = df[mask].copy()
    rest   = df[~mask].copy()
    print(f"Terms: {FILTER_GO_TERMS}")
    print(f"Matching: {mask.sum():,} / {len(df):,}  ({100 * mask.mean():.1f}%)")

    _, p_mw = stats.mannwhitneyu(
        subset["p(LLPS)"].dropna(), rest["p(LLPS)"].dropna(), alternative="two-sided"
    )

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for s, color, label in [
        (rest["p(LLPS)"].dropna(),   "#999999",     f"Other  (n={len(rest):,})"),
        (subset["p(LLPS)"].dropna(), C["Membrane"], f"{FILTER_LABEL}  (n={len(subset):,})"),
    ]:
        axes[0].hist(s, bins=40, density=True, alpha=0.4, color=color, label=label)
        sns.kdeplot(s, ax=axes[0], color=color, linewidth=2)
    axes[0].axvline(CUTOFF, color="#e74c3c", linestyle="--", linewidth=1, alpha=0.7)
    axes[0].set(xlabel="p(LLPS)", ylabel="Density",
                title=f"{FILTER_LABEL} vs rest  (Mann-Whitney p={p_mw:.2e})")
    axes[0].legend()

    long_go = pd.concat([
        subset["p(LLPS)"].dropna().rename("pLLPS").to_frame().assign(Group=FILTER_LABEL),
        rest["p(LLPS)"].dropna().rename("pLLPS").to_frame().assign(Group="Other"),
    ])
    pal = {FILTER_LABEL: C["Membrane"], "Other": "#999999"}
    sns.violinplot(data=long_go, x="Group", y="pLLPS", inner=None, cut=0,
                   hue="Group", palette=pal, legend=False, ax=axes[1], alpha=0.6)
    sns.boxplot(data=long_go, x="Group", y="pLLPS", width=0.18, showfliers=False,
                hue="Group", palette=pal, legend=False, ax=axes[1],
                boxprops={"facecolor": "white", "edgecolor": "#333333", "linewidth": 0.9},
                whiskerprops={"color": "#555555"}, capprops={"color": "#555555"},
                medianprops={"color": "#111111", "linewidth": 1.5})
    axes[1].set(ylim=(0, 1.1), xlabel="", ylabel="p(LLPS)", title="Violin + box")
    plt.tight_layout()
    plt.show()


# %% Fig 10 -- pLLPS distribution: experimental DB proteins vs background
C_EXP  = "#8E44AD"
C_BACK = "#7F8C8D"

if "in_exp_db" not in df.columns:
    print("Skipping Fig 10/11: run wrangle_llps_dbs.py first to add in_exp_db column.")
else:
    exp = df[df["in_exp_db"]]["p(LLPS)"].dropna()
    bg  = df[~df["in_exp_db"]]["p(LLPS)"].dropna()
    _, p_mw_exp = stats.mannwhitneyu(exp, bg, alternative="two-sided")

    long_exp = pd.concat([
        exp.rename("pLLPS").to_frame().assign(Group=f"Experimental DB  (n={len(exp):,})"),
        bg.rename("pLLPS").to_frame().assign(Group=f"Background  (n={len(bg):,})"),
    ])
    pal_exp = {f"Experimental DB  (n={len(exp):,})": C_EXP, f"Background  (n={len(bg):,})": C_BACK}

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].hist(bg,  bins=40, density=True, alpha=0.40, color=C_BACK, label=f"Background  (n={len(bg):,})")
    axes[0].hist(exp, bins=40, density=True, alpha=0.55, color=C_EXP,  label=f"Experimental DB  (n={len(exp):,})")
    sns.kdeplot(bg,  ax=axes[0], color=C_BACK, linewidth=2)
    sns.kdeplot(exp, ax=axes[0], color=C_EXP,  linewidth=2)
    axes[0].axvline(CUTOFF, color="#333333", linestyle="--", linewidth=1, alpha=0.6, label=f"cutoff {CUTOFF}")
    axes[0].set(xlabel="p(LLPS)", ylabel="Density", title="pLLPS: experimental DB vs background")
    axes[0].legend()

    sns.violinplot(data=long_exp, x="Group", y="pLLPS", inner=None, cut=0,
                   hue="Group", palette=pal_exp, legend=False, ax=axes[1], alpha=0.6)
    sns.boxplot(data=long_exp, x="Group", y="pLLPS", width=0.18, showfliers=False,
                hue="Group", palette=pal_exp, legend=False, ax=axes[1],
                boxprops={"facecolor": "white", "edgecolor": "#333333", "linewidth": 0.9},
                whiskerprops={"color": "#555555"}, capprops={"color": "#555555"},
                medianprops={"color": "#111111", "linewidth": 1.5})
    axes[1].set(ylim=(0, 1.1), xlabel="", ylabel="p(LLPS)",
                title=f"Violin + box  (Mann-Whitney p = {p_mw_exp:.2e})")

    plt.suptitle("pLLPS scores: experimental LLPS DB proteins vs rest", fontsize=13)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig10_exp_db_distribution.png", dpi=300)
    plt.show()


# %% Fig 11 -- ROC curve + pLLPS class enrichment for experimental DB
if "in_exp_db" in df.columns:
    valid    = df["p(LLPS)"].notna() & df["in_exp_db"].notna()
    y_true   = df.loc[valid, "in_exp_db"].astype(int).values
    y_score  = df.loc[valid, "p(LLPS)"].values
    n_pos, n_neg = y_true.sum(), (len(y_true) - y_true.sum())

    thresholds = np.sort(np.unique(y_score))[::-1]
    tpr_pts, fpr_pts = [0.0], [0.0]
    for t in thresholds:
        pred = y_score >= t
        tpr_pts.append((pred & y_true.astype(bool)).sum() / n_pos)
        fpr_pts.append((pred & ~y_true.astype(bool)).sum() / n_neg)
    tpr_pts.append(1.0)
    fpr_pts.append(1.0)
    auc = float(np.trapezoid(tpr_pts, fpr_pts))

    class_stats = (
        df.groupby("pLLPS_Class")["in_exp_db"]
        .agg(n="count", n_exp=lambda s: s.sum())
        .reindex(["Low", "Medium", "High"])
        .reset_index()
    )
    class_stats["pct_exp"] = 100 * class_stats["n_exp"] / class_stats["n"]
    overall_pct_exp = 100 * df["in_exp_db"].sum() / len(df)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].plot(fpr_pts, tpr_pts, color=C_EXP, linewidth=2, label=f"pLLPS  (AUC = {auc:.3f})")
    axes[0].plot([0, 1], [0, 1], color="#aaaaaa", linestyle="--", linewidth=1, label="Random")
    axes[0].fill_between(fpr_pts, tpr_pts, alpha=0.10, color=C_EXP)
    axes[0].set(xlabel="False Positive Rate", ylabel="True Positive Rate",
                title="ROC: pLLPS as predictor of experimental LLPS")
    axes[0].legend()

    bar_colors = [C["Low"], C["Medium"], C["High"]]
    axes[1].bar(class_stats["pLLPS_Class"], class_stats["pct_exp"],
                color=bar_colors, alpha=0.85, edgecolor="white")
    axes[1].axhline(overall_pct_exp, color="#333333", linestyle="--", linewidth=1.5,
                    label=f"Overall  {overall_pct_exp:.1f}%")
    for i, row in enumerate(class_stats.itertuples()):
        axes[1].text(i, row.pct_exp + 0.3, f"n={int(row.n_exp)}", ha="center", fontsize=9, color="#333333")
    axes[1].set(xlabel="pLLPS class", ylabel="% in experimental DB",
                title="Enrichment of experimental LLPS proteins by pLLPS class")
    axes[1].legend()

    plt.suptitle("pLLPS as a predictor of experimental LLPS", fontsize=13)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig11_exp_db_roc_enrichment.png", dpi=300)
    plt.show()


# %% Fig 12 -- in_exp_db class enrichment: membrane vs non-membrane
if "in_exp_db" in df.columns:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)

    for ax, mask, label, color in [
        (axes[0], df["Is_Membrane"],  "Membrane",     C["Membrane"]),
        (axes[1], ~df["Is_Membrane"], "Non-Membrane", C["Non-Membrane"]),
    ]:
        sub = df[mask]
        stats_ = (
            sub.groupby("pLLPS_Class")["in_exp_db"]
            .agg(n="count", n_exp="sum")
            .reindex(["Low", "Medium", "High"])
            .reset_index()
        )
        stats_["pct"] = 100 * stats_["n_exp"] / stats_["n"]
        overall = 100 * sub["in_exp_db"].sum() / len(sub)

        ax.bar(stats_["pLLPS_Class"], stats_["pct"],
               color=[C["Low"], C["Medium"], C["High"]], alpha=0.85, edgecolor="white")
        ax.axhline(overall, color="#333333", linestyle="--", linewidth=1.5,
                   label=f"Overall  {overall:.1f}%")
        for i, row in enumerate(stats_.itertuples()):
            ax.text(i, row.pct + 0.05, f"n={int(row.n_exp)}", ha="center", fontsize=9, color="#333333")
        ax.set(xlabel="pLLPS class", ylabel="% in experimental DB",
               title=f"{label}  (n={len(sub):,}, {sub['in_exp_db'].sum()} in exp DBs)")
        ax.legend()

    plt.suptitle("Experimental DB enrichment by pLLPS class: membrane vs non-membrane", fontsize=13)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig12_exp_db_membrane_vs_nonmembrane.png", dpi=300)
    plt.show()


# %% Summary
print("Figures written to", FIG_DIR.resolve())
for f in sorted(FIG_DIR.glob("*.png")):
    print(" ", f.name)
