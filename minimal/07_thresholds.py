"""
STEP 7 -- published cutoffs: who does each tool actually call positive?

An AUROC says how well a tool RANKS. It says nothing about what happens when you use
the tool the way a paper tells you to: apply its published cutoff and read off a
yes/no. This script does that, for the 10 tools that have a defensible published
cutoff (config.PUBLISHED_CUTOFF explains why the other 8 do not).

Two results:

  A. Pass rate per tool. These tools are all nominally answering one question, so
     they ought to broadly agree on how many membrane proteins pass. They do not.

  B. The agreement shortlist. Per protein, what fraction of its available cutoffs
     does it clear? Reported at two levels, >=75% (strict) and >=50% (permissive),
     as a pair of bounds rather than one cutoff pretending to be definitive.

The denominator is PER PROTEIN, not a fixed 10: four tools fail to score some
proteins, and a protein is judged only against the cutoffs it could be evaluated on.

Run:  python 07_thresholds.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C

m = pd.read_csv(C.BUILD / "master.csv")

cutoffs = {c: t for c, t in C.PUBLISHED_CUTOFF.items() if c in m.columns}
print(f"{len(cutoffs)} of {len(C.PREDICTORS)} tools have a published cutoff\n")

# --------------------------------------------------------------- A. per-tool pass rate
rows = []
for col, thr in cutoffs.items():
    s = m[col]
    scored = s.notna()
    n_pass = int((s >= thr).sum())          # NaN >= x is False, so this is safe
    rows.append({
        "tool": C.NICE[col], "cutoff": thr,
        "n_scored": int(scored.sum()), "n_pass": n_pass,
        "pass_pct": round(100 * n_pass / scored.sum(), 1),
    })

per_tool = pd.DataFrame(rows).sort_values("pass_pct", ascending=False)

print("A. HOW MANY OF THE SET EACH TOOL CALLS POSITIVE AT ITS OWN CUTOFF")
print("-" * 66)
print(per_tool.to_string(index=False))

spread = per_tool.pass_pct.max() / per_tool.pass_pct.min()
print(f"\n   pass rate ranges {per_tool.pass_pct.min()}% to {per_tool.pass_pct.max()}%"
      f"  -- a {spread:.0f}-fold spread across tools")
print("   answering nominally the same question on the same proteins.")

# ------------------------------------------------------- B. per-protein agreement
# Count, for each protein, how many of ITS evaluable cutoffs it clears.
passes = pd.DataFrame({col: m[col] >= thr for col, thr in cutoffs.items()})
evaluable = pd.DataFrame({col: m[col].notna() for col in cutoffs})

m["n_pass"] = passes.sum(axis=1)
m["n_eval"] = evaluable.sum(axis=1)
m["frac_pass"] = m.n_pass / m.n_eval

print("\n\nB. PER-PROTEIN AGREEMENT ACROSS THE 10 CUTOFFS")
print("-" * 66)
print("   proteins clearing n cutoffs:")
for n, cnt in m.n_pass.value_counts().sort_index().items():
    print(f"     {n:>2d} of ~10   {cnt:>4d}")

strict = m[m.frac_pass >= 0.75]
permissive = m[m.frac_pass >= 0.50]
print(f"\n   >=75% of its cutoffs (strict)      {len(strict):>4d} proteins")
print(f"   >=50% of its cutoffs (permissive)  {len(permissive):>4d} proteins")
print(f"   median protein clears {100 * m.frac_pass.median():.0f}% of its cutoffs")

# Is the strict group just the high-consensus-rank group, or is it also a
# low-disagreement group? Two different claims, so test both.
rest = m[m.frac_pass < 0.75]
u_rank = stats.mannwhitneyu(strict.mean_rank.dropna(), rest.mean_rank.dropna())
u_sd = stats.mannwhitneyu(strict.rank_sd.dropna(), rest.rank_sd.dropna())
print(f"\n   strict group vs rest, consensus rank:  "
      f"{strict.mean_rank.median():.3f} vs {rest.mean_rank.median():.3f}  P = {u_rank.pvalue:.1e}")
print(f"   strict group vs rest, cross-tool SD:   "
      f"{strict.rank_sd.median():.3f} vs {rest.rank_sd.median():.3f}  P = {u_sd.pvalue:.2f}")
print("   -> high-confidence is a RANK phenomenon, not a low-disagreement one:")
print("      the tools argue about the top proteins as much as about the rest.")

# Topology skew in the shortlist. Given the documented hydropathy and length
# confounds, this is better read as a property of the predictors than of biology.
#
# HEADS UP -- this number depends on how you count membrane segments, and the project
# has two counts that disagree:
#   n_transmem (here, and the canonical master's `topology` column) counts TRANSMEM
#     features only.
#   TMD_count (in rank_consensus_table.csv) counts TRANSMEM *and* INTRAMEM.
# They differ for ~8% of the set. The shortlist is where it bites: TGO1 (Q5JRA6) has
# one TRANSMEM plus one INTRAMEM, so it is single-pass by the first count and
# multi-pass by the second. This script therefore reports 15 of 15 single-pass, while
# a version built on TMD_count reports 14 of 15. Both are defensible; they are
# answering slightly different questions ("how many times does it cross the membrane"
# vs "how many membrane-associated segments does it have"). State which you used.
t = m[m.topology != "intramembrane-only"]
ts = t[t.frac_pass >= 0.75]
tab = [[int((ts.topology == "single-pass").sum()), int((ts.topology == "multi-pass").sum())],
       [int((t[t.frac_pass < 0.75].topology == "single-pass").sum()),
        int((t[t.frac_pass < 0.75].topology == "multi-pass").sum())]]
fisher = stats.fisher_exact(tab)
print(f"\n   strict group topology: {tab[0][0]} single-pass / {tab[0][1]} multi-pass"
      f"   vs rest {tab[1][0]}/{tab[1][1]}   Fisher P = {fisher.pvalue:.4f}")

per_tool.to_csv(C.TABLES / "per_tool_pass_rates.csv", index=False)
cols = ["acc", "Entry_name", "n_pass", "n_eval", "frac_pass", "mean_rank", "rank_sd",
        "topology", "n_transmem", "Length"]
(strict.sort_values("frac_pass", ascending=False)[[c for c in cols if c in strict.columns]]
 .to_csv(C.TABLES / "shortlist_75pct.csv", index=False))
(permissive.sort_values("frac_pass", ascending=False)[[c for c in cols if c in permissive.columns]]
 .to_csv(C.TABLES / "shortlist_50pct.csv", index=False))

print(f"\nwrote 3 tables to {C.TABLES}")
