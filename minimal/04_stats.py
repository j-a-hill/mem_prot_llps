"""
STEP 4 -- the statistics.

Three questions, in order:

  A. Do the predictors rank single-pass proteins above multi-pass ones?
  B. Is that a topology effect, or is it hydropathy wearing topology's clothes?
  C. Which tools drive it, one tool at a time?

Every test reports an EFFECT SIZE with a confidence interval, not just a P value,
and the assumption behind the test is checked rather than assumed. Results go to
minimal/tables/*.csv.

Run:  python 04_stats.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C

RNG = np.random.default_rng(0)      # fixed seed: reruns give identical CIs

m = pd.read_csv(C.BUILD / "master.csv")

# Flip the sign of the backwards-running tools (see config.INVERTED) so that for
# every predictor in this script, higher = more likely to phase separate. Without
# this, LLPhyScore's direction reads backwards against all the others.
for col in C.INVERTED:
    if col in m.columns:
        m[col] = -m[col]

# Drop the intramembrane-only proteins: they span no membrane, so 'single vs
# multi-pass' is not defined for them. 12 proteins.
m = m[m.topology != "intramembrane-only"].copy()
single = m[m.topology == "single-pass"]
multi = m[m.topology == "multi-pass"]
print(f"comparing {len(single)} single-pass vs {len(multi)} multi-pass proteins\n")


# ----------------------------------------------------------------- helper functions

def rank_biserial(a, b, n_boot=2000):
    """
    Effect size for a Mann-Whitney test, with a bootstrap 95% CI.

    It is the probability that a random member of `a` outranks a random member of
    `b`, rescaled to -1..+1.  0 means no difference. +1 means every a beats every b.
    This is the number to quote -- a P value only tells you the difference is not
    zero, not how big it is.
    """
    a, b = np.asarray(a), np.asarray(b)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]

    def rb(x, y):
        u = stats.mannwhitneyu(x, y, alternative="two-sided").statistic
        return 2 * u / (len(x) * len(y)) - 1

    point = rb(a, b)
    boot = [rb(RNG.choice(a, len(a), replace=True),
               RNG.choice(b, len(b), replace=True)) for _ in range(n_boot)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
    return point, lo, hi, p


def bh(pvals):
    """Benjamini-Hochberg FDR correction, for when we test 18 tools at once."""
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    n = len(p)
    adj = np.empty(n)
    adj[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


# ============================================================ A. the headline test
print("A. CONSENSUS RANK BY TOPOLOGY")
print("-" * 62)

a, b = single.mean_rank.dropna(), multi.mean_rank.dropna()
rb, lo, hi, p = rank_biserial(a, b)

print(f"  single-pass median   {a.median():.3f}   (n={len(a)})")
print(f"  multi-pass  median   {b.median():.3f}   (n={len(b)})")
print(f"  Mann-Whitney         p = {p:.2e}")
print(f"  rank-biserial        {rb:+.3f}   95% CI [{lo:+.3f}, {hi:+.3f}]")

# A Mann-Whitney test reads as a shift in location only if the two groups have
# similar spread. If the spreads differ, a 'significant' result can be a difference
# in shape instead. Fligner-Killeen tests that, and Brunner-Munzel is the version
# of the test that does not need equal spread.
fk = stats.fligner(a, b).pvalue
bm = stats.brunnermunzel(a, b).pvalue
print(f"\n  assumption checks")
print(f"    equal spread (Fligner-Killeen)  p = {fk:.3f}"
      f"   {'ok' if fk > 0.05 else 'SPREADS DIFFER -- read with care'}")
print(f"    Brunner-Munzel (spread-free)    p = {bm:.2e}"
      f"   {'agrees' if (bm < 0.05) == (p < 0.05) else 'DISAGREES'}")


# ================================================== B. is hydropathy the real cause
# Single-pass and multi-pass proteins differ in hydropathy, and hydropathy is what
# most of these predictors are sensitive to. So the topology difference may be a
# hydropathy difference in disguise.
#
# The test: split the proteins into five hydropathy bins. Inside a bin the two
# groups have similar hydropathy, so if topology matters on its own the gap should
# survive within bins. If it vanishes, hydropathy was doing the work.
print("\n\nB. IS IT HYDROPATHY RATHER THAN TOPOLOGY?")
print("-" * 62)

d = m.dropna(subset=["mean_rank", "whole_hydropathy"]).copy()
d["bin"] = pd.qcut(d.whole_hydropathy, 5, labels=[f"Q{i}" for i in range(1, 6)])

rows = []
for name, g in d.groupby("bin", observed=True):
    s, mu = g[g.topology == "single-pass"].mean_rank, g[g.topology == "multi-pass"].mean_rank
    if len(s) >= 5 and len(mu) >= 5:
        pv = stats.mannwhitneyu(s, mu).pvalue
        rows.append({
            "hydropathy_bin": name, "n_single": len(s), "n_multi": len(mu),
            "median_single": s.median(), "median_multi": mu.median(),
            "difference": s.median() - mu.median(), "p": pv,
        })
strat = pd.DataFrame(rows)
print(strat.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

n_favouring_single = int((strat.difference > 0).sum())
print(f"\n  bins where single-pass still ranks higher: {n_favouring_single} of {len(strat)}")
print("  (all 5 -> topology effect is real; ~half -> the whole-set gap was hydropathy)")

# The same question asked a second way: does the topology-rank correlation survive
# once hydropathy is held constant? Partial Spearman, done by correlating the
# residuals of each variable after regressing hydropathy out.
d["is_single"] = (d.topology == "single-pass").astype(int)


def partial_spearman(x, y, z):
    """Spearman correlation of x and y with the effect of z removed from both."""
    xr, yr, zr = (stats.rankdata(v) for v in (x, y, z))
    rx = xr - np.polyval(np.polyfit(zr, xr, 1), zr)
    ry = yr - np.polyval(np.polyfit(zr, yr, 1), zr)
    return stats.spearmanr(rx, ry)


raw = stats.spearmanr(d.is_single, d.mean_rank)
adj = partial_spearman(d.is_single, d.mean_rank, d.whole_hydropathy)
print(f"\n  topology vs rank, unadjusted        rho = {raw.statistic:+.3f}  p = {raw.pvalue:.2e}")
print(f"  topology vs rank, hydropathy held   rho = {adj.statistic:+.3f}  p = {adj.pvalue:.2e}")

strat.to_csv(C.TABLES / "hydropathy_stratified.csv", index=False)


# ==================================================== C. one tool at a time
# The consensus rank averages 18 tools. If the topology gap comes from a handful of
# them, that is a fact about those tools, not about membrane biology. So test each.
print("\n\nC. TOPOLOGY EFFECT PER PREDICTOR")
print("-" * 62)

rows = []
for col in [c for c in C.PREDICTORS if c in m.columns]:
    s, mu = single[col].dropna(), multi[col].dropna()
    if len(s) < 20 or len(mu) < 20:
        continue
    rb_, lo_, hi_, p_ = rank_biserial(s, mu, n_boot=1000)
    rows.append({
        "predictor": C.NICE[col], "n_single": len(s), "n_multi": len(mu),
        "rank_biserial": rb_, "ci_lo": lo_, "ci_hi": hi_, "p": p_,
    })

per = pd.DataFrame(rows)
per["p_adj_BH"] = bh(per.p)          # 18 tests, so correct for multiplicity
per["direction"] = np.where(per.p_adj_BH > 0.05, "no difference",
                            np.where(per.rank_biserial > 0, "favours single-pass",
                                     "favours multi-pass"))
per = per.sort_values("rank_biserial", ascending=False)

print(per.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
print("\n  " + per.direction.value_counts().to_string().replace("\n", "\n  "))

per.to_csv(C.TABLES / "per_predictor_topology.csv", index=False)


# ------------------------------------------------------------------ headline summary
summary = pd.DataFrame([{
    "n_single": len(a), "n_multi": len(b),
    "median_single": a.median(), "median_multi": b.median(),
    "mannwhitney_p": p, "rank_biserial": rb, "ci_lo": lo, "ci_hi": hi,
    "fligner_p": fk, "brunnermunzel_p": bm,
    "rho_unadjusted": raw.statistic, "rho_hydropathy_held": adj.statistic,
    "bins_favouring_single": n_favouring_single, "n_bins": len(strat),
}])
summary.to_csv(C.TABLES / "headline_stats.csv", index=False)

print(f"\nwrote 3 tables to {C.TABLES}")
