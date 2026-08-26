"""
STEP 9 -- derive two of the shipped inputs/ tables from the files upstream of them.

Not part of the everyday run; run_all.py skips it. Its job is provenance: it shows how
two files in inputs/ are computed, so they are not black boxes, and it CHECKS its own
output against the shipped copies rather than asking you to trust it.

  rank_consensus_table.csv   <- predictor_comparison.csv
  clean_masks.csv            <- leakage_map.csv

The chain has one more step above these, and it cannot run here:

  Predictors_whole_genome_sets/   ~780 MB of raw per-tool proteome scores, produced by
                                  running each of the 18 predictors locally -- some
                                  need a GPU, some are not redistributable
    -> background_scored.csv      (wrangle_background.py, in the full project)
    -> predictor_comparison.csv   (the per-set slice of that)

So this starts one level down, from tables that ARE shipped. That is the honest
boundary: below it is arithmetic you can read, above it is a re-run of the predictors.

Run:  python 09_derive_inputs.py
"""

import sys
from pathlib import Path

import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C


def normalised_ranks(scores):
    """
    Turn each predictor's raw scores into a 0..1 rank, highest score = 1.0.

    Ranking is what makes 18 tools comparable at all -- their scales are wildly
    different (probabilities, log-likelihood ratios, unbounded physics scores).
    Ranked within each tool, they share one axis.

    Each column is ranked over the proteins that tool actually scored, so a tool that
    failed on 30 proteins is not penalised for the gap; the per-protein average then
    uses however many tools scored that protein.
    """
    out = pd.DataFrame(index=scores.index, dtype=float)
    for col in [c for c in C.PREDICTORS if c in scores.columns]:
        v = scores[col]
        ok = v.notna()
        vals = v[ok].to_numpy()
        if len(vals) == 0:
            continue
        # rankdata gives rank 1 to the SMALLEST value, so negate to put the highest
        # score at rank 1, then rescale so rank 1 -> 1.0 and the last -> 0.0.
        raw = stats.rankdata(-vals, method="average")
        out.loc[ok, C.NICE[col]] = (1.0 - (raw - 1) / (len(vals) - 1)
                                    if len(vals) > 1 else 0.5)
    return out


def three_state(series):
    """
    Normalise a leakage column to True / False / NA, mapping 'unknown' to NA.

    THIS FUNCTION IS THE POINT of the clean-mask code. The leakage columns are
    three-state -- a protein is in a tool's training set, is not, or nobody could
    determine it -- and the states arrive in MIXED TYPES: most columns hold real
    booleans, but LLPhyScore's hold the STRINGS 'True' and 'unknown', because its
    training set is keyed by gene/construct and mostly unmappable to UniProt.

    Two ways to get this wrong, both of which have bitten this project:
      * `col == False` to find clean proteins. In pandas `'unknown' == False` is
        False, so unknown rows fall out of the clean set too -- which collapsed
        LLPhyScore's clean positives to zero.
      * `col != True` on the raw column. The string `'True' != True` is True, so two
        genuinely-leaked proteins (P08908, O60500) get called clean.
    Normalising the types FIRST, then asking "confirmed leaked?", avoids both.
    """
    def one(v):
        s = str(v).strip()
        if v is True or s == "True":
            return True
        if v is False or s == "False":
            return False
        return pd.NA
    return series.map(one)


def clean_mask(leak, tool):
    """
    True where a protein is usable as a test case for this tool.

    Clean unless CONFIRMED leaked. 'unknown' (now NA) means nobody could tell, which
    is not evidence of leakage, so it counts as clean -- otherwise a single
    unverifiable tool disqualifies the whole set.
    """
    pos = leak[f"pos_{tool}"].fillna(False).astype(bool)
    neg = leak[f"neg_{tool}"].fillna(False).astype(bool)
    return ~(pos | neg)


def check(name, mine, shipped, tol=1e-9):
    """Compare a derived table against the shipped one; report the worst cell."""
    common = [c for c in mine.columns if c in shipped.columns]
    worst, bad = 0.0, []
    for col in common:
        a, b = mine[col].reindex(shipped.index), shipped[col]
        if pd.api.types.is_bool_dtype(b) or b.dtype == object:
            n = int((a.astype(bool) != b.astype(bool)).sum())
            if n:
                bad.append(f"{col} ({n} rows)")
        else:
            d = (pd.to_numeric(a, errors="coerce")
                 - pd.to_numeric(b, errors="coerce")).abs().max()
            if not pd.isna(d):
                worst = max(worst, d)
                if d > tol:
                    bad.append(f"{col} (max delta {d:.2e})")
    verdict = "matches the shipped copy" if not bad else "DIFFERS: " + "; ".join(bad)
    print(f"  {name:22s} {len(common):2d} cols, {len(mine):3d} rows -> {verdict}")
    if not bad and worst:
        print(f"  {'':22s} worst numeric difference {worst:.1e} (floating-point noise)")
    return not bad


print("DERIVING inputs/ TABLES FROM THEIR UPSTREAM FILES")
print("-" * 74)

# ------------------------------------------------- 1. the consensus rank table
scores = pd.read_csv(C.SCORES).set_index("Entry")

# The shipped rank matrix has the LLPhyScore sign flip ALREADY APPLIED -- its ranks are
# on the "higher = more LLPS-like" convention, like every other tool. Flip here too, or
# LLPhyScore's ranks come out exactly inverted (a 1.00 difference on every protein).
for col in C.INVERTED:
    if col in scores.columns:
        scores[col] = -scores[col]

ranks = normalised_ranks(scores)
consensus = pd.DataFrame({
    "mean_rank": ranks.mean(axis=1, skipna=True),
    "rank_sd": ranks.std(axis=1, ddof=1, skipna=True),
    "n_predictors": ranks.notna().sum(axis=1),
})
consensus.index.name = "UniProt_ID"
ok_rank = check("rank_consensus_table", consensus,
                pd.read_csv(C.CONSENSUS).set_index("UniProt_ID"))

# ------------------------------------------------------------ 2. the clean masks
leak_path = C.INPUTS / "leakage_map.csv"
if not leak_path.exists():
    print(f"  clean_masks            skipped, {leak_path.name} not shipped in inputs/")
    ok_masks = True
else:
    leak = pd.read_csv(leak_path)
    for col in leak.columns:
        if col.startswith(("pos_", "neg_")):
            leak[col] = three_state(leak[col])

    tools = sorted({c[4:] for c in leak.columns if c.startswith("pos_")})

    masks = pd.DataFrame({"UniProt_ID": leak.UniProt_ID})
    # Globally clean = not confirmed leaked for ANY tool, including tools none of our
    # 18 score columns use.
    glob = pd.Series(True, index=leak.index)
    for t in tools:
        glob &= clean_mask(leak, t)
    masks["clean_global"] = glob.values

    for score_col, tool in C.SCORE_TO_TOOL.items():
        masks[f"clean_{score_col}"] = (True if tool is None
                                       else clean_mask(leak, tool).values)

    ok_masks = check("clean_masks", masks.set_index("UniProt_ID"),
                     pd.read_csv(C.CLEAN_MASKS).set_index("UniProt_ID"))

print()
if ok_rank and ok_masks:
    print("Both derivations reproduce the shipped tables, so the arithmetic behind them")
    print("is the arithmetic in this file -- none of it is taken on trust.")
else:
    print("A derivation DIFFERS from its shipped table. Do not treat this script as an")
    print("explanation of the shipped copy until the difference is understood.")
    sys.exit(1)
