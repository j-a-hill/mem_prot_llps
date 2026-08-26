"""
STEP 6 -- the background comparison. This is the paper's central result.

Everything up to here compared membrane LLPS proteins against each other. This script
compares them against the rest of the proteome, which is the only way to ask the
question the paper is actually about.

Two things get computed:

  A. THE OFFSET. Take the LLPS proteins out, leaving plain background proteins, and
     ask: can a tool's score alone tell a membrane protein from a soluble one? Score
     as AUROC. 0.5 means it cannot. Below 0.5 means it scores membrane proteins
     systematically LOWER -- a standing penalty that has nothing to do with phase
     separation, because none of these proteins are annotated as phase-separating.

  B. WHAT THE OFFSET DOES. Score the 473 membrane LLPS positives against a soluble
     background, then against a membrane background, and take the difference. A tool
     with a big offset looks far better on the membrane background -- not because it
     got better, but because the offset it applies to every membrane protein is now
     applied to the negatives too, and cancels.

The two numbers turn out to be almost perfectly anticorrelated. That is the finding:
a tool's apparent skill on membrane proteins is largely predicted by its bias against
them, which is a property of the comparison, not of the biology.

Run:  python 06_background.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C

bg = pd.read_csv(C.BACKGROUND)

# Guard against the legacy 427-protein files (see the note in config.py).
clean = pd.read_csv(C.CLEAN_MASKS)
assert len(clean) == 475, f"clean_masks has {len(clean)} rows, expected the 475-row canonical file"

for col in C.INVERTED:                      # higher = more LLPS-like, for every tool
    if col in bg.columns:
        bg[col] = -bg[col]

is_mem = bg.is_membrane_protein == 1
is_pos = bg.label_mem_bg == 1               # the 473 membrane LLPS positives

print(f"background table: {len(bg)} scored human proteins")
print(f"  membrane LLPS positives          {int(is_pos.sum()):>7d}")
print(f"  membrane background (no LLPS)    {int((is_mem & ~is_pos).sum()):>7d}")
print(f"  soluble background  (no LLPS)    {int((~is_mem & ~is_pos).sum()):>7d}")


def auroc(pos, neg):
    """
    Probability that a random positive scores above a random negative.

    Computed from ranks, which handles ties correctly -- a naive threshold sweep does
    not. 0.5 is chance. Returns NaN if either group is empty.
    """
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    r = stats.rankdata(np.concatenate([pos, neg]))
    return (r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


rows = []
for col in [c for c in C.PREDICTORS if c in bg.columns]:
    s = bg[col]
    mem_bg = s[is_mem & ~is_pos]            # membrane, not phase-separating
    sol_bg = s[~is_mem & ~is_pos]           # soluble, not phase-separating
    positives = s[is_pos]

    # A. the offset: membrane vs soluble among NON-LLPS proteins only
    offset = auroc(mem_bg, sol_bg)

    # B. the same positives, scored against each background
    a_soluble = auroc(positives, sol_bg)
    a_membrane = auroc(positives, mem_bg)

    rows.append({
        "tool": C.NICE[col],
        "offset_auroc": offset,
        "auroc_vs_soluble": a_soluble,
        "auroc_vs_membrane": a_membrane,
        "auroc_shift": a_membrane - a_soluble,
    })

res = pd.DataFrame(rows).sort_values("offset_auroc").reset_index(drop=True)

print("\nA. DOES THE SCORE ALONE SEPARATE MEMBRANE FROM SOLUBLE PROTEINS?")
print("   (LLPS proteins excluded from both groups, so 0.5 = no membrane bias)")
print("-" * 74)
print(res[["tool", "offset_auroc"]].to_string(index=False, float_format=lambda x: f"{x:.3f}"))

n_below = int((res.offset_auroc < 0.5).sum())
print(f"\n   {n_below} of {len(res)} tools score membrane proteins BELOW soluble ones")
print(f"   most biased: {res.iloc[0].tool} at {res.iloc[0].offset_auroc:.3f}")

print("\n\nB. THE SAME POSITIVES, SCORED AGAINST TWO DIFFERENT BACKGROUNDS")
print("-" * 74)
print(res[["tool", "auroc_vs_soluble", "auroc_vs_membrane", "auroc_shift"]]
      .to_string(index=False, float_format=lambda x: f"{x:+.3f}"))

# How many tools cross the 0.5 line -- look useless on one background and skilful on
# the other. This is the flip the paper leads on.
flip = res[(res.auroc_vs_soluble < 0.5) & (res.auroc_vs_membrane > 0.5)]
print(f"\n   {len(flip)} tools sit BELOW chance on the soluble background and ABOVE it")
print(f"   on the membrane background: {', '.join(flip.tool)}")

rho = stats.spearmanr(res.offset_auroc, res.auroc_shift)
print(f"\n   offset vs shift:  Spearman rho = {rho.statistic:+.3f}  (P = {rho.pvalue:.1e})")
print("   i.e. how much a tool gains from the membrane background is almost entirely")
print("   predicted by how much it penalises membrane proteins in the first place.")

res.to_csv(C.TABLES / "background_offset.csv", index=False)
print(f"\nwrote {C.TABLES / 'background_offset.csv'}")
