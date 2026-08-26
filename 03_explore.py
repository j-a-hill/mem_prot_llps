"""
STEP 3 -- look at the master table before doing any statistics.

Prints what the table contains, how complete each column is, and how the databases
overlap. Nothing here is a result; it is the pass you make to catch a bad join or an
empty column before it turns into a figure.

Run:  python 03_explore.py
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C

pd.set_option("display.width", 130)

m = pd.read_csv(C.BUILD / "master.csv")

print(f"master.csv: {m.shape[0]} proteins x {m.shape[1]} columns\n")

# ------------------------------------------------------------------- composition
print("TOPOLOGY")
print(m.topology.value_counts().to_string(), "\n")

print("MEMBRANE-SPANNING SEGMENTS (n_transmem)")
print(m.n_transmem.describe()[["min", "50%", "max"]].to_string(), "\n")

print("HOW MANY DATABASES EACH PROTEIN IS IN")
print(m.n_databases.value_counts().sort_index().to_string(), "\n")

db_cols = [c for c in m.columns if c.startswith("in_")]
print("PER-DATABASE COVERAGE")
for c in db_cols:
    n = int(m[c].sum())
    print(f"  {c[3:]:<10s} {n:>4d}  ({100 * n / len(m):.1f}%)")

# Which databases contribute proteins nobody else has. A database that contributes
# none is not adding independent evidence to this set.
print("\nPROTEINS UNIQUE TO ONE DATABASE")
for c in db_cols:
    n = int((m[c] & (m.n_databases == 1)).sum())
    print(f"  {c[3:]:<10s} {n:>4d}")

# ------------------------------------------------------------- the role/type words
# Each database uses its own vocabulary. They are printed separately on purpose --
# 'member', 'Client' and 'PS-other' are not the same claim.
print("\nROLE WORDS, VERBATIM PER DATABASE")
for col in ["CDCODE_role", "PhaSepDB_class", "DrLLPS_type", "PhaSePro_partner_dep"]:
    if col in m:
        vc = m[col].value_counts(dropna=False)
        print(f"\n  {col}")
        for k, v in vc.items():
            print(f"    {str(k):<24s} {v:>4d}")

# ---------------------------------------------------------------- completeness
print("\n\nCOLUMN COMPLETENESS (columns with missing values only)")
miss = m.isna().sum()
miss = miss[miss > 0].sort_values(ascending=False)
if len(miss) == 0:
    print("  no missing values")
else:
    for k, v in miss.items():
        print(f"  {k:<26s} {v:>4d} missing  ({100 * v / len(m):.1f}%)")

# ------------------------------------------------------------- predictor coverage
present = [c for c in C.PREDICTORS if c in m.columns]
print(f"\nPREDICTOR SCORE COVERAGE ({len(present)} tools)")
for c in present:
    n = int(m[c].notna().sum())
    flag = "  <-- incomplete" if n < len(m) else ""
    print(f"  {C.NICE[c]:<14s} {n:>4d} of {len(m)}{flag}")

print("\nCONSENSUS RANK")
print(m.mean_rank.describe().to_string())
