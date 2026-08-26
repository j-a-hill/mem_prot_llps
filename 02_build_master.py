"""
STEP 2 -- wrangle the five databases into one master table.

The logic, in order:

  1. Read each database. Keep human entries only. Pull out the accession and that
     database's own role/type word.
  2. Union the five accession lists  ->  every human protein any database calls
     phase-separating.
  3. Keep the ones UniProt annotates as membrane proteins (>=1 TRANSMEM or
     INTRAMEM feature). That intersection IS the study set.
  4. Join on the precomputed predictor scores, sequence features and consensus rank.
  5. Write minimal/build/master.csv, one row per protein.

Every filtering step prints how many proteins went in and how many came out, so the
set size is derived in front of you rather than asserted.

Run:  python 02_build_master.py
"""

import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))   # find config.py
import config as C

R = C.RAW


def show(label, n):
    print(f"  {label:<52s} {n:>6d}")


def read_table(path):
    """
    Read a delimited text file without caring whether it is comma- or tab-separated.

    Needed because the same logical file arrives in different formats depending on
    where it came from: UniProt's API serves TSV while the stored snapshot is CSV, and
    DrLLPS ships a tab-separated file named .csv (with one malformed line).
    """
    return pd.read_csv(path, sep=None, engine="python", on_bad_lines="skip")


def read_json_records(path):
    """Read a JSON file that holds a list of records, wrapped or bare."""
    obj = json.load(open(path))
    if isinstance(obj, dict) and "data" in obj:
        return pd.DataFrame(obj["data"])
    return pd.DataFrame(obj)


print("READING DATABASES (human entries only)")
print("-" * 62)

# ------------------------------------------------------------------ 1. CD-CODE
# A JSON list of records, one per protein-in-a-condensate.
cdcode = read_json_records(R / "cdcode_proteins.json")
cdcode = cdcode[cdcode.species_name == "Homo sapiens"]
cdcode = cdcode.rename(columns={"uniprot_id": "acc"})
# How many distinct condensates each protein was seen in.
cdcode_n = cdcode.groupby("acc").biomolecular_condensate_count.max()
show("CD-CODE human proteins", cdcode.acc.nunique())

# CD-CODE role words (driver / member / unknown) live in separate per-protein crawls,
# because the bulk protein list does not carry them. Three files were crawled at
# different times and they overlap, so we merge them; later files win.
# This is an enrichment of the table, not an inclusion criterion.


def cdcode_role_word(entry):
    """Pull one role word out of a CD-CODE role record. 'driver' outranks the rest."""
    if not isinstance(entry, dict):
        return None
    r = entry.get("roles") or entry.get("functional_types") or []
    r = list(r.keys()) if isinstance(r, dict) else list(r)
    if not r:
        return None
    return "driver" if "driver" in r else sorted(r)[0]


role_files = [
    C.SOURCES_DIR / "cdcode_protein_roles.json",
    C.SOURCES_DIR / "cdcode_addition_driver_annotations.json",
    C.SOURCES_DIR / "cdcode_roles_detail.json",
]
role_map = {}
for f in role_files:
    for a, v in json.load(open(f)).items():
        w = cdcode_role_word(v)
        if w:
            role_map[a] = w
cdcode_role = pd.Series(role_map, name="CDCODE_role")
show("  ... with a CD-CODE role word on record", len(role_map))

# ----------------------------------------------------------------- 2. PhaSepDB
# Two files, because they carry different things. The summary export gives membership
# and MLO types; the API gives the PS-self / PS-other class.
phasepdb = read_table(R / "phasepdb.csv").rename(columns={"UniProt ID": "acc"})
show("PhaSepDB human proteins", phasepdb.acc.nunique())

# The class_ field, read from the API export -- one row per curated PMID entry, so a
# protein studied in several papers appears several times. See config.PHASEPDB_CLASS
# for why this must not be scraped out of the prose summary column instead.
cls = read_table(R / "phasepdb_class.csv")
per_entry = cls.groupby("uniprot_id")["class_"].agg(
    PhaSepDB_n_PS_self_entries=lambda s: int((s == "PS-self").sum()),
    PhaSepDB_n_PS_other_entries=lambda s: int((s == "PS-other").sum()),
)

# Report the literal tags a protein carries. A protein with both is written
# "PS-self;PS-other" rather than being collapsed into an invented third category.
per_entry["PhaSepDB_class"] = [
    "PS-self;PS-other" if a and b else "PS-self" if a else "PS-other"
    for a, b in zip(per_entry.PhaSepDB_n_PS_self_entries,
                    per_entry.PhaSepDB_n_PS_other_entries)
]
# A PhaSepDB driver call = at least one entry curated as PS-self.
per_entry["PhaSepDB_is_driver"] = per_entry.PhaSepDB_n_PS_self_entries >= 1

show("  ... with a class_ tag on the API", int(per_entry.index.isin(phasepdb.acc).sum()))

# ------------------------------------------------------------------- 3. DrLLPS
# Tab-separated despite the .csv name. One row per protein-condensate pair, so a
# protein can appear several times; we keep the strongest role it is ever given.
drllps = read_table(R / "drllps.tsv").rename(columns={"UniProt ID": "acc"})
drllps = drllps[drllps.Species == "Homo sapiens"]
rank = {"Scaffold": 0, "Regulator": 1, "Client": 2}   # Scaffold is the strongest
drllps_type = (drllps.assign(r=drllps["LLPS Type"].map(rank))
                     .sort_values("r")
                     .drop_duplicates("acc")
                     .set_index("acc")["LLPS Type"]
                     .rename("DrLLPS_type"))
show("DrLLPS human proteins", drllps.acc.nunique())

# ----------------------------------------------------------------- 4. PhaSePro
# JSON keyed by accession. partner_dep tells you whether it needs a partner.
phasepro = pd.DataFrame(json.load(open(R / "phasepro.json")).values())
phasepro = phasepro[phasepro.organism.str.contains("Homo sapiens", na=False)]
phasepro = phasepro.rename(columns={"accession": "acc"})
phasepro_dep = phasepro.set_index("acc")["partner_dep"].rename("PhaSePro_partner_dep")
show("PhaSePro human proteins", phasepro.acc.nunique())

# ------------------------------------------------------------------- 5. LLPSDB
# An .xls of in-vitro experiments. No role field, so it contributes membership only.
llpsdb = pd.read_excel(R / "llpsdb.xls").rename(columns={"Uniprot ID": "acc"})
llpsdb = llpsdb[llpsdb.Species == "Homo sapiens"]
show("LLPSDB human proteins", llpsdb.acc.nunique())


# --------------------------------------------------------- 6. union the five lists
print("\nBUILDING THE SET")
print("-" * 62)

members = {
    "CD-CODE": set(cdcode.acc),
    "PhaSepDB": set(phasepdb.acc),
    "DrLLPS": set(drllps.acc),
    "PhaSePro": set(phasepro.acc),
    "LLPSDB": set(llpsdb.acc),
}

union = set().union(*members.values())
# Drop anything that is not a well-formed UniProt accession (blanks, '-', notes).
union = {a for a in union if isinstance(a, str) and re.fullmatch(r"[A-Z0-9]{6,10}", a)}
show("union of all five databases", len(union))


# ------------------------------------------------- 7. keep the membrane proteins
# UniProt writes membrane annotation as feature strings. Counting the number of
# 'TRANSMEM' occurrences gives the number of membrane-spanning segments.
tm = read_table(R / "uniprot_tm.csv").rename(columns={"Entry": "acc"})
tm["n_transmem"] = tm.Transmembrane.fillna("").str.count("TRANSMEM")
tm["n_intramem"] = tm.Intramembrane.fillna("").str.count("INTRAMEM")

show("proteins with UniProt membrane annotation", len(tm))
show("  ... of which have >=1 TRANSMEM", int((tm.n_transmem > 0).sum()))
show("  ... of which have >=1 INTRAMEM", int((tm.n_intramem > 0).sum()))

missing = len(union - set(tm.acc))
show("union members absent from the TM file", missing)

membrane = tm[(tm.n_transmem > 0) | (tm.n_intramem > 0)]
keep = union & set(membrane.acc)
show("MEMBRANE LLPS PROTEINS (the study set)", len(keep))


# ------------------------------------------------------------ 8. assemble the table
master = membrane[membrane.acc.isin(keep)][["acc", "n_transmem", "n_intramem"]].copy()

# Topology class: one spanning segment vs several. Proteins with only an
# intramembrane feature span no membrane at all, so they get their own label.
master["topology"] = [
    "intramembrane-only" if t == 0 else "single-pass" if t == 1 else "multi-pass"
    for t in master.n_transmem
]

# Which databases each protein came from, and how many.
for db, accs in members.items():
    master[f"in_{db}"] = master.acc.isin(accs)
master["n_databases"] = master[[f"in_{d}" for d in members]].sum(axis=1)

# Each database's own role word, left verbatim -- no harmonising into shared
# categories here, because the words do not mean the same thing across databases.
master = (master
          .join(cdcode_role, on="acc")
          .join(cdcode_n.rename("CDCODE_n_condensates"), on="acc")
          .merge(phasepdb[["acc", "MLO Types", "Gene Names"]].drop_duplicates("acc"),
                 on="acc", how="left")
          .join(per_entry[["PhaSepDB_class", "PhaSepDB_is_driver",
                           "PhaSepDB_n_PS_self_entries",
                           "PhaSepDB_n_PS_other_entries"]], on="acc")
          .join(drllps_type, on="acc")
          .join(phasepro_dep, on="acc"))


# ------------------------------- 9. join the predictor scores and sequence features
print("\nJOINING PRECOMPUTED PREDICTOR OUTPUT")
print("-" * 62)

scores = pd.read_csv(C.SCORES).rename(columns={"Entry": "acc"})
score_cols = ["acc"] + [c for c in C.PREDICTORS if c in scores.columns]
show("predictor score rows available", len(scores))

feats = pd.read_csv(C.FEATURES).rename(columns={"UniProt_ID": "acc"})
cons = pd.read_csv(C.CONSENSUS).rename(columns={"UniProt_ID": "acc"})

master = (master
          .merge(scores[score_cols], on="acc", how="left")
          .merge(feats, on="acc", how="left")
          .merge(cons[["acc", "mean_rank", "rank_sd", "n_predictors",
                       "Entry_name", "Length", "leakage_any"]],
                 on="acc", how="left"))

show("proteins with a consensus rank", int(master.mean_rank.notna().sum()))
show("proteins with no consensus rank", int(master.mean_rank.isna().sum()))


# ---------------------------------------------- 10. compare against the paper's set
canon = set(pd.read_csv(C.CANONICAL).UniProt_ID)
print("\nCHECK AGAINST THE PUBLISHED 475-PROTEIN SET")
print("-" * 62)
show("this rebuild", len(keep))
show("published set", len(canon))
show("in both", len(keep & canon))
show("published only (not recovered here)", len(canon - keep))
show("rebuild only (new)", len(keep - canon))

# Name the misses and say why, so the gap is explained rather than left hanging.
# Two causes, both to do with the snapshots being dated:
#   - absent from the LLPS databases:  the published set used an earlier CD-CODE
#     crawl that still listed them; the current snapshot does not.
#   - absent from the UniProt TM file: that file covers REVIEWED entries only, so
#     an unreviewed accession has no membrane annotation to match on.
for acc in sorted(canon - keep):
    if acc not in union:
        why = "not in any current database snapshot"
    elif acc not in set(tm.acc):
        why = "no UniProt membrane annotation (unreviewed entry)"
    else:
        why = "in the databases and TM file but failed the membrane filter"
    print(f"    {acc}: {why}")

master = master.sort_values("acc").reset_index(drop=True)
out = C.BUILD / "master.csv"
master.to_csv(out, index=False)
print(f"\nwrote {out}   {master.shape[0]} rows x {master.shape[1]} columns")
