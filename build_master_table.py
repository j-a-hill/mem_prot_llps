"""
build_master_table.py

Phase 0a: Build the provenance master table for all experimental LLPS proteins.

One row per UniProt accession (canonical). Columns:
  UniProt_ID, Protein_Name, Gene_Names, Organism,
  Source_DB (PhasePDB / LLPSDB / both),
  LLPS_mode (autonomous / partner_dependent / both / unclear),
  evidence_in_vitro (bool), evidence_in_vivo (bool),
  role (driver / client / both / unclear),
  MLO_Types,
  is_membrane, TMD_count, TMD_class (none/single/multi),
  Seq_Length, Disordered_fraction (if available),
  LLPSDB_designed (bool — True = synthetic/designed protein, exclude from benchmark)

LLPS mode, evidence type, and role are parsed best-effort from PhasePDB free text.
LLPSDB provides no structured evidence-type column — flagged as unknown.

Output: output/master_table.csv
"""

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parent
OUT  = ROOT / "output"

# ── Load source databases ─────────────────────────────────────────────────────

phasepdb = pd.read_csv(ROOT / "exp_db" / "phasepdb_summary_database_2026-06-15.csv")
phasepdb = phasepdb.rename(columns={"UniProt ID": "UniProt_ID"})

llpsdb = pd.read_excel(ROOT / "exp_db" / "protein_LLPSDB.xls")
llpsdb = llpsdb.rename(columns={"Uniprot ID": "UniProt_ID"})
llpsdb_human = llpsdb[llpsdb["Species"].str.contains("Homo sapiens", na=False)].copy()

fd = pd.read_csv(OUT / "full_dataset.csv")[["Entry", "TMD_count", "Length"]]
fd = fd.rename(columns={"Entry": "UniProt_ID"})

# ── Text-parsing helpers ──────────────────────────────────────────────────────

def _any(text, patterns):
    if not isinstance(text, str):
        return False
    return any(re.search(p, text, re.IGNORECASE) for p in patterns)

AUTONOMOUS_PATTERNS = [
    r"PS-self",
    r"class_:\s*PS-self",
    r"self.driven.*phase sep",
    r"autonomous.*phase sep",
    r"self.assembl.*phase sep",
    r"self.phase.separat",
]
PARTNER_PATTERNS = [
    r"PS-client",
    r"PS-other",
    r"partner.dependent",
    r"partner.induced",
    r"\bclient\b",
    r"scaffold.dependent",
    r"co.phase.separat",   # co-phase separation with a partner
]
INVITRO_PATTERNS  = [r"\bin vitro\b"]
INVIVO_PATTERNS   = [r"\bin vivo\b", r"\bin cells?\b", r"\bin living\b", r"within the (nucleus|cytoplasm|cell)"]
DRIVER_PATTERNS   = [
    r"\bscaffold\b", r"\bdriver\b", r"nucleat.*condensate",
    r"drives.*phase sep", r"PS-self",
    r"self.driven.*phase sep",
]
CLIENT_PATTERNS   = [
    r"\bclient\b", r"\bpassenger\b",
    r"recruited (to|into).*condensate",
    r"partitioned into",
    r"enriched.*within.*condensate",
]

def parse_phasepdb_row(row):
    # Combine the two main free-text fields
    text = " ".join(filter(None, [
        row.get("Overall Phase Separation Overview", ""),
        row.get("Intrinsic Determinants & Regulations", ""),
    ]))

    is_auto   = _any(text, AUTONOMOUS_PATTERNS)
    is_part   = _any(text, PARTNER_PATTERNS)
    is_iv     = _any(text, INVITRO_PATTERNS)
    is_vivo   = _any(text, INVIVO_PATTERNS)
    is_driver = _any(text, DRIVER_PATTERNS)
    is_client = _any(text, CLIENT_PATTERNS)

    if is_auto and is_part:
        mode = "both"
    elif is_auto:
        mode = "autonomous"
    elif is_part:
        mode = "partner_dependent"
    else:
        mode = "unclear"

    if is_driver and is_client:
        role = "both"
    elif is_driver:
        role = "driver"
    elif is_client:
        role = "client"
    else:
        role = "unclear"

    return pd.Series({
        "LLPS_mode":         mode,
        "evidence_in_vitro": is_iv,
        "evidence_in_vivo":  is_vivo,
        "role":              role,
    })

# ── Parse PhasePDB ────────────────────────────────────────────────────────────

parsed = phasepdb.apply(parse_phasepdb_row, axis=1)
phasepdb_parsed = pd.concat([
    phasepdb[["UniProt_ID", "Protein Name", "Gene Names", "Organism", "MLO Types"]].rename(
        columns={"Protein Name": "Protein_Name", "Gene Names": "Gene_Names", "MLO Types": "MLO_Types"}
    ),
    parsed,
], axis=1)
phasepdb_parsed["Source_DB"] = "PhasePDB"
phasepdb_parsed["LLPSDB_designed"] = False

# ── Build LLPSDB records ──────────────────────────────────────────────────────

llpsdb_records = []
for _, r in llpsdb_human.iterrows():
    uid = str(r["UniProt_ID"]).strip() if pd.notna(r["UniProt_ID"]) else None
    if not uid or uid == "nan":
        continue
    llpsdb_records.append({
        "UniProt_ID":        uid,
        "Protein_Name":      r.get("Protein name", ""),
        "Gene_Names":        r.get("Gene name", ""),
        "Organism":          "Homo sapiens",
        "MLO_Types":         r.get("Localization", ""),
        "LLPS_mode":         "unclear",        # no structured field in LLPSDB
        "evidence_in_vitro": False,
        "evidence_in_vivo":  False,
        "role":              "unclear",
        "Source_DB":         "LLPSDB",
        "LLPSDB_designed":   str(r.get("Protein type (N/D)", "N")).strip().upper() == "D",
    })
llpsdb_df = pd.DataFrame(llpsdb_records)

# ── Merge: PhasePDB + LLPSDB → one row per UniProt_ID ────────────────────────

all_ids = set(phasepdb_parsed["UniProt_ID"]) | set(llpsdb_df["UniProt_ID"])
rows = []

for uid in all_ids:
    in_ppdb  = uid in set(phasepdb_parsed["UniProt_ID"])
    in_llps  = uid in set(llpsdb_df["UniProt_ID"])

    if in_ppdb and in_llps:
        base = phasepdb_parsed[phasepdb_parsed["UniProt_ID"] == uid].iloc[0].to_dict()
        base["Source_DB"] = "LLPSDB+PhasePDB"
        # Trust PhasePDB parsed fields; keep LLPSDB_designed=False (natural)
    elif in_ppdb:
        base = phasepdb_parsed[phasepdb_parsed["UniProt_ID"] == uid].iloc[0].to_dict()
    else:
        base = llpsdb_df[llpsdb_df["UniProt_ID"] == uid].iloc[0].to_dict()

    rows.append(base)

master = pd.DataFrame(rows)

# ── Add TMD annotation ────────────────────────────────────────────────────────

master = master.merge(fd, on="UniProt_ID", how="left")
master["TMD_count"]  = pd.to_numeric(master["TMD_count"], errors="coerce").fillna(0).astype(int)
master["is_membrane"] = master["TMD_count"] > 0
master["TMD_class"] = master["TMD_count"].apply(
    lambda n: "none" if n == 0 else ("single" if n == 1 else "multi")
)

# ── Column order ──────────────────────────────────────────────────────────────

cols = [
    "UniProt_ID", "Protein_Name", "Gene_Names", "Organism",
    "Source_DB", "LLPSDB_designed",
    "LLPS_mode", "evidence_in_vitro", "evidence_in_vivo", "role",
    "MLO_Types",
    "is_membrane", "TMD_count", "TMD_class",
    "Length",
]
master = master[[c for c in cols if c in master.columns]]
master = master.sort_values(["is_membrane", "UniProt_ID"], ascending=[False, True]).reset_index(drop=True)

# ── Save ──────────────────────────────────────────────────────────────────────

master.to_csv(OUT / "master_table.csv", index=False)
print(f"Saved output/master_table.csv  ({len(master)} rows × {len(master.columns)} cols)")

# ── Summary ───────────────────────────────────────────────────────────────────

print("\n── Source breakdown ──────────────────────────────────────────────────")
print(master["Source_DB"].value_counts().to_string())

print("\n── LLPS mode (PhasePDB proteins only) ───────────────────────────────")
ppdb_only = master[master["Source_DB"].isin(["PhasePDB", "LLPSDB+PhasePDB"])]
print(ppdb_only["LLPS_mode"].value_counts().to_string())

print("\n── Evidence (PhasePDB proteins only) ────────────────────────────────")
print(f"  in vitro only:  {(ppdb_only['evidence_in_vitro'] & ~ppdb_only['evidence_in_vivo']).sum()}")
print(f"  in vivo only:   {(~ppdb_only['evidence_in_vitro'] & ppdb_only['evidence_in_vivo']).sum()}")
print(f"  both:           {(ppdb_only['evidence_in_vitro'] & ppdb_only['evidence_in_vivo']).sum()}")
print(f"  neither parsed: {(~ppdb_only['evidence_in_vitro'] & ~ppdb_only['evidence_in_vivo']).sum()}")

print("\n── Role (PhasePDB proteins only) ────────────────────────────────────")
print(ppdb_only["role"].value_counts().to_string())

print("\n── Membrane subset ──────────────────────────────────────────────────")
mem = master[master["is_membrane"]]
print(f"  Total membrane: {len(mem)}")
print(f"  TMD class:"); print(f"  {mem['TMD_class'].value_counts().to_dict()}")
print(f"  LLPS mode:"); print(f"  {mem['LLPS_mode'].value_counts().to_dict()}")
print(f"  LLPSDB designed (exclude): {mem['LLPSDB_designed'].sum()}")
print(f"  Evidence in vitro: {mem['evidence_in_vitro'].sum()}")
print(f"  Evidence in vivo:  {mem['evidence_in_vivo'].sum()}")

print("\n── LLPSDB designed proteins (should exclude from benchmark) ─────────")
designed = master[master["LLPSDB_designed"]]
print(designed[["UniProt_ID", "Protein_Name"]].to_string(index=False))
