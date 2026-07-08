"""
wrangle_exp_db.py

Combines two experimental LLPS databases (PhasePDB + LLPSDB) into a single
unified DataFrame of human proteins with experimental LLPS evidence, then
joins all predictor scores from Predictors_whole_genome_sets/.

Outputs:
  output/exp_db_combined.csv          — 882 unique human proteins, unified metadata
  output/predictor_comparison_all.csv — all 882 proteins × predictor scores
  output/predictor_comparison_mem.csv — 60 membrane proteins × predictor scores (subset)
"""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parent
EXP  = ROOT / "exp_db"
PRED = ROOT / "Predictors_whole_genome_sets"
OUT  = ROOT / "output"
OUT.mkdir(exist_ok=True)

# ── 1. Load and harmonise each source ────────────────────────────────────────

pdb = pd.read_csv(EXP / "phasepdb_summary_database_2026-06-15.csv")
pdb = pdb.rename(columns={
    "UniProt ID":   "UniProt_ID",
    "Protein Name": "Protein_Name",
    "Gene Names":   "Gene_Names",
    "MLO Types":    "MLO_Types",
})
pdb = pdb[["UniProt_ID", "Protein_Name", "Gene_Names", "MLO_Types"]].copy()
pdb["Source"] = "PhasePDB"

llps = pd.read_excel(EXP / "protein_LLPSDB.xls")
llps = llps[llps["Species"] == "Homo sapiens"].copy()
llps = llps.rename(columns={
    "Uniprot ID":   "UniProt_ID",
    "Protein name": "Protein_Name",
    "Gene name":    "Gene_Names",
    "Localization": "MLO_Types",   # closest equivalent column
})
llps = llps[["UniProt_ID", "Protein_Name", "Gene_Names", "MLO_Types"]].copy()
llps["Source"] = "LLPSDB"

# ── 2. Union with deduplication ───────────────────────────────────────────────
# Proteins in both DBs: keep PhasePDB row (richer MLO annotation), note both sources.

combined = pd.concat([pdb, llps], ignore_index=True)

# For proteins present in both, collapse rows: merge Source strings, keep first non-null values
def _agg(grp):
    sources = "+".join(sorted(grp["Source"].unique()))
    return pd.Series({
        "Protein_Name": grp["Protein_Name"].dropna().iloc[0] if grp["Protein_Name"].notna().any() else None,
        "Gene_Names":   grp["Gene_Names"].dropna().iloc[0]   if grp["Gene_Names"].notna().any()   else None,
        "MLO_Types":    grp["MLO_Types"].dropna().iloc[0]    if grp["MLO_Types"].notna().any()     else None,
        "Source":       sources,
    })

combined = combined.groupby("UniProt_ID", sort=False).apply(_agg).reset_index()

print(f"Combined experimental DB: {len(combined)} unique human proteins")
print(f"Source breakdown:\n{combined['Source'].value_counts().to_string()}")

# ── 3. Flag membrane proteins ─────────────────────────────────────────────────

mem = pd.read_csv(EXP / "membrane_exp_db_matches.csv")
mem_ids = set(mem["Entry"])
combined["is_membrane"] = combined["UniProt_ID"].isin(mem_ids)
print(f"\nMembrane proteins in combined DB: {combined['is_membrane'].sum()}")

combined.to_csv(OUT / "exp_db_combined.csv", index=False)
print(f"Saved output/exp_db_combined.csv")

# ── 4. Load all predictors ────────────────────────────────────────────────────

print("\nLoading predictors...")

picnic = pd.read_csv(PRED / "PICNIC-9606-data.csv")[["Uniprot ID", "PICNIC score", "PICNIC GO score"]]
picnic = picnic.rename(columns={"Uniprot ID": "UniProt_ID", "PICNIC score": "PICNIC_score", "PICNIC GO score": "PICNIC_GO_score"})
picnic = picnic.groupby("UniProt_ID", sort=False)[["PICNIC_score", "PICNIC_GO_score"]].max().reset_index()

psap = pd.read_excel(PRED / "PSAP_full_proteome_excluding_training_set.xlsx")[["Uniprot.ID", "PSAP_score"]]
psap = psap.rename(columns={"Uniprot.ID": "UniProt_ID"})
psap = psap.drop_duplicates(subset="UniProt_ID")

psphunter = pd.read_excel(PRED / "PSPHunter_full_proteome.xlsx", skiprows=1)[["Uniprot", "Probability"]]
psphunter = psphunter.rename(columns={"Uniprot": "UniProt_ID", "Probability": "PSPHunter_prob"})

pspire = pd.read_csv(PRED / "PSPire_Homo_sapiens_phos_scores.csv")[["Uniprot_ID", "Score", "Include_IDRs"]]
pspire = pspire.rename(columns={"Uniprot_ID": "UniProt_ID", "Score": "PSPire_score", "Include_IDRs": "PSPire_includes_IDRs"})

pspspredict_cols = ["Uniprot", "PDL", "catGRANULE", "PLAAC", "PScore",
                    "ESpritz-DisProt", "SEG", "DeepPhase", "SaPS-10fea",
                    "PdPS-10fea", "Psphunter", "PSPire"]
pspspredict = pd.read_excel(PRED / "PSPspredict_full_proteome.xlsx", sheet_name="PredictAll")[pspspredict_cols]
pspspredict = pspspredict.rename(columns={
    "Uniprot": "UniProt_ID",
    **{c: f"PSPspredict_{c.replace('-', '_').replace(' ', '_')}" for c in pspspredict_cols[1:]}
})

# ── 5. Join predictors onto combined DB ──────────────────────────────────────

df = combined.copy()
for pred_df in [picnic, psap, psphunter, pspire, pspspredict]:
    df = df.merge(pred_df, on="UniProt_ID", how="left")

score_cols = [c for c in df.columns if c not in
              ["UniProt_ID", "Protein_Name", "Gene_Names", "MLO_Types", "Source", "is_membrane"]]

print("\nCoverage (non-null) across all 882 proteins:")
print(df[score_cols].notna().sum().to_string())

df.to_csv(OUT / "predictor_comparison_all.csv", index=False)
print(f"\nSaved output/predictor_comparison_all.csv  ({len(df)} rows × {len(df.columns)} cols)")

# ── 6. Membrane-only subset ───────────────────────────────────────────────────

# Merge in membrane metadata (FuzDrop p(LLPS), pLLPS_Class, TMD_count, etc.)
mem_meta = mem.rename(columns={"Entry": "UniProt_ID"})[
    ["UniProt_ID", "Entry name", "p(LLPS)", "pLLPS_Class", "TMD_count",
     "Length", "Functional_Categories", "Compartment"]
]
df_mem = df[df["is_membrane"]].merge(mem_meta, on="UniProt_ID", how="left")
df_mem.to_csv(OUT / "predictor_comparison_mem.csv", index=False)
print(f"Saved output/predictor_comparison_mem.csv  ({len(df_mem)} rows × {len(df_mem.columns)} cols)")
