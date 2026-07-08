"""
wrangle_background.py

Builds a proteome-wide scored dataset for ROC analysis by joining all
predictor whole-genome files, then assigning label columns.

Label columns
-------------
  label_all       — 1 if in 882 experimental LLPS proteins, else 0
  label_mem       — 1 if in 60 membrane experimental proteins, else 0
                    (negatives = everything NOT in 882 experimental set)
  label_mem_fixed — same as label_mem but excludes all non-membrane
                    non-experimental proteins from negatives
  label_mem_bg    — 1 if in 60 membrane experimental proteins
                    0 if membrane protein (TMD>0) NOT in 882 exp set
                    NaN otherwise (excluded from analysis)
  label_mem_singlepass_bg — membrane positives vs single-pass membrane background
  label_mem_multipass_bg  — membrane positives vs multi-pass membrane background

Predictor score columns
------------------------
  PICNIC_score, PICNIC_GO_score   — PICNIC RF + AF2 (+ GO annotations)        [raw]
  PSAP_score                      — Random Forest, AA composition             [raw]
  PSPHunter_prob                  — Ensemble ML, multi-DB trained             [raw]
  PSPire_score                    — XGBoost + AF2 structure features          [raw]
  FuzDrop_pLLPS                   — Biophysical theory (FuzPred + Espritz)    [raw]
  catGRANULE_score                — Physicochemical disorder/RNA-binding      [raw]
  PLAAC_NLLR                      — Prion-like HMM log-likelihood ratio       [raw]
  PScore_score                    — π–π contact frequency score               [raw]
  ESpritz_score                   — Disorder prediction (BRNN)                [raw]
  SEG_score                       — Low-complexity / sequence entropy         [raw]
  DeepPhase_score                 — ML with biophysical features (task 1: LLPS vs non-LLPS) [raw]
  SaPS_score                      — PhaSePred scaffold score (10-feature GB)  [raw]
  PdPS_score                      — PhaSePred driver score (10-feature GB)    [raw]
  PDL_score                       — ProtT5 + KmerConv prion-like domain model [percentile-
                                     normalised; only form available, but monotonic so
                                     AUROC/MaxMCC are unaffected]
  RY_score                        — Arg+Tyr sequence fraction (rule-based baseline) [raw]
  ParSe2_score                    — ParSe v2 PS potential, polymer-scaling/helix theory [raw]
  LLPhyScore_score                — 8 biophysical-interaction feature sum (pretrained) [raw]

Universe: PICNIC whole-genome file (~20,447 human proteins).
TMD annotation: full_dataset.csv (20,366 proteins, TMD_count column).

Output:
  output/background_scored.csv
"""

from pathlib import Path
import json

import pandas as pd

ROOT = Path(__file__).parent
PRED = ROOT / "Predictors_whole_genome_sets"
OUT  = ROOT / "output"

# ── Universe: PICNIC whole-genome file ───────────────────────────────────────

picnic = pd.read_csv(PRED / "PICNIC-9606-data.csv")
picnic = picnic.groupby("Uniprot ID", sort=False)[["PICNIC score", "PICNIC GO score"]].max().reset_index()
picnic = picnic.rename(columns={
    "Uniprot ID":      "UniProt_ID",
    "PICNIC score":    "PICNIC_score",
    "PICNIC GO score": "PICNIC_GO_score",
})
print(f"Universe (PICNIC): {len(picnic)} proteins")

df = picnic.copy()

# ── PSAP ──────────────────────────────────────────────────────────────────────

print("Loading PSAP...")
psap = pd.read_excel(PRED / "PSAP_full_proteome_excluding_training_set.xlsx")
psap = psap.drop_duplicates("Uniprot.ID").rename(columns={"Uniprot.ID": "UniProt_ID"})
df = df.merge(psap[["UniProt_ID", "PSAP_score"]], on="UniProt_ID", how="left")
print(f"  PSAP matched: {df['PSAP_score'].notna().sum()}")

# ── PSPHunter ─────────────────────────────────────────────────────────────────

print("Loading PSPHunter...")
psphunter = pd.read_excel(PRED / "PSPHunter_full_proteome.xlsx", skiprows=1)
psphunter = psphunter.rename(columns={"Uniprot": "UniProt_ID", "Probability": "PSPHunter_prob"})
df = df.merge(psphunter[["UniProt_ID", "PSPHunter_prob"]], on="UniProt_ID", how="left")
print(f"  PSPHunter matched: {df['PSPHunter_prob'].notna().sum()}")

# ── PSPire ────────────────────────────────────────────────────────────────────

print("Loading PSPire...")
pspire = pd.read_csv(PRED / "PSPire_Homo_sapiens_phos_scores.csv")
pspire = pspire.rename(columns={"Uniprot_ID": "UniProt_ID", "Score": "PSPire_score"})
df = df.merge(pspire[["UniProt_ID", "PSPire_score"]], on="UniProt_ID", how="left")
print(f"  PSPire matched: {df['PSPire_score'].notna().sum()}")

# ── FuzDrop ───────────────────────────────────────────────────────────────────

print("Loading FuzDrop...")
fuzdrop = pd.read_excel(PRED / "FuzDrop_ful_proteome.xlsx", sheet_name="Homo sapiens")
fuzdrop = fuzdrop.rename(columns={"Entry": "UniProt_ID", "p(LLPS)": "FuzDrop_pLLPS"})
df = df.merge(fuzdrop[["UniProt_ID", "FuzDrop_pLLPS"]], on="UniProt_ID", how="left")
print(f"  FuzDrop matched: {df['FuzDrop_pLLPS'].notna().sum()}")

# ── PhaSePred JSON — raw scores for 8 predictors ─────────────────────────────
# Keys: UniProt_ID → dict with nested predictor entries.
# catGRANULE['single'], PLAAC['NLLR'], PScore['single'],
# ESpritz-DisProt['single'], SEG['single'], DeepPhase['single'],
# PhaSePred['SaPS-10fea'], PhaSePred['PdPS-10fea']

print("Loading PhaSePred JSON (737 MB, this may take a moment)...")
with open(PRED / "PhaSePred_human_reviewed.json") as f:
    phasepred_raw = json.load(f)

def _extract_phasepred(raw):
    rows = []
    for uid, entry in raw.items():
        row = {"UniProt_ID": uid}
        try:
            row["catGRANULE_score"] = entry["catGRANULE"]["single"]
        except (KeyError, TypeError):
            row["catGRANULE_score"] = None
        try:
            row["PLAAC_NLLR"] = entry["PLAAC"]["NLLR"]
        except (KeyError, TypeError):
            row["PLAAC_NLLR"] = None
        try:
            row["PScore_score"] = entry["PScore"]["single"]
        except (KeyError, TypeError):
            row["PScore_score"] = None
        try:
            row["ESpritz_score"] = entry["ESpritz-DisProt"]["single"]
        except (KeyError, TypeError):
            row["ESpritz_score"] = None
        try:
            row["SEG_score"] = entry["SEG"]["single"]
        except (KeyError, TypeError):
            row["SEG_score"] = None
        # DeepPhase.single from PhaSePred JSON has incomplete coverage (~11k/20k);
        # prefer the dedicated DeePhase CSV which covers 20k proteins.
        # Slot is left empty here and merged separately below.
        try:
            row["SaPS_score"] = entry["PhaSePred"]["SaPS-10fea"]
        except (KeyError, TypeError):
            row["SaPS_score"] = None
        try:
            row["PdPS_score"] = entry["PhaSePred"]["PdPS-10fea"]
        except (KeyError, TypeError):
            row["PdPS_score"] = None
        rows.append(row)
    return pd.DataFrame(rows)

phasepred = _extract_phasepred(phasepred_raw)
del phasepred_raw  # free ~737 MB

phasepred_cols = ["catGRANULE_score", "PLAAC_NLLR", "PScore_score",
                  "ESpritz_score", "SEG_score", "SaPS_score", "PdPS_score"]
df = df.merge(phasepred, on="UniProt_ID", how="left")
print(f"  PhaSePred matched: {df['catGRANULE_score'].notna().sum()} (catGRANULE)")

# ── DeePhase (dedicated CSV — better coverage than PhaSePred JSON) ────────────

print("Loading DeePhase...")
deephase = pd.read_csv(PRED / "DeePhase_Swissprot_all_predictions.csv")
deephase = deephase.rename(columns={
    "Uniprot_ID":         "UniProt_ID",
    "prediction_phys_1":  "DeepPhase_score",
})
df = df.merge(deephase[["UniProt_ID", "DeepPhase_score"]], on="UniProt_ID", how="left")
print(f"  DeePhase matched: {df['DeepPhase_score'].notna().sum()}")

# ── PDL (percentile-normalised; raw probability unavailable) ─────────────────

print("Loading PDL...")
pdl = pd.read_excel(PRED / "PSPspredict_full_proteome.xlsx", sheet_name="PredictAll")
pdl = pdl.rename(columns={"Uniprot": "UniProt_ID", "PDL": "PDL_score"})
df = df.merge(pdl[["UniProt_ID", "PDL_score"]], on="UniProt_ID", how="left")
print(f"  PDL matched: {df['PDL_score'].notna().sum()}")

# ── R+Y composition (rule-based baseline, no training) ───────────────────────
# Arg+Tyr fraction — pi-pi/cation-pi residues implicated in LLPS scaffolding.

print("Computing R+Y composition...")

def _parse_fasta(path):
    seqs = {}
    uid, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if uid is not None:
                    seqs[uid] = "".join(chunks)
                parts = line[1:].split("|")
                uid = parts[1].strip() if len(parts) >= 2 else line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if uid is not None:
        seqs[uid] = "".join(chunks)
    return seqs

proteome_seqs = _parse_fasta(PRED / "Seq2Phase_swiss_prot_human_220916.fasta")
ry = pd.DataFrame([
    {"UniProt_ID": uid, "RY_score": (seq.count("R") + seq.count("Y")) / len(seq)}
    for uid, seq in proteome_seqs.items() if len(seq) > 0
])
df = df.merge(ry, on="UniProt_ID", how="left")
print(f"  R+Y matched: {df['RY_score'].notna().sum()}")

# ── ParSe v2 (rule-based, ported from the web tool's own JS — see score_parse2.py) ──

print("Loading ParSe v2...")
parse2 = pd.read_csv(OUT / "parse2_scores.csv")
df = df.merge(parse2, on="UniProt_ID", how="left")
print(f"  ParSe v2 matched: {df['ParSe2_score'].notna().sum()}")

# ── LLPhyScore (pretrained, run via their standalone package — see run_llphyscore.sh) ──

print("Loading LLPhyScore...")
llphyscore = pd.read_csv(OUT / "llphyscore_scores.csv")
df = df.merge(llphyscore, on="UniProt_ID", how="left")
print(f"  LLPhyScore matched: {df['LLPhyScore_score'].notna().sum()}")

# ── Coerce all score columns to numeric ───────────────────────────────────────

score_cols = [
    "PICNIC_score", "PICNIC_GO_score", "PSAP_score", "PSPHunter_prob",
    "PSPire_score", "FuzDrop_pLLPS",
    "catGRANULE_score", "PLAAC_NLLR", "PScore_score", "ESpritz_score",
    "SEG_score", "DeepPhase_score", "SaPS_score", "PdPS_score", "PDL_score",
    "RY_score", "ParSe2_score", "LLPhyScore_score",
]
for c in score_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

print(f"\nScored universe: {len(df)} proteins, {len(df.columns)} columns")

# ── TMD annotation ────────────────────────────────────────────────────────────

fd = pd.read_csv(ROOT / "output" / "full_dataset.csv")[["Entry", "TMD_count"]]
fd = fd.rename(columns={"Entry": "UniProt_ID"})
df = df.merge(fd, on="UniProt_ID", how="left")
df["TMD_count"] = pd.to_numeric(df["TMD_count"], errors="coerce").fillna(0).astype(int)
df["is_membrane_protein"] = df["TMD_count"] > 0
print(f"Membrane proteins in universe: {df['is_membrane_protein'].sum()}")

# ── Labels ────────────────────────────────────────────────────────────────────

exp_ids = set(pd.read_csv(ROOT / "output" / "exp_db_combined.csv")["UniProt_ID"])
mem_ids = set(pd.read_csv(ROOT / "exp_db" / "membrane_exp_db_matches.csv")["Entry"])

df["label_all"] = df["UniProt_ID"].isin(exp_ids).astype(int)
df["label_mem"] = df["UniProt_ID"].isin(mem_ids).astype(int)

is_mem_pos = df["UniProt_ID"].isin(mem_ids)
is_not_exp = ~df["UniProt_ID"].isin(exp_ids)
is_mem_bg  = df["is_membrane_protein"] & is_not_exp
is_single  = (df["TMD_count"] == 1) & is_not_exp
is_multi   = (df["TMD_count"] >  1) & is_not_exp

df["label_mem_bg"] = pd.array([pd.NA] * len(df), dtype="Int8")
df.loc[is_mem_pos, "label_mem_bg"] = 1
df.loc[is_mem_bg,  "label_mem_bg"] = 0

df["label_mem_fixed"] = pd.array([pd.NA] * len(df), dtype="Int8")
df.loc[is_mem_pos,               "label_mem_fixed"] = 1
df.loc[is_not_exp & ~is_mem_pos, "label_mem_fixed"] = 0

df["label_mem_singlepass_bg"] = pd.array([pd.NA] * len(df), dtype="Int8")
df.loc[is_mem_pos, "label_mem_singlepass_bg"] = 1
df.loc[is_single,  "label_mem_singlepass_bg"] = 0

df["label_mem_multipass_bg"] = pd.array([pd.NA] * len(df), dtype="Int8")
df.loc[is_mem_pos, "label_mem_multipass_bg"] = 1
df.loc[is_multi,   "label_mem_multipass_bg"] = 0

def _counts(col):
    v = df[col].dropna()
    return int((v == 1).sum()), int((v == 0).sum())

print("\nLabel counts:")
for col, desc in [
    ("label_all",                "all experimental LLPS vs proteome"),
    ("label_mem",                "membrane LLPS vs proteome (contaminated)"),
    ("label_mem_fixed",          "membrane LLPS vs proteome (fixed)"),
    ("label_mem_bg",             "membrane LLPS vs membrane background"),
    ("label_mem_singlepass_bg",  "membrane LLPS vs single-pass background"),
    ("label_mem_multipass_bg",   "membrane LLPS vs multi-pass background"),
]:
    p, n = _counts(col)
    print(f"  {col:<26} {desc:<42} pos={p:>4}  neg={n:>6}")

print("\nScore coverage (non-null in full universe):")
for c in score_cols:
    print(f"  {c:<22} {df[c].notna().sum():>6}")

df.to_csv(OUT / "background_scored.csv", index=False)
print(f"\nSaved output/background_scored.csv ({len(df)} rows × {len(df.columns)} cols)")
