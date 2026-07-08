"""
export_membrane_scores_table.py

Exports the 60 experimentally-confirmed membrane LLPS proteins with every
predictor's raw score, plus per-tool training-set leakage flags and the
predictor category legend, to a single XLSX for sharing outside the repo.

Output: output/membrane_scores_table.xlsx
  Sheet "Scores"   — UniProt_ID, gene/protein metadata, 18 raw predictor scores
  Sheet "Leakage"  — clean_<col> flags (True = not in that tool's own training set)
  Sheet "Categories" — predictor -> target type -> mechanism legend
"""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parent
OUT  = ROOT / "output"

# Predictor order/columns and categorisation, kept in sync with plot_roc.py.
PREDICTORS = {
    "PICNIC":      "PICNIC_score",
    "PICNIC (GO)": "PICNIC_GO_score",
    "PSPire":      "PSPire_score",
    "PSPHunter":   "PSPHunter_prob",
    "SaPS":        "SaPS_score",
    "PdPS":        "PdPS_score",
    "PDL":         "PDL_score",
    "LLPhyScore":  "LLPhyScore_score",
    "PScore":      "PScore_score",
    "R+Y":         "RY_score",
    "ParSe2":      "ParSe2_score",
    "FuzDrop":     "FuzDrop_pLLPS",
    "PSAP":        "PSAP_score",
    "DeepPhase":   "DeepPhase_score",
    "catGRANULE":  "catGRANULE_score",
    "PLAAC":       "PLAAC_NLLR",
    "ESpritz":     "ESpritz_score",
    "SEG":         "SEG_score",
}

IDR_PROXY_TOOLS = {"PLAAC", "ESpritz", "SEG"}

MECHANISM = {
    "PICNIC": "AF2 structure + sequence", "PICNIC (GO)": "AF2 structure + sequence",
    "PSPire": "AF2 structure + sequence",
    "PSPHunter": "Multi-feature ML ensemble", "SaPS": "Multi-feature ML ensemble",
    "PdPS": "Multi-feature ML ensemble", "PDL": "Multi-feature ML ensemble",
    "LLPhyScore": "Multi-feature ML ensemble",
    "PScore": "Pi-pi / cation-pi", "R+Y": "Pi-pi / cation-pi",
    "ParSe2": "Sticker-spacer / polymer theory", "FuzDrop": "Sticker-spacer / polymer theory",
    "PSAP": "IDR / disorder", "DeepPhase": "IDR / disorder", "catGRANULE": "IDR / disorder",
    "PLAAC": "IDR / disorder", "ESpritz": "IDR / disorder", "SEG": "IDR / disorder",
}

NOTES = {
    "SaPS": "PhaSePred's own score (Chen et al., PNAS 2022), sourced from PhaSePred's raw "
            "JSON. Shares PhaSePred's model + catGRANULE/PLAAC/PScore/ESpritz/SEG/DeepPhase "
            "as input features — not independent of those tools.",
    "PdPS": "PhaSePred's own score (Chen et al., PNAS 2022), sourced from PhaSePred's raw "
            "JSON. Shares PhaSePred's model + catGRANULE/PLAAC/PScore/ESpritz/SEG/DeepPhase "
            "as input features — not independent of those tools.",
    "PDL":  "Wenbin Li's own tool (\"Protein Dual-model Language\", LLM+KmerConv, published "
            "as PSPsPredict on GitHub) — not PhaSePred. His benchmark compilation spreadsheet "
            "re-publishes PhaSePred's SaPS/PdPS (percentile-ranked) alongside his own PDL "
            "score, but our SaPS_score/PdPS_score come from PhaSePred's raw JSON, not that "
            "spreadsheet, so PDL has no numeric link to SaPS/PdPS in this table.",
    "LLPhyScore": "Negative leakage is checkable (LLPhyScore's 2000-protein generic-human-"
            "proteome training negatives are UniProt-keyed) — 7/60 membrane benchmark "
            "proteins are confirmed training negatives. The positive (true-positive) set is "
            "keyed by construct/gene name; gene-symbol matching found zero hits, but exact/"
            "fragment sequence matching against the human proteome resolved 66/305 (22%) of "
            "the training sequences to real UniProt IDs — 2/60 membrane benchmark proteins "
            "(O60500/Nephrin, P08908/HTR1A) are confirmed training positives this way. The "
            "remaining 239 unresolved sequences (fusions/non-human orthologs/synthetic "
            "repeats) stay genuinely unknown rather than falsely 'confirmed clean'.",
}

# ── Load data ──────────────────────────────────────────────────────────────────

master = pd.read_csv(OUT / "master_table.csv")
mem = master[master["is_membrane"]].copy()
print(f"Membrane benchmark: {len(mem)} proteins")

scored = pd.read_csv(OUT / "background_scored.csv")
mem = mem.merge(scored[["UniProt_ID"] + [c for c in PREDICTORS.values() if c in scored.columns]],
                 on="UniProt_ID", how="left")

if "ParSe2_score" not in mem.columns:
    parse2 = pd.read_csv(OUT / "parse2_scores.csv")
    mem = mem.merge(parse2, on="UniProt_ID", how="left")

if "LLPhyScore_score" not in mem.columns:
    mem["LLPhyScore_score"] = pd.NA

masks = pd.read_csv(OUT / "clean_masks.csv")
mem_leak = mem[["UniProt_ID", "Gene_Names"]].merge(masks, on="UniProt_ID", how="left")

# ── Sheet 1: Scores ────────────────────────────────────────────────────────────

meta_cols = ["UniProt_ID", "Gene_Names", "Protein_Name", "TMD_count", "TMD_class",
             "Length", "LLPS_mode", "role", "MLO_Types", "Source_DB"]
score_cols = list(PREDICTORS.values())
scores_sheet = mem[meta_cols + score_cols].rename(
    columns={v: f"{k} ({v})" for k, v in PREDICTORS.items()})

# ── Sheet 2: Leakage flags ─────────────────────────────────────────────────────

clean_cols = ["clean_global"] + [f"clean_{c}" for c in score_cols if f"clean_{c}" in mem_leak.columns]
leakage_sheet = mem_leak[["UniProt_ID", "Gene_Names"] + clean_cols]

# ── Sheet 3: Categories ────────────────────────────────────────────────────────

categories_sheet = pd.DataFrame([
    {
        "Predictor": p,
        "Score column": col,
        "Target type": "IDR proxy (not LLPS-trained)" if p in IDR_PROXY_TOOLS else "LLPS-specific",
        "Mechanism": MECHANISM[p],
        "Note": NOTES.get(p, ""),
    }
    for p, col in PREDICTORS.items()
])

# ── Write ──────────────────────────────────────────────────────────────────────

out_path = OUT / "membrane_scores_table.xlsx"
with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
    scores_sheet.to_excel(writer, sheet_name="Scores", index=False)
    leakage_sheet.to_excel(writer, sheet_name="Leakage", index=False)
    categories_sheet.to_excel(writer, sheet_name="Categories", index=False)

print(f"Saved {out_path}  ({len(scores_sheet)} proteins x {len(score_cols)} predictors)")
