"""
build_topology_scores_master.py

Merges all topology scoring approaches into a single long-form table.

Output: output/topology_scores_master.csv
  Columns: UniProt_ID, tool, approach, region, seg_idx, score

Approach labels:
  region_mean  - per-residue mean within topology region (full-protein context)
                 only for per-residue tools: catGRANULE, PLAAC, PScore, ESpritz, SEG, ParSe2
  whole_db     - whole-protein score from precomputed database
  whole        - whole-protein score from fresh local run
  concatenated - score for all region residues joined into one sequence
  segment      - score for each individual contiguous span
"""

import pandas as pd
from pathlib import Path

ROOT = Path(__file__).parent
OUT  = ROOT / "output"
RAW  = OUT / "topology_scores_raw"
TOPO = OUT / "topology_domain_scores.csv"

REGION_MEAN_TOOLS = ["catGRANULE", "PLAAC", "PScore", "ESpritz", "SEG", "ParSe2"]
REGIONS = ["Cytoplasmic", "Transmembrane", "Extracellular/Lumenal"]

# ── 1. region_mean + whole_db from topology_domain_scores.csv ──────────────
topo = pd.read_csv(TOPO)
rm_rows = []

for _, row in topo.iterrows():
    uid = row["UniProt_ID"]
    for tool in REGION_MEAN_TOOLS:
        whole_col = f"{tool}_whole_official"
        if whole_col in row.index and pd.notna(row[whole_col]):
            rm_rows.append({"UniProt_ID": uid, "tool": tool,
                            "approach": "whole_db", "region": "Whole",
                            "seg_idx": None, "score": row[whole_col]})
        for region in REGIONS:
            col = f"{tool}_{region}"
            if col in row.index and pd.notna(row[col]):
                rm_rows.append({"UniProt_ID": uid, "tool": tool,
                                "approach": "region_mean", "region": region,
                                "seg_idx": None, "score": row[col]})

rm_df = pd.DataFrame(rm_rows)
print(f"region_mean/whole_db rows (per-residue tools): {len(rm_df)}")

# ── 1b. Concat/segment rows derived from full-protein per-residue arrays ──────
per_res_cs_path = RAW / "per_residue_tools_concat_segment.csv"
per_res_cs = pd.read_csv(per_res_cs_path)
print(f"per-residue concat/segment rows: {len(per_res_cs)}")
rm_df = pd.concat([rm_df, per_res_cs], ignore_index=True)

# ── 2. all_approaches.csv files ────────────────────────────────────────────
ALL_TOOLS = [
    "ParSe2", "LLPhyScore", "PSPsPredict", "PSAP",
    "FuzDrop", "DeePhase", "PSPHunter",
    "PSPire", "PICNIC",
]

aa_dfs = []
for tool in ALL_TOOLS:
    fpath = RAW / f"{tool}_all_approaches.csv"
    if not fpath.exists():
        print(f"  [skip] {tool} — file not found")
        continue
    df = pd.read_csv(fpath)
    df["tool"] = tool

    # ── normalise score column name ──
    if "score" not in df.columns:
        if "score_mean" in df.columns:
            # ParSe2: per-residue mean (length-normalised, better for region comparison)
            df["score"] = df["score_mean"]
        elif "score_per_residue" in df.columns:
            # LLPhyScore: length-normalised per-residue
            df["score"] = df["score_per_residue"]
        elif "raw_score" in df.columns:
            df["score"] = df["raw_score"]
        else:
            print(f"  [warn] {tool} — cannot find score column in {list(df.columns)}")
            continue

    # ── normalise PSPHunter approach/region labels ──
    if tool == "PSPHunter":
        # Whole.fasta was run with approach='concat' and region='Whole'
        df.loc[(df["approach"] == "concat") & (df["region"] == "Whole"), "approach"] = "whole"
        df["approach"] = df["approach"].replace({"concat": "concatenated", "segments": "segment"})
        df["region"]   = df["region"].replace({"Extracellular_Lumenal": "Extracellular/Lumenal"})

    aa_dfs.append(df[["UniProt_ID", "tool", "approach", "region", "seg_idx", "score"]])
    counts = df.groupby("approach")["UniProt_ID"].count().to_dict()
    print(f"  {tool}: {len(df)} rows — {counts}")

aa_df = pd.concat(aa_dfs, ignore_index=True) if aa_dfs else pd.DataFrame()

# ── 3. Merge and save ──────────────────────────────────────────────────────
master = pd.concat([rm_df, aa_df], ignore_index=True)
master["seg_idx"] = pd.to_numeric(master["seg_idx"], errors="coerce")

print(f"\nTotal rows: {len(master)}")
print("\nCoverage by tool × approach:")
pivot = master.groupby(["tool", "approach"])["UniProt_ID"].count().unstack(fill_value=0)
print(pivot.to_string())

out_path = OUT / "topology_scores_master.csv"
master.to_csv(out_path, index=False)
print(f"\nSaved → {out_path.relative_to(ROOT)}")
