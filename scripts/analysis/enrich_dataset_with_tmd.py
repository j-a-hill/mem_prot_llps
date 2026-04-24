"""
Enrich the protein dataset with transmembrane domain counts from UniProt.

Fetches the Transmembrane and Intramembrane annotations for every protein
in the source XLSX, counts domain spans, and writes an updated dataset.

Usage
-----
    python scripts/analysis/enrich_dataset_with_tmd.py

Outputs
-------
    data/uniprot_tm_cache.csv      raw TM annotations (reused on reruns)
    results/full_dataset.csv       updated dataset with TMD_count column
    dashboard/full_dataset.csv     copy for the Shinylive dashboard
"""

from pathlib import Path

from llps import add_tmd_count, fetch_uniprot_tm_annotations, load_llps_data

ROOT = Path(__file__).resolve().parents[2]

df = load_llps_data(str(ROOT / "Human Phase separation data.xlsx"))
entry_ids = df["Entry"].dropna().tolist()

tm_df = fetch_uniprot_tm_annotations(
    entry_ids,
    cache_path=str(ROOT / "data" / "uniprot_tm_cache.csv"),
)

df = add_tmd_count(df, tm_annotations=tm_df)

for dest in (ROOT / "results" / "full_dataset.csv", ROOT / "dashboard" / "full_dataset.csv"):
    df.to_csv(dest, index=False)
    print(f"Saved: {dest}")
