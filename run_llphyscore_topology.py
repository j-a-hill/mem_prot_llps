"""
Run LLPhyScore on per-segment and concatenated topology FASTA files.
Saves raw CSVs and builds a merged summary CSV.
"""
import subprocess
import shutil
import sys
from pathlib import Path

import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO = Path("/home/jake/mem_prot_llps/mem_prot_llps")
LLPHY_DIR = REPO / "external_tools/LLPhyScore/standalone_package/LLPhyScore"
FASTA_DIR = REPO / "output/topology_fasta"
SEG_DIR   = FASTA_DIR / "segments"
OUT_DIR   = REPO / "output/topology_scores_raw"
SCRATCHPAD = Path(
    "/tmp/claude-1000/-home-jake-mem-prot-llps-mem-prot-llps/"
    "edb83b3c-8ad9-4357-af68-f6bcf9096c19/scratchpad"
)

OUT_DIR.mkdir(parents=True, exist_ok=True)
SCRATCHPAD.mkdir(parents=True, exist_ok=True)

PYTHON = str(REPO / ".venv/bin/python3")

# ── FASTA file registry ────────────────────────────────────────────────────────
# (input_fasta, output_csv_name, approach, region)
FILES = [
    # segments
    (SEG_DIR / "Cytoplasmic_segments.fasta",
     "LLPhyScore_Cytoplasmic_segments_raw.csv",
     "segment", "Cytoplasmic"),
    (SEG_DIR / "Transmembrane_segments.fasta",
     "LLPhyScore_Transmembrane_segments_raw.csv",
     "segment", "Transmembrane"),
    (SEG_DIR / "Extracellular_Lumenal_segments.fasta",
     "LLPhyScore_Extracellular_Lumenal_segments_raw.csv",
     "segment", "Extracellular/Lumenal"),
    # concatenated
    (FASTA_DIR / "Cytoplasmic.fasta",
     "LLPhyScore_Cytoplasmic_concat_raw.csv",
     "concatenated", "Cytoplasmic"),
    (FASTA_DIR / "Transmembrane.fasta",
     "LLPhyScore_Transmembrane_concat_raw.csv",
     "concatenated", "Transmembrane"),
    (FASTA_DIR / "Extracellular_Lumenal.fasta",
     "LLPhyScore_Extracellular_Lumenal_concat_raw.csv",
     "concatenated", "Extracellular/Lumenal"),
    # whole
    (FASTA_DIR / "Whole.fasta",
     "LLPhyScore_Whole_raw.csv",
     "whole", "Whole"),
]

SLUG_TO_REGION = {
    "cyto":  "Cytoplasmic",
    "tm":    "Transmembrane",
    "extra": "Extracellular/Lumenal",
}


# ── Helper: parse sequence lengths from a FASTA ───────────────────────────────
def parse_fasta_lengths(fasta_path: Path) -> dict[str, int]:
    """Return {header: length} for every sequence in the FASTA."""
    lengths = {}
    current_tag = None
    current_len = 0
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_tag is not None:
                    lengths[current_tag] = current_len
                current_tag = line[1:]  # strip leading '>'
                current_len = 0
            else:
                current_len += len(line)
        if current_tag is not None:
            lengths[current_tag] = current_len
    return lengths


# ── Helper: run LLPhyScore on one file ────────────────────────────────────────
def run_llphyscore(fasta_path: Path, out_csv: Path) -> pd.DataFrame | None:
    """Run LLPhyScore and return the resulting DataFrame (or None on failure)."""
    tmp_out = SCRATCHPAD / out_csv.name
    cmd = [
        PYTHON, "LLPhyScore_standalone.py",
        "-i", str(fasta_path),
        "-s", "raw",
        "-o", str(tmp_out),
    ]
    print(f"\n{'='*60}")
    print(f"Running LLPhyScore on: {fasta_path.name}")
    result = subprocess.run(
        cmd,
        cwd=str(LLPHY_DIR),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"  ERROR (return code {result.returncode}):")
        print(result.stderr[-2000:])
        return None
    if not tmp_out.exists():
        print(f"  ERROR: output file not created at {tmp_out}")
        return None
    df = pd.read_csv(tmp_out)
    # copy to permanent location
    shutil.copy2(tmp_out, out_csv)
    print(f"  Scored {len(df)} sequences → {out_csv.name}")
    return df


# ── Helper: parse UniProt ID from tag ─────────────────────────────────────────
def uid_from_tag(tag: str, approach: str) -> str:
    """Extract UniProt accession from FASTA header tag."""
    if approach == "segment":
        # format: {uid}|{region_slug}|seg{i}
        return tag.split("|")[0]
    else:
        # plain UniProt ID (no pipes), but guard against sp|UID|NAME format
        parts = tag.split("|")
        if len(parts) >= 2:
            return parts[1]
        return parts[0]


# ── Helper: parse seg_idx from tag ────────────────────────────────────────────
def seg_idx_from_tag(tag: str) -> int | float:
    """Extract segment index from segment header, or NaN."""
    try:
        # format: {uid}|{region_slug}|seg{i}
        last = tag.split("|")[-1]
        return int(last.replace("seg", ""))
    except (ValueError, IndexError):
        return float("nan")


# ── Main ──────────────────────────────────────────────────────────────────────
all_rows = []

for fasta_path, out_name, approach, region in FILES:
    out_csv = OUT_DIR / out_name
    if not fasta_path.exists():
        print(f"  SKIPPING (file not found): {fasta_path}")
        continue

    df_raw = run_llphyscore(fasta_path, out_csv)
    if df_raw is None:
        continue

    # parse lengths from FASTA
    lengths = parse_fasta_lengths(fasta_path)

    for _, row in df_raw.iterrows():
        tag = row["tag"]
        raw_score = row["8-feature sum"]
        length = lengths.get(tag, float("nan"))
        uid = uid_from_tag(tag, approach)

        if approach == "segment":
            parts = tag.split("|")
            slug = parts[1] if len(parts) > 1 else ""
            actual_region = SLUG_TO_REGION.get(slug, region)
            s_idx = seg_idx_from_tag(tag)
        else:
            actual_region = region
            s_idx = float("nan")

        score_per_res = raw_score / length if length and length > 0 else float("nan")

        all_rows.append({
            "UniProt_ID":       uid,
            "approach":         approach,
            "region":           actual_region,
            "seg_idx":          s_idx,
            "raw_score":        raw_score,
            "length":           length,
            "score_per_residue": score_per_res,
            "tag":              tag,  # keep original tag for debugging
        })

# ── Build merged summary ───────────────────────────────────────────────────────
summary = pd.DataFrame(all_rows)
summary_path = OUT_DIR / "LLPhyScore_all_approaches.csv"
summary.to_csv(summary_path, index=False)

print(f"\n{'='*60}")
print(f"Merged summary written: {summary_path}")
print(f"Total rows: {len(summary)}")
print(f"\nBreakdown by approach × region:")
print(summary.groupby(["approach", "region"]).size().to_string())
print("\nDone.")
