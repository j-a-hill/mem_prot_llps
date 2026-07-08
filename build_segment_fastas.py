"""
build_segment_fastas.py

Generates per-segment FASTA files for the 60 membrane LLPS proteins.
Each contiguous annotated topology span becomes a separate FASTA entry,
unlike the concatenated-region files in output/topology_fasta/ which join
all spans of the same type into one chimeric sequence.

Outputs
-------
  output/topology_fasta/segments/Cytoplasmic_segments.fasta
  output/topology_fasta/segments/Transmembrane_segments.fasta
  output/topology_fasta/segments/Extracellular_Lumenal_segments.fasta
  output/topology_fasta/segments/segment_manifest.csv

Header format: >{uid}|{region_slug}|seg{i}
  region_slug: cyto | tm | extra
  i: 0-indexed segment number within that protein
"""

import re
import statistics
from pathlib import Path

import pandas as pd

# ── Paths ────────────────────────────────────────────────────────────────────

ROOT        = Path(__file__).parent
PRED        = ROOT / "Predictors_whole_genome_sets"
OUT         = ROOT / "output"
SEG_DIR     = OUT / "topology_fasta" / "segments"
TOPO_CACHE  = ROOT / "data" / "uniprot_topology_cache.csv"
EXP_DB      = ROOT / "exp_db" / "membrane_exp_db_matches.csv"
FASTA_SRC   = PRED / "Seq2Phase_swiss_prot_human_220916.fasta"

SEG_DIR.mkdir(parents=True, exist_ok=True)

# Region slug mapping
REGION_SLUG = {
    "Cytoplasmic":        "cyto",
    "Transmembrane":      "tm",
    "Extracellular/Lumenal": "extra",
}
REGION_FILE = {
    "Cytoplasmic":           SEG_DIR / "Cytoplasmic_segments.fasta",
    "Transmembrane":         SEG_DIR / "Transmembrane_segments.fasta",
    "Extracellular/Lumenal": SEG_DIR / "Extracellular_Lumenal_segments.fasta",
}

# ── Load positive protein IDs ─────────────────────────────────────────────────

pos_ids = pd.read_csv(EXP_DB)["Entry"].tolist()
print(f"Loaded {len(pos_ids)} positive protein IDs from {EXP_DB.relative_to(ROOT)}")

# ── Load sequences ────────────────────────────────────────────────────────────

seqs: dict[str, str] = {}
uid: str | None = None
chunks: list[str] = []
with open(FASTA_SRC) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if uid is not None:
                seqs[uid] = "".join(chunks)
            parts = line[1:].split("|")
            uid = parts[1].strip() if len(parts) >= 2 else None
            chunks = []
        else:
            chunks.append(line)
    if uid is not None:
        seqs[uid] = "".join(chunks)

print(f"Loaded {len(seqs)} sequences from {FASTA_SRC.name}")
found = sum(1 for u in pos_ids if u in seqs)
print(f"  {found}/{len(pos_ids)} positive IDs have a sequence")

# ── Load topology cache ────────────────────────────────────────────────────────

if not TOPO_CACHE.exists():
    raise SystemExit(f"{TOPO_CACHE} not found — run extract_topology_scores.py first.")

topo_raw = (
    pd.read_csv(TOPO_CACHE)
    .rename(columns={"Entry": "UniProt_ID",
                     "Topological domain": "topo_dom_raw",
                     "Transmembrane": "transmem_raw"})
    .set_index("UniProt_ID")
)

# ── Feature parsing (same regex as extract_topology_fasta.py) ─────────────────

FEATURE_RE = re.compile(r'(TOPO_DOM|TRANSMEM)\s+(\d+)\.\.(\d+);\s*/note="([^"]*)"')


def parse_features(raw_text) -> list[tuple[str, int, int, str]]:
    """Return list of (feature_type, start_1based, end_1based, note)."""
    if pd.isna(raw_text):
        return []
    return [
        (m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))
        for m in FEATURE_RE.finditer(str(raw_text))
    ]


def bucket_topo_dom(note: str) -> str:
    return "Cytoplasmic" if "ytoplasmic" in note else "Extracellular/Lumenal"


def get_segments(uid: str, seq: str) -> list[tuple[str, int, int, str]]:
    """
    Return ordered list of (region, start_1based, end_1based, subsequence)
    for every annotated span.  Spans are returned in genomic order.
    """
    if uid not in topo_raw.index:
        return []

    row = topo_raw.loc[uid]
    spans: list[tuple[str, int, int]] = []

    for _, s, e, note in parse_features(row["topo_dom_raw"]):
        e_clipped = min(e, len(seq))
        if s <= e_clipped:
            spans.append((bucket_topo_dom(note), s, e_clipped))

    for _, s, e, note in parse_features(row["transmem_raw"]):
        e_clipped = min(e, len(seq))
        if s <= e_clipped:
            spans.append(("Transmembrane", s, e_clipped))

    # Sort by start position
    spans.sort(key=lambda x: x[1])

    return [(region, s, e, seq[s - 1:e]) for region, s, e in spans]


# ── Build per-region segment lists ────────────────────────────────────────────

# region -> list of (uid, seg_idx_within_protein, start, end, seq)
region_segments: dict[str, list[tuple[str, int, int, int, str]]] = {
    r: [] for r in REGION_SLUG
}

manifest_rows: list[dict] = []

for uid in pos_ids:
    seq = seqs.get(uid)
    if seq is None:
        continue

    all_spans = get_segments(uid, seq)

    # Track per-region segment index for this protein
    region_counter: dict[str, int] = {r: 0 for r in REGION_SLUG}

    for region, start, end, subseq in all_spans:
        if region not in REGION_SLUG:
            continue  # shouldn't happen, but be safe

        seg_idx = region_counter[region]
        region_counter[region] += 1

        region_segments[region].append((uid, seg_idx, start, end, subseq))

        manifest_rows.append({
            "UniProt_ID":  uid,
            "region":      region,
            "seg_idx":     seg_idx,
            "start_1based": start,
            "end_1based":   end,
            "length":       len(subseq),
            "sequence":     subseq,
            "short":        len(subseq) < 5,
        })

# ── Write FASTA files ─────────────────────────────────────────────────────────

for region, fpath in REGION_FILE.items():
    slug = REGION_SLUG[region]
    with open(fpath, "w") as f:
        for uid, seg_idx, _start, _end, subseq in region_segments[region]:
            header = f">{uid}|{slug}|seg{seg_idx}"
            f.write(f"{header}\n{subseq}\n")
    n_segs = len(region_segments[region])
    n_prots = len({uid for uid, *_ in region_segments[region]})
    print(f"Wrote {fpath.relative_to(ROOT)}  "
          f"({n_segs} segments from {n_prots} proteins)")

# ── Write segment manifest ────────────────────────────────────────────────────

manifest_df = pd.DataFrame(manifest_rows, columns=[
    "UniProt_ID", "region", "seg_idx", "start_1based", "end_1based",
    "length", "sequence", "short",
])
manifest_path = SEG_DIR / "segment_manifest.csv"
manifest_df.to_csv(manifest_path, index=False)
print(f"\nWrote {manifest_path.relative_to(ROOT)}  ({len(manifest_df)} rows)")

# ── Summary statistics ────────────────────────────────────────────────────────

print("\n── Summary statistics ──────────────────────────────────────────────────")
for region in REGION_SLUG:
    segs = region_segments[region]
    if not segs:
        print(f"\n{region}:  no segments found")
        continue

    prots_with_segs = {uid for uid, *_ in segs}
    lengths = [len(subseq) for _, _, _, _, subseq in segs]
    short_count = sum(1 for l in lengths if l < 5)

    print(f"\n{region}:")
    print(f"  proteins with ≥1 segment : {len(prots_with_segs)}")
    print(f"  total segments           : {len(segs)}")
    print(f"  min / median / max length: "
          f"{min(lengths)} / {statistics.median(lengths):.0f} / {max(lengths)}")
    if short_count:
        print(f"  short segments (<5 aa)   : {short_count}")

print("\nDone.")
