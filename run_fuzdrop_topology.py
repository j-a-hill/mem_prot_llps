"""
run_fuzdrop_topology.py

Run FuzDrop on per-segment topology FASTAs and on concatenated region FASTAs
for the 60 membrane LLPS proteins.

Outputs:
  output/topology_scores_raw/FuzDrop_segments_raw.csv
    columns: UniProt_ID, region, seg_idx, length, fuzdrop_score

  output/topology_scores_raw/FuzDrop_concat_raw.csv
    columns: UniProt_ID, region, length, fuzdrop_score
"""

import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT           = Path(__file__).parent
PRED           = ROOT / "Predictors_whole_genome_sets"
OUT            = ROOT / "output"
RAW_OUT        = OUT / "topology_scores_raw"
FASTA_DIR      = OUT / "topology_fasta"
TOPO_CACHE     = ROOT / "data" / "uniprot_topology_cache.csv"
PHASEPRED_JSON = PRED / "PhaSePred_human_reviewed.json"
FUZDROP_BIN    = ROOT / "external_tools" / "FuzDrop" / "FuzDrop"
SEG_MANIFEST   = FASTA_DIR / "segments" / "segment_manifest.csv"

RAW_OUT.mkdir(parents=True, exist_ok=True)

# ── ESpritz header (copied verbatim from extract_topology_fasta.py) ──────────
ESPRITZ_HEADER = (
    "*" * 108 + "\n\n"
    "Licensed to: Prof. MONIKA FUXREITER (University of Padova) located in Padova, Italy. \n"
    "This license is for non-commercial use only. Please see LICENSE file for details \n"
    "(https://biocomputingup.it/assets/data/LICENSE.txt)\n\n"
    "Contact silvio.tosatto@unipd.it for commercial licensing details.\n\n\n"
    + "*" * 108 + "\n"
)

FEATURE_RE = re.compile(r'(TOPO_DOM|TRANSMEM)\s+(\d+)\.\.(\d+);\s*/note="([^"]*)"')

# ── Load topology cache for rebuild_region_labels (concat approach) ──────────
topo_raw = pd.read_csv(TOPO_CACHE).rename(columns={
    "Entry": "UniProt_ID",
    "Topological domain": "topo_dom_raw",
    "Transmembrane": "transmem_raw",
}).set_index("UniProt_ID")


def parse_features(raw_text):
    if pd.isna(raw_text):
        return []
    return [(m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))
            for m in FEATURE_RE.finditer(str(raw_text))]


def bucket_topo_dom(note):
    return "Cytoplasmic" if "ytoplasmic" in note else "Extracellular/Lumenal"


def build_region_labels(uid, length):
    labels = np.full(length, None, dtype=object)
    if uid not in topo_raw.index:
        return labels
    row = topo_raw.loc[uid]
    for _, s, e, note in parse_features(row["topo_dom_raw"]):
        labels[s - 1:min(e, length)] = bucket_topo_dom(note)
    for _, s, e, note in parse_features(row["transmem_raw"]):
        labels[s - 1:min(e, length)] = "Transmembrane"
    return labels


# ── Helper: run FuzDrop in a temp dir, return score or None ─────────────────

def run_fuzdrop(stem: str, fasta_content: str, espritz_content: str, work_dir: Path):
    """Write FASTA + espritz, run FuzDrop, parse p(LLPS), return float or None."""
    fasta_path   = work_dir / f"{stem}.fasta"
    espritz_path = work_dir / f"{stem}.espritz"
    res_path     = work_dir / f"{stem}_res.txt"

    fasta_path.write_text(fasta_content)
    espritz_path.write_text(espritz_content)

    try:
        result = subprocess.run(
            [str(FUZDROP_BIN.resolve()), fasta_path.name, espritz_path.name],
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=120,
        )
    except subprocess.TimeoutExpired:
        return None
    finally:
        fasta_path.unlink(missing_ok=True)
        espritz_path.unlink(missing_ok=True)

    score = None
    if result.returncode == 0 and res_path.exists():
        last_line = res_path.read_text().strip().splitlines()[-1]
        try:
            score = float(last_line.split("=")[-1].strip())
        except ValueError:
            pass
        res_path.unlink(missing_ok=True)

    return score


def make_espritz_content(labels_iter, scores_iter):
    body = "\n".join(f"{l}\t{s}" for l, s in zip(labels_iter, scores_iter))
    return ESPRITZ_HEADER + body + "\n"


# ── Load PhaSePred JSON once ─────────────────────────────────────────────────

print("Loading PhaSePred JSON (737 MB) — please wait...")
with open(PHASEPRED_JSON) as f:
    phasepred = json.load(f)
print(f"  Loaded {len(phasepred):,} entries.")

# ── Load segment manifest ────────────────────────────────────────────────────

manifest = pd.read_csv(SEG_MANIFEST)
print(f"\nSegment manifest: {len(manifest)} segments across "
      f"{manifest['UniProt_ID'].nunique()} proteins")

# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: Per-segment FuzDrop
# ═══════════════════════════════════════════════════════════════════════════════

work_dir = Path(tempfile.mkdtemp(prefix="fuzdrop_seg_"))
print(f"\n[Segments] Working dir: {work_dir}")

seg_rows = []
n_ok = n_skip = n_fail = 0

for i, row in manifest.iterrows():
    uid       = row["UniProt_ID"]
    region    = row["region"]
    seg_idx   = int(row["seg_idx"])
    start0    = int(row["start_1based"]) - 1   # 0-based inclusive
    end0      = int(row["end_1based"])          # 0-based exclusive  (= end_1based)
    length    = int(row["length"])
    sequence  = row["sequence"]
    is_short  = bool(row["short"])

    if (i % 20) == 0:
        print(f"  [{i}/{len(manifest)}] {uid} {region} seg{seg_idx}  "
              f"(len={length}{'  SHORT' if is_short else ''})")

    # Slice ESpritz data
    esp = phasepred.get(uid, {}).get("ESpritz-DisProt", {})
    esp_labels_all = esp.get("label", "").split(",")
    esp_scores_all = esp.get("residue", "").split(",")

    full_seq = phasepred.get(uid, {}).get("Sequence", "")
    if not full_seq or len(esp_labels_all) != len(full_seq) or len(esp_scores_all) != len(full_seq):
        n_skip += 1
        continue

    esp_labels_seg = esp_labels_all[start0:end0]
    esp_scores_seg = esp_scores_all[start0:end0]

    if len(esp_labels_seg) != length:
        # ESpritz length mismatch for this slice
        n_skip += 1
        continue

    region_slug = region.replace("/", "_")
    stem = f"{uid}_{region_slug}_seg{seg_idx}"

    fasta_content   = f">{uid}|{region_slug}|seg{seg_idx}\n{sequence}\n"
    espritz_content = make_espritz_content(esp_labels_seg, esp_scores_seg)

    score = run_fuzdrop(stem, fasta_content, espritz_content, work_dir)

    if score is not None:
        seg_rows.append({
            "UniProt_ID":    uid,
            "region":        region,
            "seg_idx":       seg_idx,
            "length":        length,
            "fuzdrop_score": score,
        })
        n_ok += 1
    else:
        if not is_short:
            print(f"    [FAIL] {uid} {region} seg{seg_idx} len={length}")
        n_fail += 1

shutil.rmtree(work_dir, ignore_errors=True)

seg_df = pd.DataFrame(seg_rows)
seg_out = RAW_OUT / "FuzDrop_segments_raw.csv"
seg_df.to_csv(seg_out, index=False)
print(f"\n[Segments] Done: {n_ok} scored, {n_fail} failed, {n_skip} skipped (no ESpritz data)")
print(f"  Saved {len(seg_df)} rows to {seg_out.relative_to(ROOT)}")

# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: Concat region FASTAs
# ═══════════════════════════════════════════════════════════════════════════════

REGION_TO_FILE = {
    "Cytoplasmic":         FASTA_DIR / "Cytoplasmic.fasta",
    "Transmembrane":       FASTA_DIR / "Transmembrane.fasta",
    "Extracellular/Lumenal": FASTA_DIR / "Extracellular_Lumenal.fasta",
}

def parse_fasta(fasta_path):
    """Return dict of {uid: sequence} from a FASTA file."""
    seqs = {}
    uid_cur, chunks = None, []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if uid_cur is not None:
                    seqs[uid_cur] = "".join(chunks)
                uid_cur = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if uid_cur is not None:
        seqs[uid_cur] = "".join(chunks)
    return seqs


work_dir = Path(tempfile.mkdtemp(prefix="fuzdrop_concat_"))
print(f"\n[Concat] Working dir: {work_dir}")

concat_rows = []
n_ok = n_skip = n_fail = 0
total_count = 0

for region, fasta_path in REGION_TO_FILE.items():
    if not fasta_path.exists():
        print(f"  [skip] {fasta_path.name} not found")
        continue

    region_seqs = parse_fasta(fasta_path)
    print(f"\n  {region}: {len(region_seqs)} sequences in {fasta_path.name}")

    for j, (uid, concat_seq) in enumerate(region_seqs.items()):
        if j % 20 == 0:
            print(f"    [{j}/{len(region_seqs)}] {uid} len={len(concat_seq)}")

        # Get full-protein ESpritz arrays
        esp = phasepred.get(uid, {}).get("ESpritz-DisProt", {})
        esp_labels_all = esp.get("label", "").split(",")
        esp_scores_all = esp.get("residue", "").split(",")
        full_seq = phasepred.get(uid, {}).get("Sequence", "")

        if not full_seq or len(esp_labels_all) != len(full_seq) or len(esp_scores_all) != len(full_seq):
            n_skip += 1
            continue

        # Build topology mask to select correct residue subset
        topo_labels = build_region_labels(uid, len(full_seq))
        mask = [lab == region for lab in topo_labels]

        # Verify: masked sequence must match concat FASTA sequence
        masked_seq = "".join(c for c, m in zip(full_seq, mask) if m)
        if masked_seq != concat_seq:
            # Topology cache / sequence mismatch — use len(concat_seq) as fallback
            # by matching by position count rather than bailing out
            print(f"    [WARN] {uid} masked_seq != concat_seq  "
                  f"(masked={len(masked_seq)}, fasta={len(concat_seq)}) — skipping")
            n_skip += 1
            continue

        sub_labels = [l for l, m in zip(esp_labels_all, mask) if m]
        sub_scores = [s for s, m in zip(esp_scores_all, mask) if m]

        region_slug = region.replace("/", "_")
        stem = f"{uid}_{region_slug}_concat"

        fasta_content   = f">{uid}\n{concat_seq}\n"
        espritz_content = make_espritz_content(sub_labels, sub_scores)

        score = run_fuzdrop(stem, fasta_content, espritz_content, work_dir)

        if score is not None:
            concat_rows.append({
                "UniProt_ID":    uid,
                "region":        region,
                "length":        len(concat_seq),
                "fuzdrop_score": score,
            })
            n_ok += 1
        else:
            print(f"    [FAIL] {uid} {region} len={len(concat_seq)}")
            n_fail += 1

        total_count += 1

shutil.rmtree(work_dir, ignore_errors=True)

concat_df = pd.DataFrame(concat_rows)
concat_out = RAW_OUT / "FuzDrop_concat_raw.csv"
concat_df.to_csv(concat_out, index=False)
print(f"\n[Concat] Done: {n_ok} scored, {n_fail} failed, {n_skip} skipped (no ESpritz data)")
print(f"  Saved {len(concat_df)} rows to {concat_out.relative_to(ROOT)}")

print("\n=== SUMMARY ===")
print(f"  Segments: {len(seg_df)} rows  -> {seg_out}")
print(f"  Concat:   {len(concat_df)} rows  -> {concat_out}")
