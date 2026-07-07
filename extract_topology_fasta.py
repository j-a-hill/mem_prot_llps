"""
extract_topology_fasta.py

Writes per-region sub-sequence FASTA files for the 60 membrane LLPS proteins
(Cytoplasmic / Transmembrane / Extracellular-Lumenal / Whole), using the same
UniProt topology cache as extract_topology_scores.py. Two purposes:

  1. Automates the 2 remaining tools with locally-runnable code that wasn't
     already covered by extract_topology_scores.py's per-residue extraction:
       - R+Y: trivial Arg+Tyr composition, computed directly
       - LLPhyScore: standalone package run on each region FASTA
     Both get merged into output/topology_domain_scores.csv.

  2. Leaves the FASTA files in output/topology_fasta/ for manually submitting
     to the predictors that have no local runnable code or per-residue data
     (PICNIC, PICNIC (GO), PSPire, PSPHunter, SaPS, PdPS, PDL, FuzDrop, PSAP,
     DeepPhase) — see output/topology_manual_retrieval_checklist.csv.

Outputs
-------
  output/topology_fasta/{Cytoplasmic,Transmembrane,Extracellular_Lumenal,Whole}.fasta
  output/topology_domain_scores.csv          (RY_*, LLPhyScore_* columns added)
  output/topology_manual_retrieval_checklist.csv
"""

import json
import re
import subprocess
import tempfile
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT        = Path(__file__).parent
PRED        = ROOT / "Predictors_whole_genome_sets"
OUT         = ROOT / "output"
FIG_DIR     = OUT / "figures"
FASTA_DIR   = OUT / "topology_fasta"
TOPO_CACHE  = ROOT / "data" / "uniprot_topology_cache.csv"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
LLPHYSCORE_DIR = ROOT / "external_tools" / "LLPhyScore" / "standalone_package" / "LLPhyScore"
FUZDROP_DIR = ROOT / "external_tools" / "FuzDrop"
FUZDROP_BIN = FUZDROP_DIR / "FuzDrop"
PHASEPRED_JSON = PRED / "PhaSePred_human_reviewed.json"

FASTA_DIR.mkdir(parents=True, exist_ok=True)

REGIONS = ["Cytoplasmic", "Transmembrane", "Extracellular/Lumenal"]
REGION_FILE_SLUG = {
    "Cytoplasmic": "Cytoplasmic", "Transmembrane": "Transmembrane",
    "Extracellular/Lumenal": "Extracellular_Lumenal", "Whole": "Whole",
}

# ── Load the 60 positives + sequences (same as extract_topology_scores.py) ──

pos_ids = pd.read_csv(ROOT / "exp_db" / "membrane_exp_db_matches.csv")["Entry"].tolist()

seqs = {}
uid, chunks = None, []
with open(PRED / "Seq2Phase_swiss_prot_human_220916.fasta") as f:
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

# ── Topology cache (must already exist from extract_topology_scores.py) ────

if not TOPO_CACHE.exists():
    raise SystemExit(f"{TOPO_CACHE} not found — run extract_topology_scores.py first.")
topo_raw = pd.read_csv(TOPO_CACHE).rename(columns={
    "Entry": "UniProt_ID", "Topological domain": "topo_dom_raw", "Transmembrane": "transmem_raw",
}).set_index("UniProt_ID")

FEATURE_RE = re.compile(r'(TOPO_DOM|TRANSMEM)\s+(\d+)\.\.(\d+);\s*/note="([^"]*)"')


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


# ── Extract contiguous sub-sequences per region (concatenated if >1 span) ──

region_seqs = {r: {} for r in REGIONS}
region_seqs["Whole"] = {}

for uid in pos_ids:
    seq = seqs.get(uid)
    if not seq:
        continue
    region_seqs["Whole"][uid] = seq
    labels = build_region_labels(uid, len(seq))
    for region in REGIONS:
        sub = "".join(c for c, lab in zip(seq, labels) if lab == region)
        if sub:
            region_seqs[region][uid] = sub

# ── Write FASTA files ────────────────────────────────────────────────────────

for region, slug in REGION_FILE_SLUG.items():
    fpath = FASTA_DIR / f"{slug}.fasta"
    with open(fpath, "w") as f:
        for uid, sub in region_seqs.get(region, {}).items():
            f.write(f">{uid}\n{sub}\n")
    n = len(region_seqs.get(region, {}))
    print(f"Wrote {fpath.relative_to(ROOT)}  ({n} sequences)")

# ── R+Y: trivial, computed directly per region ──────────────────────────────

ry_rows = []
for uid in pos_ids:
    row = {"UniProt_ID": uid}
    for region, slug in REGION_FILE_SLUG.items():
        sub = region_seqs.get(region, {}).get(uid)
        row[f"RY_{region}" if region != "Whole" else "RY_whole"] = (
            (sub.count("R") + sub.count("Y")) / len(sub) if sub else np.nan
        )
    ry_rows.append(row)
ry_df = pd.DataFrame(ry_rows)
print(f"\nComputed R+Y for {ry_df['RY_whole'].notna().sum()}/{len(ry_df)} proteins")

# ── LLPhyScore: run the standalone package on each region FASTA ────────────

llphyscore_cols = {}
if LLPHYSCORE_DIR.exists():
    for region, slug in REGION_FILE_SLUG.items():
        fasta_path = FASTA_DIR / f"{slug}.fasta"
        if fasta_path.stat().st_size == 0:
            continue
        out_csv = Path(tempfile.mktemp(suffix=".csv"))
        cmd = ["python3", "LLPhyScore_standalone.py", "-i", str(fasta_path.resolve()),
               "-s", "raw", "-o", str(out_csv)]
        print(f"\nRunning LLPhyScore on {slug} ({len(region_seqs.get(region, {}))} sequences)...")
        result = subprocess.run(cmd, cwd=LLPHYSCORE_DIR, capture_output=True, text=True)
        if result.returncode != 0 or not out_csv.exists():
            print(f"  [skipped] LLPhyScore failed for {slug}: {result.stderr[-500:]}")
            continue
        scores = pd.read_csv(out_csv)
        scores["UniProt_ID"] = scores["tag"].str.split("|").str[1].fillna(scores["tag"])
        # "8-feature sum" is a SUM over residues (scales with sequence length) — divide by
        # length so a short region's score is comparable to the whole protein's, same fix
        # as catGRANULE/PScore/ParSe2 above.
        lengths = scores["UniProt_ID"].map(lambda u: len(region_seqs.get(region, {}).get(u, "")))
        colname = f"LLPhyScore_{region}" if region != "Whole" else "LLPhyScore_whole"
        llphyscore_cols[colname] = (scores.set_index("UniProt_ID")["8-feature sum"] /
                                     lengths.set_axis(scores["UniProt_ID"]).replace(0, np.nan))
        out_csv.unlink(missing_ok=True)
        print(f"  scored {len(scores)} sequences")
else:
    print(f"\n[skipped] LLPhyScore standalone package not found at {LLPHYSCORE_DIR} "
          f"— re-clone julie-forman-kay-lab/LLPhyScore if you want this filled in.")

llphyscore_df = pd.DataFrame(llphyscore_cols).reset_index().rename(columns={"index": "UniProt_ID"}) \
    if llphyscore_cols else pd.DataFrame({"UniProt_ID": pos_ids})

# ── Load PhaSePred JSON once for SaPS/PdPS whole scores + FuzDrop's ESpritz input ──

phasepred = {}
if PHASEPRED_JSON.exists():
    print("\nLoading PhaSePred JSON for SaPS/PdPS whole scores + FuzDrop's ESpritz-DisProt input...")
    with open(PHASEPRED_JSON) as f:
        phasepred = json.load(f)

# ── SaPS / PdPS: whole-protein score is already in the PhaSePred JSON ──────
# PhaSePred's gradient-boosted models combine several whole-protein-level
# features (e.g. phosphosite count) that can't be recomputed from a residue
# mask, so the Cytoplasmic/Transmembrane/Extracellular splits still need
# manual submission to predict.phasep.pro — only "whole" is free here.

saps_pdps_df = pd.DataFrame([
    {
        "UniProt_ID": uid,
        "SaPS_whole": phasepred.get(uid, {}).get("PhaSePred", {}).get("SaPS-10fea", np.nan),
        "PdPS_whole": phasepred.get(uid, {}).get("PhaSePred", {}).get("PdPS-10fea", np.nan),
    }
    for uid in pos_ids
])
print(f"  SaPS/PdPS: filled whole-protein score for "
      f"{saps_pdps_df['SaPS_whole'].notna().sum()}/{len(pos_ids)} proteins from local JSON")

# ── FuzDrop: local Linux binary (Fuxreiter lab), region-sliced ESpritz input ──
# FuzDrop expects a paired .espritz file (the "NMR" Espritz flavour from a
# separate Perl program). We reuse the ESpritz-DisProt residue arrays already
# cached in the PhaSePred JSON instead of installing a second disorder
# predictor — same data source extract_topology_scores.py uses for the other
# residue-level tools, just repackaged into FuzDrop's expected format. The
# header below is copied verbatim from FuzDrop's own example data: the binary
# hard-codes skipping that many lines before reading scores, so it must be
# present even though its content (a license notice) is otherwise inert.
ESPRITZ_HEADER = (
    "*" * 108 + "\n\n"
    "Licensed to: Prof. MONIKA FUXREITER (University of Padova) located in Padova, Italy. \n"
    "This license is for non-commercial use only. Please see LICENSE file for details \n"
    "(https://biocomputingup.it/assets/data/LICENSE.txt)\n\n"
    "Contact silvio.tosatto@unipd.it for commercial licensing details.\n\n\n"
    + "*" * 108 + "\n"
)

fuzdrop_rows = {region: {} for region in REGION_FILE_SLUG}
if FUZDROP_BIN.exists() and phasepred:
    work_dir = Path(tempfile.mkdtemp(prefix="fuzdrop_"))
    n_run = 0
    for uid in pos_ids:
        seq = seqs.get(uid)
        esp = phasepred.get(uid, {}).get("ESpritz-DisProt", {})
        esp_labels = esp.get("label", "").split(",")
        esp_scores = esp.get("residue", "").split(",")
        if not seq or len(esp_labels) != len(seq) or len(esp_scores) != len(seq):
            continue
        topo_labels = build_region_labels(uid, len(seq))

        for region in REGION_FILE_SLUG:
            mask = [True] * len(seq) if region == "Whole" else [lab == region for lab in topo_labels]
            sub_seq = "".join(c for c, m in zip(seq, mask) if m)
            if not sub_seq:
                continue
            sub_labels = (l for l, m in zip(esp_labels, mask) if m)
            sub_scores = (s for s, m in zip(esp_scores, mask) if m)

            fasta_path = work_dir / f"{uid}.fasta"
            espritz_path = work_dir / f"{uid}.espritz"
            res_path = work_dir / f"{uid}_res.txt"
            fasta_path.write_text(f">{uid}\n{sub_seq}\n")
            espritz_path.write_text(
                ESPRITZ_HEADER + "\n".join(f"{l}\t{s}" for l, s in zip(sub_labels, sub_scores)) + "\n"
            )
            result = subprocess.run([str(FUZDROP_BIN.resolve()), fasta_path.name, espritz_path.name],
                                     cwd=work_dir, capture_output=True, text=True)
            if result.returncode == 0 and res_path.exists():
                last_line = res_path.read_text().strip().splitlines()[-1]
                try:
                    fuzdrop_rows[region][uid] = float(last_line.split("=")[-1].strip())
                    n_run += 1
                except ValueError:
                    pass
            for p in (fasta_path, espritz_path, res_path):
                p.unlink(missing_ok=True)
    print(f"  FuzDrop: scored {n_run} protein x region combinations")
else:
    print(f"\n[skipped] FuzDrop binary not found at {FUZDROP_BIN} or PhaSePred JSON missing "
          f"— download FuzDrop_linux.zip from fuxreiterlab.github.io if you want this filled in.")

fuzdrop_cols_named = {
    (f"FuzDrop_{region}" if region != "Whole" else "FuzDrop_whole"): pd.Series(scores)
    for region, scores in fuzdrop_rows.items() if scores
}
fuzdrop_df = pd.DataFrame(fuzdrop_cols_named).reset_index().rename(columns={"index": "UniProt_ID"}) \
    if fuzdrop_cols_named else pd.DataFrame({"UniProt_ID": pos_ids})

# ── Merge into topology_domain_scores.csv ───────────────────────────────────

topo_scores_path = OUT / "topology_domain_scores.csv"
topo_df = pd.read_csv(topo_scores_path)
for df_new in (ry_df, llphyscore_df, saps_pdps_df, fuzdrop_df):
    new_cols = [c for c in df_new.columns if c != "UniProt_ID"]
    topo_df = topo_df.drop(columns=[c for c in new_cols if c in topo_df.columns], errors="ignore")
    topo_df = topo_df.merge(df_new, on="UniProt_ID", how="left")
topo_df.to_csv(topo_scores_path, index=False)
print(f"\nUpdated {topo_scores_path.relative_to(ROOT)} with RY_*, LLPhyScore_*, SaPS/PdPS_whole, "
      f"and FuzDrop_* columns")

# ── Regenerate the paired delta summary + figure with all tools now present ──

import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

ALL_TOOLS = sorted({c[:-6] for c in topo_df.columns if c.endswith("_whole")})

summary_rows = []
for tool in ALL_TOOLS:
    for region in REGIONS:
        if f"{tool}_{region}" not in topo_df.columns:
            # whole-only tool (e.g. SaPS/PdPS — no region-split data available)
            summary_rows.append({"Tool": tool, "Region": region, "n": 0,
                                  "mean_delta": np.nan, "median_delta": np.nan,
                                  "frac_region_higher": np.nan, "wilcoxon_p": np.nan})
            continue
        sub = topo_df[[f"{tool}_whole", f"{tool}_{region}"]].dropna()
        n = len(sub)
        if n < 3:
            summary_rows.append({"Tool": tool, "Region": region, "n": n,
                                  "mean_delta": np.nan, "median_delta": np.nan,
                                  "frac_region_higher": np.nan, "wilcoxon_p": np.nan})
            continue
        delta = sub[f"{tool}_{region}"] - sub[f"{tool}_whole"]
        try:
            _, p = wilcoxon(delta)
        except ValueError:
            p = np.nan
        summary_rows.append({
            "Tool": tool, "Region": region, "n": n,
            "mean_delta": round(float(delta.mean()), 4),
            "median_delta": round(float(delta.median()), 4),
            "frac_region_higher": round(float((delta > 0).mean()), 3),
            "wilcoxon_p": round(float(p), 4) if not np.isnan(p) else np.nan,
        })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(OUT / "topology_score_deltas.csv", index=False)
print("\nUpdated output/topology_score_deltas.csv")
print(summary_df.to_string(index=False))

REGION_COLOR = {"Cytoplasmic": "#0072B2", "Transmembrane": "#888888", "Extracellular/Lumenal": "#e74c3c"}

fig, ax = plt.subplots(figsize=(14, 5.5))
x = np.arange(len(ALL_TOOLS))
bar_w = 0.25
for i, region in enumerate(REGIONS):
    sub = summary_df[summary_df["Region"] == region].set_index("Tool").reindex(ALL_TOOLS)
    offset = (i - 1) * bar_w
    ax.bar(x + offset, sub["mean_delta"], bar_w, color=REGION_COLOR[region], alpha=0.85, label=region)
    for j, r in enumerate(sub.itertuples()):
        if pd.notna(r.wilcoxon_p) and r.wilcoxon_p < 0.05:
            ax.text(x[j] + offset, r.mean_delta + (0.01 if r.mean_delta >= 0 else -0.01),
                    "*", ha="center", va="bottom" if r.mean_delta >= 0 else "top",
                    fontsize=11, color=REGION_COLOR[region])
ax.axhline(0, color="#333", linewidth=0.8)
ax.set_xticks(x)
ax.set_xticklabels(ALL_TOOLS, rotation=20, ha="right")
ax.set_ylabel("Mean Δ (region mean score − whole-protein score)")
ax.set_title("Does the region score higher than the whole protein?  (* = Wilcoxon p<0.05)")
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(FIG_DIR / "topology_score_deltas.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Updated output/figures/topology_score_deltas.png")

# ── Manual-retrieval checklist for the remaining tools ──────────────────────
# PICNIC, PICNIC (GO), PSPire, PSPHunter, PSAP, and DeepPhase were already
# completed via manual web-submission batches (see bin_fasta_for_submission.py)
# and are merged into topology_domain_scores.csv directly, so they're not
# tracked here anymore. PDL (PSPsPredict) has local code at
# external_tools/PSPsPredict but is parked pending a safe way to handle very
# long sequences on this machine — see training_set_provenance.md.

MANUAL_TOOLS = {
    "SaPS": "PhaSePred web server: predict.phasep.pro (accepts FASTA upload). Whole-protein "
            "score already filled from local PhaSePred_human_reviewed.json — only the 3 "
            "region FASTA below need manual submission.",
    "PdPS": "PhaSePred web server: predict.phasep.pro (accepts FASTA upload, same submission "
            "as SaPS). Whole-protein score already filled from local PhaSePred_human_reviewed"
            ".json — only the 3 region FASTA below need manual submission.",
}

checklist_rows = []
for tool, note in MANUAL_TOOLS.items():
    for region in REGIONS:
        checklist_rows.append({
            "Predictor": tool,
            "Region": region,
            "FASTA file": f"output/topology_fasta/{REGION_FILE_SLUG[region]}.fasta",
            "Submission method": note,
            "Score (paste here)": "",
        })

checklist_df = pd.DataFrame(checklist_rows)
checklist_path = OUT / "topology_manual_retrieval_checklist.csv"
checklist_df.to_csv(checklist_path, index=False)
print(f"Saved {checklist_path.relative_to(ROOT)}  "
      f"({len(MANUAL_TOOLS)} predictors x {len(REGIONS)} regions)")
