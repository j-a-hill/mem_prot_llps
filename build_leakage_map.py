"""
build_leakage_map.py

Phase 0c: Cross-reference benchmark proteins against predictor training sets.

Tracks BOTH positive and negative leakage separately:
  - positive leakage: benchmark protein was a training POSITIVE
    → score likely inflated
  - negative leakage: benchmark protein was a training NEGATIVE
    → score likely deflated (model was trained to predict it as non-LLPS)

Reads training set files from  training_sets/  and produces:

  output/leakage_map.csv     — one row per membrane protein; pos/neg flag per tool
  output/leakage_summary.csv — per-tool counts of pos/neg/clean/contradictory proteins
  output/globally_clean.csv  — proteins not in ANY tool's training set (pos or neg)
"""

import pickle
from pathlib import Path

import pandas as pd

ROOT      = Path(__file__).parent
TRAIN_DIR = ROOT / "training_sets"
PRED      = ROOT / "Predictors_whole_genome_sets"
OUT       = ROOT / "output"

MIN_FRAGMENT_LEN = 15  # avoid spurious short-substring matches


def _parse_fasta_dict(path, key_fn=lambda header: header):
    """Parse a FASTA file into {key: sequence}; key_fn maps the raw header to a key."""
    seqs = {}
    key, chunks = None, []
    for line in Path(path).read_text().splitlines():
        line = line.rstrip()
        if line.startswith(">"):
            if key is not None:
                seqs[key] = "".join(chunks)
            key = key_fn(line[1:].strip())
            chunks = []
        else:
            chunks.append(line)
    if key is not None:
        seqs[key] = "".join(chunks)
    return seqs


def _match_by_sequence(query_seqs, target_seqs, min_fragment_len=MIN_FRAGMENT_LEN):
    """
    UniProt IDs (target_seqs keys) whose sequence exactly matches, contains,
    or is contained within any query sequence. Used when a training set is
    keyed by construct/gene name rather than UniProt ID — recovers real
    human proteins/domains/fusions among non-mappable names by sequence
    identity instead.
    """
    matched = set()
    for seq in query_seqs:
        hits = [uid for uid, tseq in target_seqs.items() if tseq == seq]
        if not hits and len(seq) >= min_fragment_len:
            hits = [uid for uid, tseq in target_seqs.items() if seq in tseq]
        if not hits:
            hits = [uid for uid, tseq in target_seqs.items()
                     if len(tseq) >= min_fragment_len and tseq in seq]
        matched.update(hits)
    return matched


# ── File format handlers ──────────────────────────────────────────────────────

def load_ids(cfg):
    """
    Return a set of UniProt IDs from a training-set config dict, or None if
    the file is missing (unknown leakage status).
    """
    if cfg is None:
        return set()

    fmt = cfg.get("fmt", "txt_lines")

    # Multiple FASTA files merged (e.g. MambaPhase scaffold + client)
    if fmt == "fasta_sp_multi":
        ids = set()
        for fname in cfg["files"]:
            fp = TRAIN_DIR / fname
            if not fp.exists():
                return None
            for line in fp.read_text().splitlines():
                if line.startswith(">"):
                    parts = line.split("|")
                    if len(parts) >= 2:
                        ids.add(parts[1].strip())
        return ids

    fpath = TRAIN_DIR / cfg["file"]
    if not fpath.exists():
        return None

    # File exists but its IDs can't be mapped to UniProt accessions (e.g. gene-name
    # headers) — leakage status is genuinely unverifiable, not "confirmed clean".
    if fmt == "unmappable":
        return None

    # Pickled list of "ACCESSION_GENE_ORGANISM"-style tags (LLPhyScore's processed
    # training tags) — UniProt accession is the substring before the first "_".
    if fmt == "pkl_tag_prefix":
        tags = pickle.load(open(fpath, "rb"))
        return {t.split("_")[0] for t in tags}

    # Construct/gene-name-keyed FASTA (e.g. LLPhyScore's true-positive set) —
    # resolve by matching sequences directly against the human proteome.
    # cfg["tags_file"]: pickled list restricting to a specific split (e.g.
    # training-only, excluding the held-out test split).
    if fmt == "seq_match_proteome":
        restrict_tags = None
        if cfg.get("tags_file"):
            tags_path = TRAIN_DIR / cfg["tags_file"]
            if not tags_path.exists():
                return None
            restrict_tags = set(pickle.load(open(tags_path, "rb")))
        seqs = _parse_fasta_dict(fpath)
        if restrict_tags is not None:
            seqs = {k: v for k, v in seqs.items() if k in restrict_tags}
        proteome_path = PRED / "Seq2Phase_swiss_prot_human_220916.fasta"
        if not proteome_path.exists():
            return None
        proteome = _parse_fasta_dict(proteome_path, key_fn=lambda h: h.split("|")[1] if "|" in h else h)
        return _match_by_sequence(seqs.values(), proteome)

    if fmt == "txt_lines":
        return {l.strip() for l in fpath.read_text().splitlines() if l.strip()}

    if fmt == "txt_space":
        return {i.strip() for i in fpath.read_text().split() if i.strip()}

    if fmt == "fasta_bare":
        # >ACCESSION (bare UniProt accession, possible trailing whitespace)
        return {l[1:].strip() for l in fpath.read_text().splitlines() if l.startswith(">")}

    if fmt == "fasta_sp":
        # >sp|ACCESSION|NAME_ORGANISM ...
        ids = set()
        for line in fpath.read_text().splitlines():
            if line.startswith(">"):
                parts = line.split("|")
                if len(parts) >= 2:
                    ids.add(parts[1].strip())
        return ids

    if fmt == "csv":
        df = pd.read_csv(fpath)
        col = cfg["col"]
        return _clean(df[col])

    if fmt == "xlsx":
        sheet  = cfg.get("sheet", 0)
        header = cfg.get("header", 0)
        col    = cfg["col"]
        df = pd.read_excel(fpath, sheet_name=sheet, header=header)
        return _clean(df[col])

    if fmt == "xlsx_filter":
        sheet  = cfg.get("sheet", 0)
        header = cfg.get("header", 0)
        col    = cfg["col"]
        df = pd.read_excel(fpath, sheet_name=sheet, header=header)
        for fcol, fvals in cfg.get("filter", {}).items():
            vals = fvals if isinstance(fvals, list) else [fvals]
            df = df[df[fcol].isin(vals)]
        return _clean(df[col])

    if fmt == "xlsx_cols":
        # Multiple columns contribute to the same set (e.g. Seq2Phase scaffold + client)
        sheet  = cfg.get("sheet", 0)
        header = cfg.get("header", 0)
        cols   = cfg["cols"]
        df = pd.read_excel(fpath, sheet_name=sheet, header=header)
        ids = set()
        for c in cols:
            if c in df.columns:
                ids.update(_clean(df[c]))
        return ids

    raise ValueError(f"Unknown format: {fmt!r}")


def _clean(series):
    return {str(v).strip() for v in series.dropna()
            if str(v).strip() not in ("", "nan")}


# ── Training set registry ─────────────────────────────────────────────────────
# Each entry: tool → {pos: cfg, neg: cfg, note: str}
# cfg = None means no explicit set (treated as empty / no leakage for that direction)

TRAINING_SETS = {
    "PSAP": {
        "pos": {"file": "PSAP_S1_training_set.xlsx", "fmt": "xlsx",
                "sheet": "Human", "col": "Uniprot_ID"},
        "extra_pos": {"O43561"},
        "neg": None,
        "note": "EMBO Mol Med 2021; 90 human proteins (paper's Table S1) + O43561 (LAT_HUMAN), "
                "present in the actual shipped model's bundled training-ID list "
                "(psap package's data/assets/uniprot_ids.txt, 91 IDs) but absent from the "
                "paper's table; no explicit negatives (proteome BG)",
    },
    "DeePhase": {
        "pos": {"file": "DeePhase_llps_plus.csv",  "fmt": "csv", "col": "Uniprot_ID"},
        "neg": {"file": "DeePhase_llps_minus.csv", "fmt": "csv", "col": "Uniprot_ID"},
        "note": "github.com/kadiliissaar/DeePhase; in vitro homotypic <100µM filter",
    },
    "PSPHunter": {
        "pos": {"file": "PSPHunter_training_pos_hPS167.txt",      "fmt": "txt_lines"},
        "neg": {"file": "PSPHunter_training_neg_non-hPS5754.txt", "fmt": "txt_space"},
        "note": "Nat Commun 15:2662; hPS167 positives + 5754 negatives",
    },
    "PSPire": {
        "pos": {"file": "PSPire_S4_training_sets.xlsx", "fmt": "xlsx_filter",
                "sheet": "Datasets", "col": "UniprotEntry",
                "filter": {"Datasets": "Training", "Type": ["ID-PSP", "noID-PSP"]}},
        "neg": {"file": "PSPire_S4_training_sets.xlsx", "fmt": "xlsx_filter",
                "sheet": "Datasets", "col": "UniprotEntry",
                "filter": {"Datasets": "Training", "Type": ["non-PSP"]}},
        "note": "Nat Commun 15:2147 Suppl Data 4; 259 pos + 8323 neg (Training split only)",
    },
    "PICNIC": {
        "pos": {"file": "PICNIC_S1_training_sets.xlsx", "fmt": "xlsx",
                "sheet": "All(train+test)_positive", "col": "uniprot_id"},
        "neg": {"file": "PICNIC_S1_training_sets.xlsx", "fmt": "xlsx",
                "sheet": "All(train+test)_negative", "col": "uniprot_id"},
        "note": "Edmond doi:10.17617/3.0Y9Q8N; CD-CODE-derived; 2142 pos + 1709 neg",
    },
    "PDL": {
        "pos": {"file": "PSPsPredict_S1_training_sets.xlsx", "fmt": "xlsx_filter",
                "sheet": "TrainData", "col": "name", "filter": {"label": [1]}},
        "neg": {"file": "PSPsPredict_S1_training_sets.xlsx", "fmt": "xlsx_filter",
                "sheet": "TrainData", "col": "name", "filter": {"label": [0]}},
        "note": "github.com/xmuzhanglab/PSPsPredict; ProtT5+KmerConv; 448 pos + 525 neg",
    },
    "SaPS": {
        "pos": {"file": "Phasepro_S2_training_sets.xlsx", "fmt": "xlsx",
                "sheet": "hSaPS", "col": "UniprotEntry", "header": 1},
        "neg": {"file": "Phasepro_S2_training_sets.xlsx", "fmt": "xlsx",
                "sheet": "hNoPS", "col": "UniprotEntry", "header": 1},
        "note": "PNAS 119:e2115369119 Suppl S2; hSaPS pos (59) + hNoPS neg (8801)",
    },
    "PdPS": {
        "pos": {"file": "Phasepro_S2_training_sets.xlsx", "fmt": "xlsx",
                "sheet": "hPdPS", "col": "UniprotEntry", "header": 1},
        "neg": {"file": "Phasepro_S2_training_sets.xlsx", "fmt": "xlsx",
                "sheet": "hNoPS", "col": "UniprotEntry", "header": 1},
        "note": "PNAS 119:e2115369119 Suppl S2; hPdPS pos (96); same negatives as SaPS",
    },
    "MolPhase": {
        "pos": {"file": "Molphase_D_EV1_training_set.xlsx", "fmt": "xlsx",
                "sheet": "Positive training set", "col": "ID"},
        "neg": None,
        "note": "EMBO J 2024; 606 positives from multi-DB; PDB-derived negatives have no UniProt IDs",
    },
    "Opt_PredLLPS": {
        "pos": {"file": "Opt_PredLLPS_LLPS_training_set.fasta",     "fmt": "fasta_bare"},
        "neg": {"file": "Opt_PredLLPS_non_LLPS_training_set.fasta", "fmt": "fasta_bare"},
        "note": "Briefings Bioinform 2024; successor to PredLLPS_PSSM; 798 pos + 798 neg",
    },
    "MambaPhase": {
        "pos": {"fmt": "fasta_sp_multi",
                "files": ["mambaphase_drllps_scaffold_clstr_Homo_sapiens.fasta",
                          "mambaphase_drllps_client_clstr_Homo_sapiens.fasta"]},
        "neg": {"file": "mambaphase_drllps_nonllps_clstr_Homo_sapiens.fasta", "fmt": "fasta_sp"},
        "note": "Briefings Bioinform 2025; ESM2+Mamba; human subset; scaffold+client=pos (2775); nonllps=neg (12562)",
    },
    "Seq2Phase": {
        "pos": {"file": "seq2phase_Supplementary_Data_training_sets.xlsx", "fmt": "xlsx_cols",
                "sheet": "S1", "header": 3, "cols": ["scaffold", "client"]},
        "neg": {"file": "seq2phase_Supplementary_Data_training_sets.xlsx", "fmt": "xlsx_cols",
                "sheet": "S1", "header": 3, "cols": ["non-LLPS"]},
        "note": "Bioinformatics Advances 2023; DrLLPS clients; human scaffold+client=pos (2775); non-LLPS=neg (12562)",
    },
    "LLPhyScore": {
        "pos": {"file": "LLPhyScore_positive_untagged_20191126_training_set.fasta",
                "fmt": "unmappable"},
        # Most of the 305 training-positive names are unmappable (constructs, fusions,
        # non-human orthologs, synthetic repeats), but ~22% resolve to real human
        # proteins/domains by exact sequence match — used as a *positive-only*
        # override below (non-matches stay "unknown", not falsely "clean").
        "pos_override": {"file": "LLPhyScore_positive_untagged_20191126_training_set.fasta",
                          "tags_file": "LLPhyScore_positive_training_tags.pkl",
                          "fmt": "seq_match_proteome"},
        "neg": {"file": "LLPhyScore_human_negatives_training_tags.pkl",
                "fmt": "pkl_tag_prefix"},
        "note": "github.com/julie-forman-kay-lab/LLPhyScore; 305 train pos (construct/gene "
                "names, not UniProt IDs — mostly unverifiable, ~22% sequence-resolved); "
                "2000 train neg (generic human-proteome sample, UniProt-accession-keyed — "
                "negative leakage IS checkable)",
    },
    # Theory-based / rule-based — no LLPS training sets:
    # FuzDrop, ParSe 2.0, catGRANULE, PLAAC, PScore, ESpritz, SEG, R+Y
}

# ── Load benchmark membrane proteins ─────────────────────────────────────────

master = pd.read_csv(OUT / "master_table.csv")
mem_proteins = master[master["is_membrane"]]["UniProt_ID"].tolist()
print(f"Benchmark: {len(mem_proteins)} membrane proteins")

# ── Load all training sets ────────────────────────────────────────────────────

pos_sets = {}   # tool → set of UniProt IDs (or None = file missing)
neg_sets = {}

print()
for tool, cfg in TRAINING_SETS.items():
    pos = load_ids(cfg["pos"])
    neg = load_ids(cfg["neg"])
    if pos is not None and cfg.get("extra_pos"):
        pos = pos | set(cfg["extra_pos"])
    pos_sets[tool] = pos
    neg_sets[tool] = neg
    pos_n   = f"{len(pos)}"  if pos  is not None else "MISSING"
    neg_n   = f"{len(neg)}"  if neg  is not None else ("—" if cfg["neg"] is None else "MISSING")
    status  = "OK" if (pos is not None) else "MISSING"
    print(f"  [{status}] {tool:<14} pos={pos_n:>6}  neg={neg_n:>6}  | {cfg['note'][:60]}")

# Positive-only overrides: confirmed positives recovered by sequence matching
# even when the primary "pos" set is unmappable. Non-matches are left exactly
# as the primary set already determined (typically "unknown") rather than
# being upgraded to "confirmed clean".
pos_overrides = {}
for tool, cfg in TRAINING_SETS.items():
    if cfg.get("pos_override"):
        override = load_ids(cfg["pos_override"])
        pos_overrides[tool] = override or set()
        print(f"  [override] {tool:<14} sequence-matched pos={len(pos_overrides[tool])}")

# ── Build leakage map ─────────────────────────────────────────────────────────

rows = []
for uid in mem_proteins:
    row = {"UniProt_ID": uid}
    for tool in TRAINING_SETS:
        p = pos_sets[tool]
        n = neg_sets[tool]
        row[f"pos_{tool}"] = (uid in p)       if p is not None else "unknown"
        row[f"neg_{tool}"] = (uid in n)       if n is not None else (
                              False            if TRAINING_SETS[tool]["neg"] is None
                              else "unknown")
        if tool in pos_overrides and uid in pos_overrides[tool]:
            row[f"pos_{tool}"] = True
    rows.append(row)

leakage_df = pd.DataFrame(rows)
leakage_df.to_csv(OUT / "leakage_map.csv", index=False)
print(f"\nSaved output/leakage_map.csv  ({len(leakage_df)} proteins × {len(leakage_df.columns)-1} columns)")

# ── Summary table ─────────────────────────────────────────────────────────────

summary = []
for tool in TRAINING_SETS:
    pc = f"pos_{tool}"
    nc = f"neg_{tool}"

    known_pos = leakage_df[leakage_df[pc] != "unknown"]
    known_neg = leakage_df[leakage_df[nc] != "unknown"]

    n_pos_leaked = (known_pos[pc] == True).sum()
    n_neg_leaked = (known_neg[nc] == True).sum()

    if pc in leakage_df.columns and nc in leakage_df.columns:
        both_known = leakage_df[(leakage_df[pc] != "unknown") & (leakage_df[nc] != "unknown")]
        n_contra = ((both_known[pc] == True) & (both_known[nc] == True)).sum()
    else:
        n_contra = 0

    # Clean = not pos-leaked AND not neg-leaked (both must be confirmed False)
    clean_mask = (leakage_df[pc] == False) & (leakage_df[nc] == False)
    clean_ids  = leakage_df[clean_mask]["UniProt_ID"].tolist()

    summary.append({
        "Tool":            tool,
        "N_pos_leaked":    int(n_pos_leaked),
        "N_neg_leaked":    int(n_neg_leaked),
        "N_contradictory": int(n_contra),
        "N_clean":         len(clean_ids),
        "Clean_IDs":       ";".join(clean_ids),
        "Note":            TRAINING_SETS[tool]["note"][:80],
    })

summary_df = pd.DataFrame(summary)
summary_df.to_csv(OUT / "leakage_summary.csv", index=False)
print("Saved output/leakage_summary.csv")

# ── Globally clean set ────────────────────────────────────────────────────────

def is_globally_clean(row):
    """
    Clean across every tool whose leakage status is actually determinable.
    "unknown" (e.g. LLPhyScore's gene-name-keyed training set) means no
    information either way, so it must not disqualify a protein — otherwise
    a single unverifiable tool collapses the global-clean set to zero for
    every predictor, not just that tool.
    """
    for tool in TRAINING_SETS:
        if row.get(f"pos_{tool}") is True:
            return False
        if row.get(f"neg_{tool}") is True:
            return False
    return True

globally_clean = leakage_df[leakage_df.apply(is_globally_clean, axis=1)]["UniProt_ID"].tolist()
pd.DataFrame({"UniProt_ID": globally_clean}).to_csv(OUT / "globally_clean.csv", index=False)
print(f"Saved output/globally_clean.csv  (N={len(globally_clean)} proteins clean across all tools)")

# ── Clean masks CSV (one row per protein, one col per tool) ──────────────────
# Maps score column names → tool name for use in analysis scripts.
# Theory-based tools (FuzDrop, catGRANULE, PLAAC, PScore, ESpritz, SEG) map to None
# meaning all 60 proteins are clean for those predictors.

SCORE_TO_TOOL = {
    "PICNIC_score":    "PICNIC",
    "PICNIC_GO_score": "PICNIC",
    "PSAP_score":      "PSAP",
    "PSPHunter_prob":  "PSPHunter",
    "PSPire_score":    "PSPire",
    "FuzDrop_pLLPS":   None,          # theory-based — no training set
    "catGRANULE_score": None,         # yeast training — no human LLPS leakage
    "PLAAC_NLLR":      None,          # yeast HMM — no LLPS labels
    "PScore_score":    None,          # biophysical scorer — no LLPS training
    "ESpritz_score":   None,          # disorder predictor — no LLPS training
    "SEG_score":       None,          # algorithmic LCR scorer
    "DeepPhase_score": "DeePhase",
    "SaPS_score":      "SaPS",
    "PdPS_score":      "PdPS",
    "PDL_score":       "PDL",
    "RY_score":        None,             # rule-based composition — no training set
    "ParSe2_score":    None,             # polymer-scaling/helix theory — no training set
    "LLPhyScore_score": "LLPhyScore",
}

masks = leakage_df[["UniProt_ID"]].copy()
masks["clean_global"] = masks["UniProt_ID"].isin(globally_clean)

for score_col, tool in SCORE_TO_TOOL.items():
    if tool is None:
        masks[f"clean_{score_col}"] = True
    else:
        pc = f"pos_{tool}"
        nc = f"neg_{tool}"
        masks[f"clean_{score_col}"] = (leakage_df[pc] == False) & (leakage_df[nc] == False)

masks.to_csv(OUT / "clean_masks.csv", index=False)
print(f"Saved output/clean_masks.csv  ({len(masks.columns)-1} mask columns)")

# ── Print summary ─────────────────────────────────────────────────────────────

print(f"\n{'Tool':<16} {'pos_leaked':>10} {'neg_leaked':>10} {'contra':>7} {'clean':>7}")
print("-" * 55)
for r in summary:
    print(f"{r['Tool']:<16} {r['N_pos_leaked']:>10} {r['N_neg_leaked']:>10} "
          f"{r['N_contradictory']:>7} {r['N_clean']:>7}")

print(f"\nGlobally clean (not in any tool's pos or neg set): {len(globally_clean)} / {len(mem_proteins)}")
if globally_clean:
    print("  " + ", ".join(globally_clean[:15]) + ("..." if len(globally_clean) > 15 else ""))

print("\nPer-tool clean counts (for existing score columns):")
print(f"  {'Score column':<35} {'Tool':<14} {'N_clean':>7}")
print("  " + "-" * 58)
for score_col, tool in SCORE_TO_TOOL.items():
    if tool is None:
        n = len(mem_proteins)
        tool_label = "(theory-based)"
    else:
        row = next(r for r in summary if r["Tool"] == tool)
        n = row["N_clean"]
        tool_label = tool
    print(f"  {score_col:<35} {tool_label:<14} {n:>7}")
