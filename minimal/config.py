"""
Paths and constants. This is the ONLY file you should need to edit if things move.

Everything else in this folder imports from here, so there is exactly one place
where a path is written down.
"""

from pathlib import Path

# ---------------------------------------------------------------- where things live

# The project root (the folder that contains v2.0/, exp_db/, data/ ...).
PROJECT = Path("/home/jake/mem_prot_llps/mem_prot_llps")

# One level up, where the catGRANULE checkouts sit (holds the DrLLPS table).
PROJECT_PARENT = PROJECT.parent

# This minimal pipeline's own folders. Created automatically if missing.
HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"          # 01 writes the source database files here
BUILD = HERE / "build"      # 02 writes master.csv here
FIGS = HERE / "figures"     # 05 writes PNGs here
TABLES = HERE / "tables"    # 03 and 04 write result tables here

for _d in (RAW, BUILD, FIGS, TABLES):
    _d.mkdir(exist_ok=True)


# ---------------------------------------------------------------- source databases
#
# Each entry says: where to download it from, and which local file to fall back on
# if the download does not work. 01_download.py reads this dict and nothing else.
#
# The local files are dated snapshots kept in the project. They are the reproducible
# option in the sense that they do not change under you -- a live download gives a
# different answer every month as the databases grow. But they are NOT all the exact
# files behind the published 475-protein set: the CD-CODE snapshot here is a later
# crawl than the one used then, and two proteins in the published set are no longer
# listed in it. 02_build_master.py reports that gap rather than hiding it. The URL is
# there so you can refresh a database deliberately when you want to.

SOURCES = {
    "cdcode": {
        "url": "https://cd-code.org/api/proteins",
        "local": PROJECT / "output/cdcode/cdcode_proteins_all.json",
        "out": "cdcode_proteins.json",
        "kind": "json_records",
        "min_records": 4000,
        "note": "CD-CODE, crawled from the public API 2026-07-16",
    },
    "phasepdb": {
        "url": "http://db.phasep.pro/download/",
        "local": PROJECT / "exp_db/phasepdb_summary_database_2026-06-15.csv",
        "out": "phasepdb.csv",
        "kind": "table",
        "min_records": 500,
        "note": "PhaSepDB summary table, downloaded 2026-06-15",
    },
    "drllps": {
        "url": "https://guolab.wchscu.cn/DrLLPS/download/LLPS.txt",
        "local": PROJECT_PARENT / "catGRANULE2.0_current/Condensates_Analysis/DrLLPS/LLPS.csv",
        "out": "drllps.tsv",
        "kind": "table",
        "min_records": 5000,
        "note": "DrLLPS bulk LLPS table (tab-separated despite the .csv name)",
    },
    "phasepro": {
        "url": "https://phasepro.elte.hu/download_full.txt",
        "local": PROJECT / "LLPS_DB_data/phasepro.txt",
        "out": "phasepro.json",
        "kind": "json_by_key",
        "min_records": 100,
        "note": "PhaSePro full download (JSON keyed by accession)",
    },
    "llpsdb": {
        "url": "http://bio-comp.org.cn/llpsdb/download.html",
        "local": PROJECT / "exp_db/protein_LLPSDB.xls",
        "out": "llpsdb.xls",
        "kind": "excel",
        "min_records": 150,
        "note": "LLPSDB natural-protein table",
    },
    # Not an LLPS database -- this is the UniProt membrane annotation that decides
    # which proteins count as membrane proteins.
    "uniprot_tm": {
        "url": ("https://rest.uniprot.org/uniprotkb/stream"
                "?query=reviewed:true+AND+organism_id:9606"
                "&fields=accession,ft_transmem,ft_intramem&format=tsv"),
        "local": PROJECT / "data/uniprot_tm_cache.csv",
        "out": "uniprot_tm.csv",
        "kind": "table",
        "min_records": 15000,
        "note": "UniProt human TRANSMEM / INTRAMEM feature strings",
    },
}


# ------------------------------------------------------- precomputed predictor input
#
# Running the 18 LLPS predictors is a separate, slow job (external tools, GPU time).
# This pipeline treats their output as an INPUT, joined on UniProt accession.

_RUN = PROJECT / "v2.0/output/figures/union475_v3/output"

SCORES = _RUN / "predictor_comparison.csv"
FEATURES = _RUN / "consensus_features.csv"
CONSENSUS = _RUN / "rank_consensus_table.csv"

# The BACKGROUND table: all 20,447 scored human proteins, not just the LLPS set.
# This is what Figures 1 and 2 need -- you cannot ask "how does this tool score a
# membrane protein relative to a soluble one" from the 472-row table alone.
#   is_membrane_protein  1 for the 5,192 membrane proteins
#   label_mem_bg         1 for the 473 membrane LLPS positives
BACKGROUND = _RUN / "background_scored.csv"

# Per-tool training-set leakage flags. A protein a tool was TRAINED on is not a fair
# test of that tool, so AUROC is reported both ways -- all positives, and clean only.
CLEAN_MASKS = _RUN / "clean_masks.csv"

# IMPORTANT: always read these from the per-run directory above, never from the
# top-level v2.0/output/. The top-level copies are legacy 427-protein files. Fast
# check: the canonical clean-mask table has 475 rows, the stale one 427.

# The canonical 475-protein set, used ONLY to check the rebuild against the paper.
CANONICAL = PROJECT / "v2.0/output/membrane_llps_master_475.csv"

# The 18 predictor score columns, as they are named in predictor_comparison.csv.
PREDICTORS = [
    "PICNIC_score", "PICNIC_GO_score", "PSAP_score", "PSPHunter_prob",
    "PSPire_score", "FuzDrop_pLLPS", "catGRANULE_score", "PLAAC_NLLR",
    "PScore_score", "ESpritz_score", "SEG_score", "SaPS_score",
    "PdPS_score", "DeepPhase_score", "PDL_score", "RY_score",
    "ParSe2_score", "LLPhyScore_score",
]

# PUBLISHED CUTOFFS -- the score above which the tool's own authors call a protein
# positive. Only 10 of the 18 have a defensible one; the other 8 are left out rather
# than given a made-up 0.5, because:
#   PSAP        its proteome range is compressed to [0, 0.14], so 0.5 is unreachable
#   ESpritz/SEG disorder and low-complexity proxies, not LLPS classifiers
#   SaPS/PdPS   PhaSePred publishes a proteome rank, not a calibrated probability
#   R+Y         a composition ratio
#   ParSe2      a region-level classifier with no whole-protein cutoff
#   LLPhyScore  unbounded and sign-flipped
PUBLISHED_CUTOFF = {
    "catGRANULE_score": 0.0,     # author
    "PLAAC_NLLR": 0.0,           # author
    "PScore_score": 4.0,         # author
    "FuzDrop_pLLPS": 0.6,        # author
    "PICNIC_score": 0.5,         # default 0.5 probability operating point
    "PICNIC_GO_score": 0.5,
    "PSPHunter_prob": 0.5,
    "PSPire_score": 0.5,
    "DeepPhase_score": 0.5,
    "PDL_score": 0.5,
}

# Short names for axis labels.
NICE = {
    "PICNIC_score": "PICNIC", "PICNIC_GO_score": "PICNIC (GO)",
    "PSAP_score": "PSAP", "PSPHunter_prob": "PSPHunter",
    "PSPire_score": "PSPire", "FuzDrop_pLLPS": "FuzDrop",
    "catGRANULE_score": "catGRANULE", "PLAAC_NLLR": "PLAAC",
    "PScore_score": "PScore", "ESpritz_score": "ESpritz",
    "SEG_score": "SEG", "SaPS_score": "SaPS", "PdPS_score": "PdPS",
    "DeepPhase_score": "DeepPhase", "PDL_score": "PDL", "RY_score": "R+Y",
    "ParSe2_score": "ParSe2", "LLPhyScore_score": "LLPhyScore",
}

# GOTCHA -- LLPhyScore runs BACKWARDS. For every other tool a higher score means
# more likely to phase separate; for LLPhyScore a LOWER score does. If you compare
# it with the others without flipping the sign you will read its direction backwards.
# 04_stats.py flips it. Published AUROC tables report it UNFLIPPED, so if you are
# reproducing one of those, turn this off.
INVERTED = ["LLPhyScore_score"]

# Two colours used everywhere, colourblind-safe.
C_SINGLE = "#0072B2"   # single-pass
C_MULTI = "#E69F00"    # multi-pass
