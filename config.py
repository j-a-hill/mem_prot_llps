"""
Paths and constants. This is the ONLY file you should need to edit if things move.

Everything else in this folder imports from here, so there is exactly one place
where a path is written down.

Every path below is RELATIVE to this folder, so the pipeline runs anywhere you
clone it -- no machine-specific paths.
"""

from pathlib import Path

# ---------------------------------------------------------------- where things live

HERE = Path(__file__).resolve().parent

SOURCES_DIR = HERE / "sources"   # the five LLPS databases + UniProt, as shipped
INPUTS = HERE / "inputs"         # precomputed predictor output (see README)

RAW = HERE / "raw"               # 01 writes the source files it will actually use here
BUILD = HERE / "build"           # 02 writes master.csv here
FIGS = HERE / "figures"          # 05 and 08 write PNGs here
TABLES = HERE / "tables"         # 04, 06 and 07 write result tables here

for _d in (RAW, BUILD, FIGS, TABLES):
    _d.mkdir(exist_ok=True)


# ---------------------------------------------------------------- source databases
#
# Each entry says: where to download it from, which shipped copy to fall back on,
# and how to check that a download is real (see 01_download.py -- a download can
# return HTTP 200 at a plausible size and still be junk).
#
# The shipped copies in sources/ are dated snapshots. They are what makes this
# folder self-contained; the URL is there so you can refresh a database when you
# want to. They are NOT all the exact files behind the published 475-protein set --
# the CD-CODE snapshot is a later crawl, which is why the rebuild lands at 472.

SOURCES = {
    "cdcode": {
        "url": "https://cd-code.org/api/proteins",
        "local": SOURCES_DIR / "cdcode_proteins_all.json",
        "out": "cdcode_proteins.json",
        "kind": "json_records",
        "min_records": 4000,
        "note": "CD-CODE, crawled from the public API 2026-07-16",
    },
    "phasepdb": {
        "url": "http://db.phasep.pro/download/",
        "local": SOURCES_DIR / "phasepdb_summary.csv",
        "out": "phasepdb.csv",
        "kind": "table",
        "min_records": 500,
        "note": "PhaSepDB summary table, downloaded 2026-06-15 (membership + MLO types)",
    },
    # PhaSepDB's PS-self / PS-other classification, from its JSON API rather than the
    # summary export. This is a separate fetch because the two carry different things:
    # the export gives membership and MLO types, the API gives the class_ field.
    # See the note on PHASEPDB_CLASS below for why the API is the only correct source.
    "phasepdb_class": {
        "url": "https://db.phasep.pro/api/proteins?limit=500",
        "local": SOURCES_DIR / "phasepdb_api_class.csv",
        "out": "phasepdb_class.csv",
        "kind": "paginated_api",
        "min_records": 3000,
        "note": "PhaSepDB API class_ field, one row per curated PMID entry",
    },
    "drllps": {
        "url": "https://guolab.wchscu.cn/DrLLPS/download/LLPS.txt",
        "local": SOURCES_DIR / "drllps_LLPS.tsv",
        "out": "drllps.tsv",
        "kind": "table",
        "min_records": 5000,
        "note": "DrLLPS bulk LLPS table (tab-separated despite the .csv name)",
    },
    "phasepro": {
        "url": "https://phasepro.elte.hu/download_full.txt",
        "local": SOURCES_DIR / "phasepro.json",
        "out": "phasepro.json",
        "kind": "json_by_key",
        "min_records": 100,
        "note": "PhaSePro full download (JSON keyed by accession)",
    },
    "llpsdb": {
        "url": "http://bio-comp.org.cn/llpsdb/download.html",
        "local": SOURCES_DIR / "protein_LLPSDB.xls",
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
        "local": SOURCES_DIR / "uniprot_tm_cache.csv",
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

# ------------------------------------------------------- precomputed predictor input
#
# Running the 18 LLPS predictors is a separate, slow job (external tools, GPU time,
# and ~780 MB of whole-proteome score files). This pipeline treats their output as an
# INPUT, joined on UniProt accession. The five files below are shipped in inputs/.
#
# They ARE regenerable, but not from anything in this folder -- the chain is
#   Predictors_whole_genome_sets/  (~780 MB of raw per-tool proteome scores)
#     -> wrangle_background.py     -> background_scored.csv
#     -> run_pipeline_475.py       -> predictor_comparison.csv, clean_masks.csv
#     -> rank_consensus.py         -> rank_consensus_table.csv, consensus_features.csv
# and those generators live in the full project, not here. See README for what to
# fetch if you need to rebuild them rather than use the shipped copies.

SCORES = INPUTS / "predictor_comparison.csv"
FEATURES = INPUTS / "consensus_features.csv"
CONSENSUS = INPUTS / "rank_consensus_table.csv"

# The BACKGROUND table: all 20,447 scored human proteins, not just the LLPS set.
# Steps 06 and 08 need this -- you cannot ask "how does this tool score a membrane
# protein relative to a soluble one" from the 472-row table alone.
#   is_membrane_protein  1 for the 5,192 membrane proteins
#   label_mem_bg         1 for the 473 membrane LLPS positives
BACKGROUND = INPUTS / "background_scored.csv"

# Per-tool training-set leakage flags. A protein a tool was TRAINED on is not a fair
# test of that tool. 06_background.py asserts this file has 475 rows: a 427-row copy
# is a superseded set from an earlier stage of the project, and mixing the two caused
# a real erratum here.
CLEAN_MASKS = INPUTS / "clean_masks.csv"

# The canonical 475-protein set, used ONLY to check the rebuild against the paper.
CANONICAL = INPUTS / "membrane_llps_master_475.csv"

# The 18 predictor score columns, as they are named in predictor_comparison.csv.
PREDICTORS = [
    "PICNIC_score", "PICNIC_GO_score", "PSAP_score", "PSPHunter_prob",
    "PSPire_score", "FuzDrop_pLLPS", "catGRANULE_score", "PLAAC_NLLR",
    "PScore_score", "ESpritz_score", "SEG_score", "SaPS_score",
    "PdPS_score", "DeepPhase_score", "PDL_score", "RY_score",
    "ParSe2_score", "LLPhyScore_score",
]

# ------------------------------------------------------------- PhaSepDB class_ field
#
# PS-self  = the protein drives its own phase separation
# PS-other = it partitions into a condensate driven by something else
#
# READ THIS FROM THE API, NEVER BY REGEX ON THE PROSE EXPORT.
#
# class_ is a structured field on PhaSepDB's JSON API and it is populated for every
# one of the 3,528 curated entries (2,489 PS-other / 1,039 PS-self, 1,873 unique
# accessions). The phasepdb_summary export is an AI-written prose summary in which the
# tag only sometimes appears inline, so scraping it finds a small fraction and makes
# the field look sparse. That sparsity is a scraping artefact, not a property of the
# database -- an earlier version of this pipeline made exactly that mistake and
# reported the tag for 12 of our 58 proteins instead of 58 of 58.
#
# class_ is per-PMID-ENTRY, not per-protein, so a protein curated in several papers
# can carry both tags. Four proteins in our set do. Report the literal tags present
# ("PS-self;PS-other") and expose the per-protein entry counts so the balance is
# visible; do not invent a "mixed" category. A protein counts as a PhaSepDB driver
# call if it has at least one PS-self entry.
PHASEPDB_CLASS = SOURCES_DIR / "phasepdb_api_class.csv"

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

# Which training-set tool determines each score column's leakage status. None means
# the tool has no training set (a physics or composition score), so it cannot leak and
# every protein is clean for it by construction. Used by 09_derive_inputs.py.
SCORE_TO_TOOL = {
    "PICNIC_score": "PICNIC",
    "PICNIC_GO_score": "PICNIC",
    "PSAP_score": "PSAP",
    "PSPHunter_prob": "PSPHunter",
    "PSPire_score": "PSPire",
    "FuzDrop_pLLPS": None,
    "catGRANULE_score": None,
    "PLAAC_NLLR": None,
    "PScore_score": None,
    "ESpritz_score": None,
    "SEG_score": None,
    "SaPS_score": "SaPS",
    "PdPS_score": "PdPS",
    "DeepPhase_score": "DeePhase",
    "PDL_score": "PDL",
    "RY_score": None,
    "ParSe2_score": None,
    "LLPhyScore_score": "LLPhyScore",
}

# Two colours used everywhere, colourblind-safe.
C_SINGLE = "#0072B2"   # single-pass
C_MULTI = "#E69F00"    # multi-pass
