# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

This project benchmarks LLPS (Liquid-Liquid Phase Separation) predictors against
human membrane proteins. It is a data analysis project, not a software library.
The primary entry points are standalone Python scripts that run top-to-bottom
(no argparse, `print` for progress) and write to `output/`.

The old Jupyter notebook / Shiny dashboard workflow lives on the `main` branch.
The `minimal` branch is the active analysis branch.

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .   # only needed for the llps/ package (used in legacy notebooks)
```

Python ≥ 3.11 required.

## Pipeline overview

Run in order on a fresh machine. Steps that need `Predictors_whole_genome_sets/`
(not committed, ~800 MB) are marked with *.

```bash
python build_master_table.py       # builds output/master_table.csv
python wrangle_exp_db.py           # * builds output/predictor_comparison*.csv
python build_leakage_map.py        # builds output/leakage_map.csv
python wrangle.py                  # builds output/full_dataset.csv
python wrangle_background.py       # * builds output/background_scored.csv
python plot.py                     # predictor comparison figures
python plot_roc.py                 # ROC/AUROC/MCC; builds output/clean_auroc_table.csv
python plot_distributions.py       # distribution comparison figures
python rank_consensus.py           # rank-based consensus analysis (Phase 3)
```

All downstream analysis (`plot_roc.py`, `rank_consensus.py`, etc.) reads from the
**committed** output CSVs, so they can be re-run without `Predictors_whole_genome_sets/`.

## Script conventions

- Scripts follow the pattern of `plot.py` / `wrangle.py`: module-level execution,
  `ROOT = Path(__file__).parent`, `print(f"[step N] ...")` for progress.
- All outputs go under `output/`, figures under `output/figures/`.
- Both `output/` and `output/figures/` are committed (except `output/figures/v*/`
  which are versioned ROC run subdirs).
- Do not add argparse. Do not add CLI wrappers. Keep scripts flat and readable.

## Key data facts

- **60 membrane LLPS proteins** — the core benchmark set, from PhasePDB + LLPSDB
- **882 experimental LLPS proteins** (all human) — background for the full proteome ROC
- **18 predictors** — scores in `output/background_scored.csv` and
  `output/predictor_comparison.csv`
- **p(LLPS) = FuzDrop score.** Do not use it as an independent predictor variable
  in analyses that already include FuzDrop.
- **pLLPS_Class (High/Medium/Low)** is an arbitrary cutoff from prior work.
  Use it only as metadata; never as an ordinal feature in statistical analyses.
- **LLPhyScore direction is inverted**: lower raw score = more LLPS-prone.
  In `rank_consensus.py` the sign is flipped before ranking. AUROC is ~0.38
  without flipping, ~0.49 after — this predictor performs poorly on membrane proteins.

## Key output files

| File | What it is |
|---|---|
| `output/master_table.csv` | 882 experimental LLPS proteins, provenance metadata |
| `output/predictor_comparison.csv` | 60 membrane proteins × all predictor scores + metadata |
| `output/background_scored.csv` | Proteome-wide score matrix; primary ROC input |
| `output/clean_auroc_table.csv` | AUROC/MCC per predictor × scenario, with leakage-aware "clean" column |
| `output/leakage_map.csv` | Per-protein positive/negative training-set leakage flags |
| `output/rank_consensus_table.csv` | Per-protein consensus ranks, AUROC-weighted mean, SD |
| `output/rank_consensus_report.md` | Full write-up of consensus analysis methods and findings |

## Score columns in background_scored.csv

```
PICNIC_score, PICNIC_GO_score, PSAP_score, PSPHunter_prob, PSPire_score,
FuzDrop_pLLPS, catGRANULE_score, PLAAC_NLLR, PScore_score, ESpritz_score,
SEG_score, SaPS_score, PdPS_score, DeepPhase_score, PDL_score, RY_score,
ParSe2_score, LLPhyScore_score
```

## What is NOT committed

- `Predictors_whole_genome_sets/` (~800 MB) — whole-genome predictor score files.
  Required only for `wrangle_exp_db.py` and `wrangle_background.py`. Expected
  file paths are documented at the top of each script.
- `external_tools/` (~700 MB) — cloned predictor repos (PICNIC, PSPire, etc.)
- `data/alphafold_pdbs/` — AlphaFold PDB files for topology analysis
- `.venv/` — virtual environment

## Topology analysis

Additional scripts (`run_picnic_pdb_segments.py`, `run_pspire_pdb_segments.py`, etc.)
run predictors on individual topology segments (TM helices, cytoplasmic tails).
These require `external_tools/` and `data/alphafold_pdbs/`. Results are in
`output/topology_scores_master.csv` and `output/figures/topology_*.png`.

## Legacy code

The `llps/` package and `llps_functions.py` shim are from the original notebook
workflow on the `main` branch. They are not used by any script on `minimal`.
Do not add new logic there; it will be out of sync.

`deprecated/` — archived notebooks, do not reference.
