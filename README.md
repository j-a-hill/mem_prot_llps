# Membrane Protein LLPS Predictor Comparison

Benchmarking and consensus analysis of LLPS (Liquid-Liquid Phase Separation)
predictors applied to human membrane proteins.

---

## What this project does

18 LLPS predictors are scored against a curated set of 60 human membrane proteins
with experimental LLPS evidence, and evaluated against proteome/membrane backgrounds.
The analysis produces:

- **ROC / AUROC / MCC benchmarks** for each predictor across five comparison scenarios
- **Training-set leakage maps** identifying which benchmark proteins appeared in
  each predictor's own training data
- **Topology-aware scoring** running predictors on individual transmembrane segments
  vs. cytoplasmic tails vs. the full protein
- **Rank-based consensus analysis** identifying which proteins are consistently
  predicted, systematically underpredicted, or show cross-predictor disagreement

---

## Quick start

```bash
# 1. Create and activate virtual environment
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate

# 2. Install dependencies
pip install -r requirements.txt

# 3. Reproduce all outputs from committed data
python build_master_table.py      # output/master_table.csv
python build_leakage_map.py       # output/leakage_map.csv
python wrangle.py                 # output/full_dataset.csv
python wrangle_exp_db.py          # output/predictor_comparison*.csv
python wrangle_background.py      # output/background_scored.csv
python plot.py                    # output/figures/predictor_*.png
python plot_roc.py                # output/figures/roc_*.png, auroc_*.png, mcc_*.png
python plot_distributions.py      # output/figures/distributions_*.png
python rank_consensus.py          # output/rank_consensus_table.csv + figures
```

Python ≥ 3.11 required.

> **Note:** Steps 3–8 can be run with only the committed data — no external tools needed.
> The predictor score files in `Predictors_whole_genome_sets/` (not committed, ~800 MB)
> are only required if you need to re-run `wrangle_exp_db.py` or `wrangle_background.py`
> from scratch. See [§ Reproducing predictor scores](#reproducing-predictor-scores).

---

## Repository layout

```
mem_prot_llps/
│
├── ── Core pipeline scripts ──────────────────────────────────────────────────
│
├── build_master_table.py          Phase 0a: build provenance table for all
│                                  experimental LLPS proteins (PhasePDB + LLPSDB)
├── wrangle_exp_db.py              Phase 0b: combine experimental DBs, join all
│                                  predictor whole-genome scores
├── build_leakage_map.py           Phase 0c: cross-reference benchmark proteins
│                                  against each predictor's training set
├── wrangle.py                     Phase 0d: annotate original dataset with GO,
│                                  location, functional categories
├── wrangle_background.py          Phase 1: build proteome-wide scored dataset
│                                  for ROC analysis (label columns + all scores)
│
├── plot.py                        Phase 2a: predictor comparison plots
├── plot_roc.py                    Phase 2b: ROC/AUROC/MCC/AUPRC analysis;
│                                  outputs clean_auroc_table.csv
├── plot_distributions.py          Phase 2c: score distribution comparisons
│
├── rank_consensus.py              Phase 3: rank-based consensus analysis
│                                  (see output/rank_consensus_report.md)
│
├── ── Topology analysis ──────────────────────────────────────────────────────
│
├── extract_topology_fasta.py      Extract per-segment FASTA sequences
├── build_segment_fastas.py        Build FASTA files for each topology region
├── run_picnic_pdb_segments.py     Run PICNIC on individual PDB segments
├── run_pspire_pdb_segments.py     Run PSPire on individual PDB segments
├── run_picnic_pdb_region.py       Run PICNIC on PDB-derived topology regions
├── run_pspire_pdb_region.py       Run PSPire on PDB-derived topology regions
├── integrate_pdb_segments.py      Merge segment scores into master table
├── integrate_picnic_pdb_region.py Merge PICNIC region scores
├── integrate_pspire_pdb_region.py Merge PSPire region scores
├── build_topology_scores_master.py Consolidate all topology scores
├── extract_topology_scores.py     Extract scores for plotting
├── plot_topology_*.py             Topology comparison figures
│
├── ── Utility scripts ────────────────────────────────────────────────────────
│
├── score_parse2.py                Score proteins with ParSe2
├── score_topology_parse2.py       Score topology segments with ParSe2
├── run_fuzdrop_topology.py        Run FuzDrop on topology segments
├── run_fuzdrop_whole.py           Run FuzDrop on whole proteins
├── run_deephase_topology.py       Run DeepPhase on topology segments
├── run_llphyscore_topology.py     Run LLPhyScore on topology segments
├── run_psphunter.py               Run PSPHunter
├── run_pspspredict_segments.py    Run PSPspredict on segments
├── export_membrane_scores_table.py Export scores as Excel table
├── bin_fasta_for_submission.py    Bin FASTAs for submission to predictors
│
├── ── Reference data ─────────────────────────────────────────────────────────
│
├── LLPS_DB_data/                  Raw database downloads
│   ├── LLPS_Natural_protein.zip   LLPSDB natural protein entries
│   └── phasepdb_summary_database_*.csv  PhasePDB snapshots
├── exp_db/                        Processed experimental database files
│   ├── phasepdb_summary_database_*.csv
│   └── protein_LLPSDB.xls
├── training_sets/                 Predictor training set files (for leakage mapping)
│   ├── PICNIC_S1_training_sets.xlsx
│   ├── PSPire_S4_training_sets.xlsx
│   └── ...
│
├── ── Outputs ────────────────────────────────────────────────────────────────
│
├── output/
│   ├── master_table.csv               882 experimental LLPS proteins, provenance
│   ├── exp_db_combined.csv            Unified PhasePDB+LLPSDB dataset
│   ├── predictor_comparison_all.csv   All 882 proteins × predictor scores
│   ├── predictor_comparison_mem.csv   60 membrane proteins × predictor scores
│   ├── predictor_comparison.csv       60 membrane proteins (with metadata)
│   ├── background_scored.csv          Proteome-wide scores (ROC input)
│   ├── clean_auroc_table.csv          AUROC/MCC/AUPRC per predictor × scenario
│   ├── matched_auroc_table.csv        AUROC with matched background sizes
│   ├── leakage_map.csv                Per-protein training-set leakage flags
│   ├── leakage_summary.csv            Per-predictor leakage counts
│   ├── globally_clean.csv             Proteins absent from all training sets
│   ├── clean_masks.csv                Clean-only masks per predictor
│   ├── rank_consensus_table.csv       Per-protein consensus ranks + summary stats
│   ├── rank_consensus_report.md       Methods + findings for consensus analysis
│   ├── training_set_provenance.md     Which proteins came from which training set
│   ├── predictor_methods_report.md    Predictor algorithm descriptions
│   ├── predictor_data_and_normalisation.md  Score normalisation notes
│   ├── topology_scores_master.csv     Per-region topology scores
│   └── figures/                       All generated figures (committed)
│
├── ── Not committed ──────────────────────────────────────────────────────────
│
├── Predictors_whole_genome_sets/  ~800 MB — whole-genome score files from each
│                                  predictor. Required only to re-run wrangle_exp_db.py
│                                  or wrangle_background.py from scratch.
├── external_tools/                ~700 MB — cloned predictor repos (PICNIC, PSPire, etc.)
├── .venv/                         Python virtual environment
└── data/alphafold_pdbs/           Downloaded AlphaFold PDB files (topology analysis)
```

---

## Pipeline in detail

### Phase 0 — Build reference tables (run once)

| Script | Input | Output | Notes |
|---|---|---|---|
| `build_master_table.py` | `LLPS_DB_data/`, `exp_db/` | `output/master_table.csv` | Parses PhasePDB + LLPSDB; one row per UniProt accession |
| `wrangle_exp_db.py` | `exp_db/`, `Predictors_whole_genome_sets/` | `output/exp_db_combined.csv`, `output/predictor_comparison_all.csv`, `output/predictor_comparison_mem.csv` | Joins all predictor whole-genome scores onto the 882 experimental proteins |
| `build_leakage_map.py` | `training_sets/`, `output/master_table.csv` | `output/leakage_map.csv`, `output/leakage_summary.csv`, `output/globally_clean.csv` | Flags which benchmark proteins were in each predictor's positive/negative training set |
| `wrangle.py` | `Human Phase separation data.xlsx` | `output/full_dataset.csv` | Annotates original dataset with UniProt GO terms, location, functional categories |
| `wrangle_background.py` | `Predictors_whole_genome_sets/`, `output/exp_db_combined.csv` | `output/background_scored.csv` | Proteome-wide score matrix with label columns for ROC analysis |

### Phase 1 — Scoring (requires external tools)

These scripts invoke the predictor tools directly and are only needed if re-running
from scratch. Outputs land in `Predictors_whole_genome_sets/` (not committed).

See `output/predictor_methods_report.md` for how each predictor was invoked
and `output/predictor_data_and_normalisation.md` for score normalisation.

### Phase 2 — Analysis

| Script | Input | Output |
|---|---|---|
| `plot.py` | `output/predictor_comparison*.csv`, `output/background_scored.csv` | `output/figures/predictor_*.png` |
| `plot_roc.py` | `output/background_scored.csv`, `output/leakage_map.csv` | `output/figures/roc_*.png`, `output/clean_auroc_table.csv`, `output/matched_auroc_table.csv` |
| `plot_distributions.py` | `output/background_scored.csv` | `output/figures/distributions_*.png` |

### Phase 3 — Rank consensus

```bash
python rank_consensus.py
```

Reads `output/predictor_comparison.csv`, `output/background_scored.csv`,
`output/clean_auroc_table.csv`, and `output/leakage_map.csv`.

Outputs:
- `output/rank_consensus_table.csv` — 60 proteins × normalised ranks + AUROC-weighted mean, SD
- `output/rank_consensus_report.md` — full methods, key findings, caveats
- `output/figures/rank_consensus_spearman.png` — predictor pairwise Spearman correlation heatmap
- `output/figures/rank_consensus_heatmap.png` — protein × predictor rank heatmap
- `output/figures/rank_consensus_scatter.png` — consensus vs disagreement scatter
- `output/figures/rank_consensus_comparison.png` — unweighted vs weighted rank comparison
- `output/figures/rank_consensus_features.png` — feature correlations (length, TMD count)

See `output/rank_consensus_report.md` for a full write-up.

---

## Predictors included

| Predictor | Score column | Notes |
|---|---|---|
| PICNIC | `PICNIC_score` | Structure-aware |
| PICNIC (GO) | `PICNIC_GO_score` | PICNIC + GO term enrichment |
| PSPire | `PSPire_score` | ML, IDR-aware |
| PSPHunter | `PSPHunter_prob` | ML |
| FuzDrop | `FuzDrop_pLLPS` | Biophysics-based; **p(LLPS) in the source dataset is the FuzDrop score** |
| catGRANULE | `catGRANULE_score` | Granule propensity |
| PLAAC | `PLAAC_NLLR` | Prion-like domain log-likelihood ratio |
| PScore | `PScore_score` | Cation-π interaction propensity |
| ESpritz | `ESpritz_score` | Intrinsic disorder fraction |
| SEG | `SEG_score` | Low-complexity sequence fraction |
| SaPS | `SaPS_score` | Sequence-based |
| PdPS | `PdPS_score` | Sequence-based |
| DeepPhase | `DeepPhase_score` | Deep learning |
| PDL | `PDL_score` | Phase diagram-based |
| R+Y | `RY_score` | Arg + Tyr composition |
| ParSe2 | `ParSe2_score` | Physicochemical |
| LLPhyScore | `LLPhyScore_score` | Physical free-energy score; **lower = more LLPS-prone** |
| PSAP | `PSAP_score` | Sequence-based ML |

All scores are available in `output/background_scored.csv` in a unified format.
For LLPhyScore the sign is flipped in `rank_consensus.py` before ranking.

---

## Reproducing predictor scores

`wrangle_exp_db.py` and `wrangle_background.py` read raw score files from
`Predictors_whole_genome_sets/`. This directory is not committed (~800 MB).

To reproduce it:
1. Download whole-proteome score files from each predictor's website/Zenodo.
2. Place them under `Predictors_whole_genome_sets/<predictor_name>/`.
3. The expected file paths are documented at the top of `wrangle_exp_db.py`.

Alternatively, work directly from the committed `output/predictor_comparison*.csv`
and `output/background_scored.csv` — all downstream analysis steps read from these.

---

## Key caveats

- **p(LLPS) is the FuzDrop score.** Do not use it as an independent variable in
  analyses that already include FuzDrop as a predictor.
- **pLLPS_Class (High/Medium/Low)** is an arbitrary cutoff from a prior version
  of the project. Treat it as metadata only; do not use as an ordinal predictor.
- **Training-set leakage** affects PICNIC, PSAP, PSPHunter, PSPire, SaPS, PdPS,
  LLPhyScore, and DeepPhase for at least some benchmark proteins. The `leakage_any`
  column in `rank_consensus_table.csv` flags affected proteins.
- **n = 60** membrane LLPS proteins limits statistical power. Feature correlations
  with BH-adjusted p < 0.05 should be treated as hypothesis-generating.

---

## Environment

Python 3.11+. All dependencies in `requirements.txt`.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .   # if using the llps/ package
```
