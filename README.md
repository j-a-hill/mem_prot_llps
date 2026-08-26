# minimal/ — the membrane-protein LLPS analysis in eight short scripts

Self-contained. Clone it, run one command, get every table and figure. Plain pandas
and matplotlib, no package to import, no framework, no machine-specific paths.

Covers five of the six main figures in the current paper arc
(`manuscript/FIGURE_ARC_REVIEW.md`, 25 Aug). The one it does not is the evidence
ladder — see *Not covered* at the bottom.

## Layout

```
minimal/
  0*.py, config.py    the pipeline (about 1,500 lines total)
  sources/    14 MB   the five LLPS databases + UniProt, as dated snapshots
  inputs/      5 MB   precomputed predictor output (see "Regenerating inputs/")
  raw/                written by 01 -- the source files it actually used
  build/              written by 02 -- master.csv
  tables/             written by 04, 06, 07
  figures/            written by 05, 08
```

`sources/` and `inputs/` are shipped so the pipeline runs on clone. Everything else
is generated.

## Run it

```bash
conda activate memllps
cd minimal
python run_all.py            # everything, about a minute
```

Or one step at a time, which is the point of splitting them up:

```bash
python 01_download.py            # get the database files -> raw/
python 02_build_master.py        # wrangle them into one table -> build/master.csv
python 03_explore.py             # look at the table before trusting it
python 04_stats.py               # topology / hydropathy tests
python 05_figures.py             # figures for 04
python 06_background.py          # the offset and the background flip
python 07_thresholds.py          # published cutoffs and the shortlist
python 08_figures_background.py  # figures for 06 and 07
```

## The files

| file | what it does |
|---|---|
| `config.py` | every path and constant. **The only file to edit if things move.** |
| `01_download.py` | fetches the six source files into `raw/`, writes `raw/PROVENANCE.csv` |
| `02_build_master.py` | parses the five databases, applies the membrane filter, joins predictor scores → `build/master.csv` |
| `03_explore.py` | prints composition, coverage and missingness. No results, just checks |
| `04_stats.py` | topology vs hydropathy, per-predictor breakdown → `tables/` |
| `05_figures.py` | figures 1–4 |
| `06_background.py` | the membrane score offset and what it does to AUROC → `tables/` |
| `07_thresholds.py` | pass rate per published cutoff, agreement shortlists → `tables/` |
| `08_figures_background.py` | figures 5–8 |
| `run_all.py` | runs 01–08 in order |

## Which figure is which paper figure

The arc review reordered the paper: the score offset opens it, and the topology
result became Figure 4 with its robustness work in supplement. Mapping:

| paper | minimal | file |
|---|---|---|
| Fig 1 | Fig 5 | `fig5_membrane_offset.png` |
| Fig 2 | Fig 6 | `fig6_background_flip.png` |
| Fig 3 | — | feature correlations, not rebuilt (see below) |
| Fig 4 | Fig 3 | `fig3_hydropathy_confound.png` |
| Fig 5 | — | evidence ladder, not rebuilt (see below) |
| Fig 6 | Figs 7–8 | `fig7_pass_rates.png`, `fig8_cutoff_agreement.png` |
| Fig S1 | Fig 1 | `fig1_database_coverage.png` |
| Fig S6 | Figs 2, 4 | `fig2_topology_consensus.png`, `fig4_per_predictor_topology.png` |

The minimal numbering is just run order; it is not the paper's.

## Where the data comes from

Five LLPS databases and one UniProt annotation file, all shipped in `sources/`:

| source | what it gives | file in `sources/` |
|---|---|---|
| CD-CODE | condensate membership, role words (driver/member) | `cdcode_proteins_all.json` + 3 role crawls |
| PhaSepDB | membership, PS-self / PS-other class | `phasepdb_summary.csv` |
| DrLLPS | membership, Scaffold / Regulator / Client type | `drllps_LLPS.tsv` |
| PhaSePro | membership, partner dependence | `phasepro.json` |
| LLPSDB | membership only (no role field) | `protein_LLPSDB.xls` |
| UniProt | TRANSMEM / INTRAMEM features → the membrane filter | `uniprot_tm_cache.csv` |

`01_download.py` tries each URL, then **parses what came back** and keeps it only if
it contains a plausible number of records. Anything else falls back to the shipped
snapshot, with the reason recorded in `raw/PROVENANCE.csv`.

The parse check is the point of the script, not a nicety. A download can return HTTP
200 at a plausible size and still be useless, and all three failure modes are live:

- **PhaSepDB** serves an HTML page from its download URL, not a table.
- **CD-CODE**'s API is paginated 25 records at a time, so one GET returns 25 of
  11,144 — a perfectly valid JSON file that would silently shrink the set by 99%.
- **UniProt** serves TSV while the snapshot is CSV.

A size check catches none of these. On a typical run only UniProt downloads cleanly
(20,431 rows, slightly ahead of the 20,366-row snapshot); the rest fall back. The set
size is unchanged either way, which is the useful thing to know.

The snapshots are stable — they do not shift under you the way a live download does.
They are **not** all the exact files behind the published 475-protein set: the CD-CODE
snapshot is a later crawl than the one used then. That is the source of the 472-vs-475
difference below, and `02_build_master.py` reports it rather than hiding it.

## Regenerating `inputs/`

The five files in `inputs/` are the output of the 18 LLPS predictors, treated here as
an input and joined on UniProt accession.

| file | rows | why it is needed |
|---|---|---|
| `predictor_comparison.csv` | 475 | the 18 scores for the LLPS set |
| `background_scored.csv` | 20,447 | **all** scored human proteins — the only way to ask how a tool scores a membrane protein relative to a soluble one. Steps 6 and 8 need this |
| `consensus_features.csv` | 475 | sequence features (hydropathy, length, composition) |
| `rank_consensus_table.csv` | 475 | the cross-tool consensus rank |
| `clean_masks.csv` | 475 | per-tool training-leakage flags |

They are regenerable, but **not from anything in this folder** — the chain starts from
~780 MB of raw per-tool proteome scores that are not in the repo:

```
Predictors_whole_genome_sets/        raw per-tool scores, ~780 MB, gitignored
  -> wrangle_background.py           -> background_scored.csv
  -> run_pipeline_475.py             -> predictor_comparison.csv, clean_masks.csv
  -> rank_consensus.py               -> rank_consensus_table.csv, consensus_features.csv
```

Two things to know before trying it. `Predictors_whole_genome_sets/` holds the outputs
of local runs of each predictor (PICNIC, PSAP, PSPHunter, PSPire, FuzDrop, DeePhase,
PDL and the rest) — some of those tools need a GPU, and a few are not redistributable,
so this step is a re-run of the tools, not a download. And the 475-era generators live
in the full project's `v2.0/` directory, which is not on this branch; the tracked
root-level `rank_consensus.py` has diverged from the `v2.0` one by ~45 lines, so use
the `v2.0` copy if you are reproducing the shipped tables exactly.

For editing the analysis — which is what this folder is for — use the shipped copies.

## How the set is built
## How the set is built

```
5 databases, human entries only
        |  union of accessions
     4,862 human proteins any database calls phase-separating
        |  keep those with >=1 UniProt TRANSMEM or INTRAMEM feature
       472 membrane LLPS proteins   <-- the study set
```

`02_build_master.py` prints the count at every stage, so the set size is derived in
front of you rather than asserted.

**This gives 472, against the published 475.** The script names all three misses and
why: two (`O15533`, `Q01726`) are in no current database snapshot — the published set
used an earlier CD-CODE crawl that still listed them; one (`I3L0A0`) is an unreviewed
UniProt entry, and the TM annotation file covers reviewed entries only. Nothing is in
the rebuild that is not in the published set. If you need exact paper numbers, read
`config.CANONICAL` instead of rebuilding.

## Things that will bite you

These are real traps in this data, each handled in the code with a comment:

- **LLPhyScore runs backwards.** Lower = more likely to phase separate, opposite to
  every other tool. `config.INVERTED` lists it and `04_stats.py` flips the sign.
  Published AUROC tables report it *unflipped*, so turn the flip off if reproducing
  one of those. Without the flip it looks like the lone tool favouring single-pass
  proteins; with it, it is the lone tool favouring multi-pass.
- **DrLLPS `LLPS.csv` is tab-separated** despite the `.csv` name, and one line has a
  stray extra field. `02_build_master.py` reads every delimited file through a
  `read_table()` helper that sniffs the separator, so the same code handles the
  snapshot and a live download without a flag.
- **Role words are not comparable across databases.** CD-CODE `member`, DrLLPS
  `Client` and PhaSepDB `PS-other` are different claims from different evidence.
  `master.csv` keeps each database's word verbatim in its own column and never
  merges them into a shared category.
- **PhaSepDB class tags are sparse in this snapshot** — present for 12 of our 58
  PhaSepDB proteins, because the tag lives in a prose summary column. Missing tag
  means *unknown*, not *not self-driving*.
- **Predictor coverage is incomplete** for PSAP (441/472), PScore (452) and
  DeepPhase (466) — they fail on short or very long sequences. The failures are
  length-dependent, not random, so a per-tool number computed on its own available
  subset is not strictly comparable to another tool's. `03_explore.py` prints the
  coverage per tool.
- **Intramembrane-only proteins (12) span no membrane**, so single-vs-multi-pass is
  undefined for them. `04_stats.py` and `05_figures.py` drop them.

- **The membrane-segment count has two definitions in this project and they
  disagree for ~8% of the set.** `n_transmem` (used here, and by the canonical
  master's `topology` column) counts TRANSMEM features only; `TMD_count` in
  `rank_consensus_table.csv` counts TRANSMEM **and** INTRAMEM. It bites on the
  shortlist: TGO1 has one of each, so the >=75% group is 15/15 single-pass by the
  first count and 14/15 by the second. Both are defensible — say which you used.

## What the analysis finds

**The headline (steps 6 and 8).** A tool's score alone separates membrane from soluble
proteins even among proteins *no database calls phase-separating*: 16 of 18 tools score
membrane proteins below soluble ones, PICNIC (GO) at AUROC 0.15. That standing penalty
is what decides the verdict. Score the same 473 positives against a soluble background
and then a membrane one, and 12 of 18 tools cross from below chance to above it — the
tools did not change, only the comparison did. And how much a tool gains from the
membrane background is almost entirely predicted by how much it penalises membrane
proteins in the first place: Spearman ρ = −0.99. This is the paper's cold open.

**Published cutoffs disagree wildly (step 7).** On identical proteins, pass rates run
from 4.5% (PLAAC) to 71.1% (catGRANULE) — a 16-fold spread across 10 tools nominally
answering one question. The median protein clears 30% of its cutoffs. Requiring >=75%
agreement gives 15 proteins; >=50% gives 129 — reported as a strict/permissive pair,
not one cutoff pretending to be definitive. The >=75% group separates on consensus rank
(0.79 vs 0.48, P = 1.2e-9) but **not** on cross-tool disagreement (P = 0.09): being
high-confidence is a rank phenomenon, and the tools argue about the top proteins as
much as about the rest.

**The topology strand (steps 4 and 5).**

1. **Predictors rank single-pass membrane proteins above multi-pass ones.** Median
   consensus rank 0.55 vs 0.44, rank-biserial +0.39 (95% CI +0.29 to +0.48),
   P = 1.2e-12. Spreads are equal (Fligner-Killeen P = 0.65) so this reads as a
   location shift, and the spread-free Brunner-Munzel test agrees.
2. **Hydropathy, not topology, is doing the work.** Within hydropathy quintiles the
   gap disappears — single-pass ranks higher in only 1 of 5 bins. Partial Spearman
   goes from +0.33 unadjusted to −0.09 with hydropathy held constant. Figure 3 shows
   why: both classes sit on one downward hydropathy–rank relationship and differ
   only in where along it they sit. Caveat: the extreme bins are badly unbalanced
   (Q1 has 82 single-pass vs 10 multi-pass, Q5 the reverse), which limits how much
   weight the stratified test can carry.
3. **12 of 18 tools show the effect individually**, 5 show nothing, and LLPhyScore
   reverses it (ranking multi-pass higher, rank-biserial −0.34). A property shared by
   most tools but not all is easier to read as a feature of the predictors than of
   membrane biology — and LLPhyScore is a useful counterexample, since it is the one
   tool whose behaviour a purely hydropathy-driven account does not predict.

## Figures

Following this project's conventions: no titles on the canvas (the claim goes in the
caption you write), no notes about how the figure was built, database provenance in
the axis label, and a colour key only where the mark alone does not say it.

- `fig1_database_coverage.png` — per-database contribution, split shared vs unique
- `fig2_topology_consensus.png` — consensus rank by topology class
- `fig3_hydropathy_confound.png` — the confound in one panel
- `fig4_per_predictor_topology.png` — the effect tool by tool, with CIs
- `fig5_membrane_offset.png` — the offset per tool, and offset vs the AUROC it buys
- `fig6_background_flip.png` — same positives, two backgrounds, one line per tool
- `fig7_pass_rates.png` — pass rate at each published cutoff
- `fig8_cutoff_agreement.png` — cutoffs cleared per protein, and the >=75% group on
  the consensus axes

## Not covered

Two things in the current paper are deliberately absent from this folder:

- **Feature correlations (paper Fig 3).** The 14-feature Spearman panel. Not rebuilt
  because the published version uses an AUROC-weighted rank target and a BH correction
  across both protein sets — reproducing it faithfully means importing machinery that
  would roughly double this pipeline.
- **The evidence ladder (paper Fig 5).** Local LR+ by the Pejaver method — adaptive
  sliding window, bootstrap 95% lower bound, monotone threshold rule. Left out by
  choice: it is the one piece that genuinely fights "simple enough to follow", and it
  needs its own careful treatment rather than a compressed version.

Both exist as verified tables in the full project (`positron/data/tables/`, files
`feature_correlations_both_sets.csv` and
`evidence_calibration_pejaver_method.csv`) — not on this branch. Drop either into
`inputs/` and add a `09_*.py` if you want to plot them without re-deriving.

## Adding your own analysis

`build/master.csv` is one row per protein, 69 columns: membership booleans per
database, each database's role word, `n_transmem`, `topology`, 18 predictor scores,
sequence features, `mean_rank`. Start a `06_myanalysis.py` with

```python
import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C, pandas as pd
m = pd.read_csv(C.BUILD / "master.csv")
```

and remember the LLPhyScore sign flip if you touch that column.
