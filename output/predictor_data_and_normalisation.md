# PSPspredict Sub-scores: Raw Data Availability & Normalisation Options

Date: 2026-06-16

---

## Part 1 — Raw Data for PSPspredict Sub-scores

The PSPspredict compilation file contains percentile-normalised scores. Raw outputs for most constituent tools are accessible via downloads or local CLI runs.

| Tool | Raw scores available? | Route | Notes |
|---|---|---|---|
| **PDL** | Yes — download | Supplementary S2, *Briefings in Bioinformatics* (Li et al. 2025) | Same paper as PSPspredict; S2 = raw probabilities before rank transform |
| **catGRANULE** | Yes — download + CLI | Zenodo (zenodo.org/doi/10.5281/zenodo.14205831); GitHub `tartaglialabIIT/catGRANULE2.0`; Suppl. Table S3 of Agostini et al. *Genome Biology* 2025 | **Caveat: this is v2.0** (adds AlphaFold2 features); PSPspredict used v1.0. Scores are not directly comparable — v1 web server still available at tartaglialab.com/catGRANULES for custom runs |
| **PLAAC** | Run locally | Java CLI at `github.com/whitehead/plaac` (MIT); accepts multi-FASTA | No pre-computed proteome download. Straightforward to run on 882 sequences |
| **PScore** | Yes — download | *eLife* 2018 Fig 5 source data 2: Excel file with human proteome scores. Python script also in supplementary | Direct download, raw log-likelihood-style continuous scores |
| **ESpritz-DisProt** | Yes — via D2P2 | `d2p2.pro` — pre-computed for 1,765 proteomes including human, downloadable flat files | Fastest route; no need to run locally. Old ESpritz web server (`old.protein.bio.unipd.it/espritz`) also accepts multi-FASTA upload |
| **SEG** | Run locally (trivial) | `apt install ncbi-seg` or `segmasker` (BLAST+); parse output for per-protein masked fraction | Runs in minutes on full proteome. No licence issues |
| **DeePhase** | Run locally | GitHub `kadiliissaar/DeePhase`, `PREDICT/` directory (CC-BY-NC, academic use) | No pre-computed proteome file in repo or paper. Python CLI on local FASTA |
| **SaPS-10fea / PdPS-10fea** | Likely — check download page | `predict.phasep.pro/download` — stated to host scores for 116,806 sequences across 18 species | Server was unreachable during lookup. Need to verify whether download exposes raw scores or binary calls |

**Priority actions if you want raw sub-scores:**
1. Download PScore human proteome file from eLife directly — already done for you in the paper.
2. Download D2P2 ESpritz-DisProt human proteome flat file.
3. Get PDL Supplementary S2 from the Briefings in Bioinformatics paper.
4. Check predict.phasep.pro/download for SaPS/PdPS raw scores.
5. Run PLAAC locally via JAR (simplest CLI of the batch).
6. Run SEG locally via segmasker — 10 minutes of work.
7. DeePhase and catGRANULE 2.0 require local Python environment setup.

---

## Part 2 — Normalisation Options

### What the field actually does

**The honest answer: the LLPS benchmarking literature doesn't normalise scores before cross-predictor comparison.** Every published benchmark (Saar et al. 2021; Briefings in Bioinformatics 2023 evaluation of 11 predictors; CSBJ 2026 systematic evaluation of 9 predictors) computes AUROC and AUPRC directly from raw scores. This is valid because **AUROC is a rank statistic** — any monotone transformation of scores leaves the ROC curve and AUC unchanged. PLAAC's log-likelihood ratios and PSAP's [0, 0.14] range produce identical AUROC values to their percentile-ranked equivalents.

The field sidesteps the normalisation problem by only ever comparing predictors via rank-based metrics, not composite scores or side-by-side visualisations.

---

### Options and verdicts

#### A. No normalisation (use raw scores)
- **When to use:** ROC curves, AUROC, AUPRC, Spearman correlations
- **Verdict: correct choice for all rank-based analyses.** Spearman correlation on raw scores is mathematically identical to Pearson on percentile-ranked scores. Never transform scores before plotting ROC curves.

---

#### B. Percentile / rank normalisation
Formula: `x' = rank(x) / n` (or mid-rank `(rank - 0.5) / n` to avoid 0 and 1)

- **When to use:** Composite/ensemble scoring; side-by-side visualisation; any situation where you need scores on a common scale
- **Pros:** Robust to outliers; handles PSAP's compressed range and PLAAC's unbounded LLR without distortion; consistent with what PSPspredict already did
- **Cons:** Destroys absolute probability meaning; sensitive to composition of the reference set; ties require handling (use average-rank)
- **Reference set matters:** Rank within the 882-protein experimental set (not just the 60 membrane proteins — too small and too narrow)
- **PSPspredict sub-scores already are proteome-wide percentile ranks** — do not re-rank them within your 882 subset, as this changes the reference and breaks comparability to any published thresholds. Use them as-is on the [0,1] scale.
- **Verdict: recommended for composite scoring and correlation heatmaps.**

---

#### C. Min-max normalisation
Formula: `x' = (x - min) / (max - min)`

- **Pros:** Intuitive [0,1] output; simple
- **Cons:** Highly sensitive to outliers (one extreme protein anchors the scale for everything else); the min/max in your 882-protein experimental set are not biologically meaningful anchors; produces a false impression of calibration
- **Especially problematic for PSAP** (max ~0.008 in your membrane set, ~0.14 in full PSAP output) and **PLAAC** (LLR has no natural maximum — one long Q/N repeat domain can dominate the range)
- **Verdict: not recommended.**

---

#### D. Z-score standardisation
Formula: `x' = (x - μ) / σ`

- **Pros:** Handles different variances; scores in units of standard deviations
- **Cons:** Assumes distributions are conceptually comparable — they are not (PSAP's compression is mechanistic, not a calibration error; z-scoring does not fix it). Output is unbounded, unsuitable for averaging. Strongly influenced by the composition of your test set
- **Verdict: acceptable only for Pearson correlations between predictors** (standard practice in any correlation analysis). Not suitable for composite scoring.

---

#### E. Calibration (Platt scaling / isotonic regression)
Maps raw scores to true posterior probabilities P(LLPS | score) using labelled data.

- **Pros:** If successful, scores from different predictors actually mean the same thing. Addresses PSAP's compression mechanistically (the random forest was trained on class-imbalanced data, biasing it toward the majority class — a sigmoid fit corrects this)
- **Cons:** Requires labelled positives/negatives. With n=60 membrane proteins, isotonic regression (many free parameters) will overfit. Platt scaling (sigmoid, 2 parameters per predictor) is borderline feasible but requires leave-one-out or 5-fold cross-validation within the 882 proteins
- **Verdict: worth attempting if the downstream goal is a combined probability score** with proper cross-validation on the 882 proteins. Not appropriate as a first step; overkill for ROC comparison.

---

### Special considerations for your dataset

**Small membrane subset (n=60):** Min-max and isotonic regression are inadvisable. Bootstrap all correlation estimates to get confidence intervals — Spearman correlations on n=60 have wide CIs. Use AUPRC alongside AUROC (more informative under class imbalance).

**Domain shift:** All predictors were trained on IDR-rich soluble proteins. Membrane proteins are systematically different. Any normalisation *within* the 60 membrane proteins corrects for this implicitly (by re-ranking among membrane proteins only), but breaks comparability to the full experimental set. Decide whether you want to characterise predictor performance in general, or specifically for membrane proteins — these require different reference sets.

**PSPspredict sub-scores are already ranked against ~20,000 human proteins.** A membrane protein at PSPspredict_PLAAC = 0.3 ranks in the 30th percentile of all human proteins for PLAAC. That is directly interpretable and comparable to other proteins. Re-ranking within 882 changes the meaning.

---

### What to actually do — decision tree

| Analysis goal | Method | Transform? |
|---|---|---|
| ROC curves / AUROC per predictor | Raw scores | None |
| AUPRC per predictor | Raw scores | None |
| Spearman correlation heatmap between predictors | Raw scores | None (Spearman = rank-based) |
| Composite/ensemble score | Percentile rank within 882-protein set | Yes — rank transform, then average |
| Visual distribution comparison | Percentile rank (or ECDF plots) | Yes — for display only |
| Calibrated probability output | Platt scaling with cross-validation | Only if needed; complex |

**Do not:** min-max normalise PSAP or PLAAC; z-score before averaging; apply isotonic regression to n=60; re-rank PSPspredict sub-scores within your subset.

---

### Field consensus — summary

There is no published consensus for cross-predictor LLPS composite scoring. The field avoids the problem by using rank-based metrics only (AUROC/Spearman). The closest to a de facto standard for composite scoring is proteome-wide percentile ranking, as implemented in PSPspredict and implied by PhaSePro's scoring approach. No paper has systematically compared normalisation strategies in this domain.
