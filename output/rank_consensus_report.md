# Rank-Based Consensus Analysis

## Methods

### Rationale
LLPS predictors output scores on heterogeneous scales (probabilities, log-likelihood
ratios, raw neural-network outputs, etc.) that cannot be meaningfully averaged
directly.  Converting each predictor's scores to normalised ranks removes
scale dependence and makes the consensus robust to outlier scores and
different score distributions.

### Score directionality
All 18 predictors are listed below.  For 17 of them, a higher raw score indicates
greater LLPS propensity.  The exception is **LLPhyScore**, for which a *lower* raw
score is more LLPS-prone (AUROC 0.376 < 0.5 when left as-is).  Its raw score is
sign-flipped before ranking so that higher rank still means more LLPS-prone.

| Predictor | AUROC weight |
| --- | --- |
| PICNIC | 0.670 |
| PICNIC (GO) | 0.698 |
| PSAP | 0.592 |
| PSPHunter | 0.899 |
| PSPire | 0.814 |
| FuzDrop | 0.663 |
| catGRANULE | 0.672 |
| PLAAC | 0.618 |
| PScore | 0.668 |
| ESpritz | 0.577 |
| SEG | 0.481 |
| SaPS | 0.719 |
| PdPS | 0.812 |
| DeepPhase | 0.583 |
| PDL | 0.580 |
| R+Y | 0.517 |
| ParSe2 | 0.680 |
| LLPhyScore (sign-flipped) | 0.488 |

### Normalisation formula
For each predictor, ranks are computed over the N valid (non-NaN) scores using
average-rank ties (scipy.stats.rankdata).  Normalised rank = 1 − (rank − 1) / (N − 1),
mapping the highest-scoring protein to 1.0 and the lowest to 0.0.  Proteins
without a score for a given predictor receive NaN for that predictor's normalised
rank and are excluded from that predictor's ranking.

### AUROC weighting
Weights for the weighted mean rank are taken from the **"Membrane vs membrane
background"** scenario in `output/clean_auroc_table.csv`, using `AUROC_clean`
(leakage-aware) when available, falling back to `AUROC_full` when `AUROC_clean`
is NaN (this affects LLPhyScore only).  Weights are normalised to sum to 1 over
non-NaN predictors per protein before taking the dot product.

### Training set leakage
`output/leakage_map.csv` records whether each protein appeared in the positive
training set of each predictor.  A protein is marked `leakage_any = True` if it
was in the positive training set of *any* predictor.  In the heatmap, individual
cells are asterisked where a protein was in the positive training set specifically
for that predictor.

---

## Key findings

### Top 10 consistently predicted proteins (highest weighted mean rank)

| Entry_name | weighted_mean_rank | rank_sd | leakage_any | TMD_count |
| --- | --- | --- | --- | --- |
| JPH2_HUMAN | 0.840 | 0.228 | True | 1 |
| CKAP4_HUMAN | 0.821 | 0.228 | True | 1 |
| EGFR_HUMAN | 0.788 | 0.194 | True | 1 |
| ERBB2_HUMAN | 0.752 | 0.127 | True | 1 |
| MAVS_HUMAN | 0.742 | 0.277 | True | 1 |
| ALK_HUMAN | 0.726 | 0.199 | True | 1 |
| LAT_HUMAN | 0.699 | 0.272 | True | 1 |
| NOTC1_HUMAN | 0.697 | 0.330 | False | 1 |
| LRP6_HUMAN | 0.685 | 0.183 | False | 1 |
| ERBB4_HUMAN | 0.671 | 0.260 | True | 1 |

### Bottom 10 systematically underpredicted proteins (lowest weighted mean rank)

| Entry_name | weighted_mean_rank | rank_sd | leakage_any | TMD_count |
| --- | --- | --- | --- | --- |
| MALL_HUMAN | 0.171 | 0.214 | False | 4 |
| NPY2R_HUMAN | 0.180 | 0.180 | False | 7 |
| TAZ_HUMAN | 0.208 | 0.138 | False | 1 |
| PAR3_HUMAN | 0.210 | 0.212 | False | 7 |
| MFSD1_HUMAN | 0.219 | 0.244 | False | 12 |
| COX7C_HUMAN | 0.276 | 0.322 | False | 1 |
| CRLS1_HUMAN | 0.282 | 0.274 | False | 5 |
| CD28_HUMAN | 0.289 | 0.221 | False | 1 |
| CAV1_HUMAN | 0.304 | 0.251 | False | 1 |
| SC6A4_HUMAN | 0.308 | 0.182 | False | 12 |

### Most discordant proteins (highest rank SD)

| Entry_name | weighted_mean_rank | rank_sd | leakage_any | TMD_count |
| --- | --- | --- | --- | --- |
| SHSA5_HUMAN | 0.486 | 0.358 | False | 1 |
| NOTC1_HUMAN | 0.697 | 0.330 | False | 1 |
| SYPH_HUMAN | 0.605 | 0.327 | True | 4 |
| COX7C_HUMAN | 0.276 | 0.322 | False | 1 |
| USH2A_HUMAN | 0.578 | 0.315 | False | 1 |
| FZD8_HUMAN | 0.602 | 0.311 | False | 7 |
| VAMP2_HUMAN | 0.444 | 0.302 | False | 1 |
| MUC1_HUMAN | 0.666 | 0.301 | False | 1 |
| AGRV1_HUMAN | 0.561 | 0.296 | False | 7 |
| CD3Z_HUMAN | 0.568 | 0.294 | False | 1 |

### Feature correlations
Spearman ρ was computed between each feature and (a) `weighted_mean_rank` and
(b) `rank_sd`, with Benjamini–Hochberg FDR correction across all 4
(feature × target) combinations.

  - **log₁₀(Length)** is positively correlated with **Weighted mean rank** (ρ = 0.53, BH-adjusted p = 0.000).
  - **TMD count** is negatively correlated with **Weighted mean rank** (ρ = -0.43, BH-adjusted p = 0.001).

See `output/figures/rank_consensus_features.png` for the full barplot with
95 % confidence intervals.

---

## Caveats

- **LLPhyScore** performs poorly on membrane proteins (AUROC 0.376 before sign
  flip, 0.624 after); it is included but interpret its contribution with caution.
- **Training-set leakage** may artificially inflate ranks for flagged proteins;
  use `AUROC_clean` weights (which use leakage-free test sets) to mitigate this.
- **n = 60** limits statistical power for feature correlations; treat non-significant
  results as inconclusive rather than evidence of no effect.
- The consensus reflects **inter-predictor agreement**, not ground truth.  High
  agreement among predictors that share training data may reflect correlated
  biases rather than true LLPS propensity.

---

## Output files

| File | Description |
| --- | --- |
| `output/rank_consensus_table.csv` | Per-protein summary (mean rank, weighted mean rank, SD, leakage flag, metadata) |
| `output/figures/rank_consensus_spearman.png` | Clustered Spearman correlation heatmap between predictor normalised ranks |
| `output/figures/rank_consensus_heatmap.png` | Protein × predictor normalised rank heatmap |
| `output/figures/rank_consensus_scatter.png` | AUROC-weighted mean rank vs rank SD scatter plot |
| `output/figures/rank_consensus_comparison.png` | Side-by-side unweighted vs weighted scatters; direct comparison panel |
| `output/figures/rank_consensus_features.png` | Feature correlations with consensus rank and disagreement |
| `output/rank_consensus_report.md` | This report |
