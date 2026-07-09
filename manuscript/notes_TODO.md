# Manuscript notes / open items

Working notes for future revisions. Not part of the manuscript itself —
kept separate so we don't have to re-render `manuscript.docx` for every
small addition.

## Deferred analysis

- **Coiled-coil content** — not yet computed. Original feature-analysis
  spec (alongside disorder / LCR / topology / cyto-tail / charge / PTM)
  wanted this checked as a candidate correlate of consensus rank. No
  dedicated tool (COILS, Marcoil, ncoils) is installed in `memllps`; would
  need either (a) installing one of those, or (b) a heuristic heptad-repeat
  / hydrophobic-moment window scan as a stand-in. Deferred by user decision
  on 2026-07-09 — revisit if a future draft wants a more complete feature
  battery, or if reviewers ask about coiled-coil-mediated self-association
  as an alternative LLPS-driving mechanism for membrane proteins (relevant
  e.g. for SNARE-family and some receptor tail proteins in the 60-set).

## Settled this session (for reference, not re-litigated)

- PTM burden added as a feature: total UniProt PTM sites ρ = +0.49
  (FDR p < 0.001), cyto-tail-restricted PTM ρ = +0.44 (FDR p = 0.001) vs
  weighted_mean_rank. Partial correlation controlling for Length:
  ρ ≈ +0.29–0.33, still FDR-significant — PTM density is not fully a length
  proxy. Folded into Fig 5 panel (e) and Results bullet list.
- Fig 6 (proteome-wide KDE, PSPHunter/PSPire/PdPS vs SEG/LLPhyScore across
  ~20,400 background proteins) added as **Supplementary Figure S3** —
  confirms the membrane-architecture score bias (ρ ≈ −0.35 to −0.41 vs
  TMD count) is a general property of the leading predictors, not an
  artifact of the curated 60-protein set. SEG is the outlier (ρ = +0.04),
  consistent with its lone topology reversal (Fig 2b) and chance-level
  ROC performance (Fig 1). LLPhyScore panels are percentile-clipped
  (0.5–99.5th pctile) because its raw score is an unbounded linear sum
  (range ≈ −19,600 to +10,300) — noted explicitly in the S3 caption.
- All reference DOIs CrossRef-verified, including the two that were
  previously missing: ParSe2 = 10.1002/pro.4756; LLPhyScore =
  10.3390/biom12081131.
- manuscript.md / manuscript.docx / references.bib / references_QC.md all
  re-saved as new artifact versions after the above changes. DOCX QC passed:
  9 embedded images, 0 unresolved `[@key]` citation markers, References
  section present.

## Settled this session (2026-07-09, consensus/discordance + region-scoring completion)

- **Fig 7 (consensus/discordance)**: 60x18 rank heatmap + mean_rank-vs-rank_sd
  scatter + strip plots for consistently-high (n=8, all single-pass) vs
  consistently-under-ranked (n=9, 6/9 multi-pass) + top-10 disagreement
  proteins. Groups saved to `output/consensus_discordance_groups.csv`.
  Topology-class split between high/under groups is itself significant
  (Fisher's exact p=9.1e-3). Disagreement (rank_sd) uncorrelated with
  mean_rank (ρ=0.14, p=0.30) and with every tested feature after FDR
  (all FDR>0.2) — genuine negative result, written up as such.
- **Fig 8 (group feature comparison)**: boxplots of 8 FDR-significant
  features (hydropathy, length, disorder_frac, LCR_frac, TMD_count,
  total_PTM, whole_FCR, cyto_frac) separating high vs under-ranked groups
  (Mann-Whitney U, BH-FDR). Full test table in
  `output/group_feature_comparison_tests.csv`.
- **Fig 9 + sequence-region completion**: built unified wide table
  (`output/predictor_region_score_comparison_wide.csv`) with whole /
  concatenated (or `pdb_region` for PICNIC & PSPire) / segment-mean scores
  side by side for all 14 region-capable tools, from
  `output/topology_scores_master.csv`. Coverage gap (PICNIC/PSPire lack
  concatenated+segment sequence-based scoring, structure-based `pdb_region`
  used instead) now explicitly tabulated in
  `output/predictor_coverage_summary.csv` and stated in Methods/Fig 9
  caption rather than left implicit. Small-multiples figure reproduces the
  same TM-suppression/SEG-reversal pattern as the coarser topology test.
- Manuscript.md updated with new Results section 4 (consensus/discordance),
  extended section 2 topology bullet (region-scoring completion), new
  Methods paragraphs (classification thresholds, region-scoring pipeline),
  Discussion bullets, and Fig 7/8/9 legends. **DOCX re-render is PENDING**
  (batched per the workflow preference below) — do it in the same pass as
  any other outstanding manuscript.md edits before next artifact save.

## Reminder for next session

- Git push: **use HTTPS + `GITHUB_TOKEN` env var**, confirmed working this
  session (SSH/port-22 is blocked by sandbox network policy — do not
  retry SSH). Pattern:
  `git -c http.extraHeader="Authorization: Basic $(printf 'x-access-token:%s' "$GITHUB_TOKEN" | base64 -w0)" push https://github.com/j-a-hill/mem_prot_llps.git minimal:minimal`
- If coiled-coil is later added: re-run the Fig 5 rebuild cell (uses the
  `battery3`/`scatter3` helpers already defined this session; extend
  `feats2` and `LAB` the same way PTM was added), then do ONE combined
  DOCX re-render rather than one per feature.
- Batch DOCX re-renders — don't re-render after every individual
  manuscript.md edit within a session.
