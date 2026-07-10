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
- **Fig 9 + sequence-region completion (superseded/extended, see next
  bullet)**: built unified wide table
  (`output/predictor_region_score_comparison_wide.csv`) with whole /
  concatenated (or `pdb_region` for PICNIC & PSPire) / segment-mean scores
  side by side for all 14 region-capable tools, from
  `output/topology_scores_master.csv`.

## Settled this session (2026-07-09/10, structural concatenation + 2-factor Fig 9 redesign)

- **Structural `pdb_concatenated` built for PICNIC and PSPire**, the missing
  piece needed to give the two structure-based tools a true analogue of
  sequence concatenation (not just `pdb_region` per-fragment scoring).
  `build_concatenated_pdbs.py` renumbers each region's (possibly disjoint)
  true AlphaFold-model fragments into one continuous chain (1..N),
  preserving true 3D coordinates — this was the user-specified "Option 1"
  structural analogue, chosen over generating a fresh AlphaFold/ESMFold
  structure from the concatenated FASTA (deferred, see below). Scored via
  `run_picnic_pdb_concatenated.py` (146/146 rows, PICNIC manual-scoring API)
  and `run_pspire_pdb_concatenated.py` (146/146 rows, PSPire.py CLI);
  merged into `topology_scores_master.csv` via `integrate_pdb_concatenated.py`.
- **Fig 9 redesigned as a 2-panel, 2-factor figure**
  (`fig9_topology_vs_splicing.png`) that cleanly separates (a) the
  topology-class biology effect, tested on unspliced (segment-mean) scores
  only (`output/predictor_region_biology_tests.csv`), from (b) the
  slicing/splicing methodological artefact, tested as spliced-vs-unspliced
  within each region (`output/predictor_splice_artifact_tests.csv`) — both
  now computed uniformly across all 14 region-capable tools, sequence- and
  structure-based alike, replacing the old single-multiples/whole-vs-region
  design that conflated the two questions.
- **Key new finding**: PICNIC and PSPire both show a strong, tool-consistent
  splicing artefact specific to the Transmembrane region — concatenating
  disjoint single-pass TM fragments into one continuous structure lowers the
  score in nearly every protein (PICNIC 0/56, PSPire 1/56 scored higher when
  concatenated; both FDR < 10⁻³) — with no comparable effect in the
  Cytoplasmic or Extracellular/Lumenal splicing tests. The topology-class
  effect (Cyto > TM) itself is unaffected by splicing (it is tested purely
  on unspliced scores) and holds for both structure tools.
- **Resolved the FuzDrop/ParSe2 coverage-gap flag** raised mid-session:
  confirmed (not a pipeline bug) — both tools have a real minimum
  sequence-length requirement (ParSe2: hard `MIN_LEN=25`; FuzDrop has a
  similar practical disorder-window floor), and single-pass Transmembrane
  segments are short (median 21 aa), so segment-level scoring
  systematically under-covers TM for these two tools specifically. Now
  stated explicitly in the Methods paragraph.
- Manuscript.md updated: Results section 2's region-scoring bullet rewritten
  around the two-factor design with concrete per-tool numbers; Methods
  region-comparison paragraph rewritten to describe both the sequence
  `concatenated` and structural `pdb_concatenated`/`pdb_region` approaches
  and the two statistical tests; Fig 9 caption rewritten to match the new
  panels; Discussion bullet softened/corrected — it previously claimed the
  region-level pipeline was "now complete," which was premature (the
  structural splicing test didn't exist yet); it now describes what was
  actually established (topology effect vs splicing artefact, separated,
  for all 14 tools). **DOCX re-render is PENDING** (batched per the workflow
  preference below) — do it in the same pass as any other outstanding
  manuscript.md edits before next artifact save.

## Deferred: fresh-structure concatenation check (Option 2)

- User's original request had two complementary structural-concatenation
  approaches: **Option 1** (renumber true AlphaFold fragments into one
  chain, keep real coordinates — chosen as primary, done this session) and
  **Option 2** (generate a *fresh* AlphaFold2 or ESMFold structure directly
  from the concatenated-region FASTA, so the fold is re-predicted for the
  spliced sequence rather than reusing per-fragment coordinates). Option 2
  is NOT started — would need to check whether an `esmfold`-type skill/tool
  is available and compute-feasible (ESMFold is lighter-weight than AF2 and
  more plausible to run in this environment); if feasible, it would give an
  orthogonal check on whether the PICNIC/PSPire TM-splicing effect found
  above is a renumbering/coordinate-discontinuity artefact specifically, or
  survives when the model is allowed to re-fold the spliced sequence.

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
