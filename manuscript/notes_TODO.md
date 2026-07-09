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

## Reminder for next session

- Git remote (SSH) push still not confirmed working from this sandbox —
  push manuscript/repo changes from your own machine.
- If coiled-coil is later added: re-run the Fig 5 rebuild cell (uses the
  `battery3`/`scatter3` helpers already defined this session; extend
  `feats2` and `LAB` the same way PTM was added), then do ONE combined
  DOCX re-render rather than one per feature.
