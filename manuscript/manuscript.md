---
title: "Phase-separation predictors systematically misread membrane proteins: a topology-resolved benchmark of eighteen tools"
author:
  - name: "[Author list to be completed]"
date: "2026"
subtitle: "DRAFT — bullet outline for internal review"
---

> **Status: working draft (bullet form).** Numbers are final and trace to
> `output/results_digest.json`; prose to be expanded at submission. Figures
> embedded inline.

# Abstract (draft points)

- LLPS predictors are trained almost entirely on **soluble** proteins but are
  increasingly applied to **membrane** proteins, whose TM segments and
  constrained cytoplasmic domains violate their assumptions.
- We benchmarked **18 predictors** against **60 human membrane proteins** with
  curated experimental LLPS evidence, scoring against a **same-class
  (membrane-vs-membrane) background** to remove the trivial soluble-vs-membrane
  signal.
- **Structure-aware models led**: PSPHunter AUROC 0.80 (0.90 leakage-corrected);
  PSPire and PdPS ≈ 0.81. **LLPhyScore**, tuned on soluble IDRs, was at
  **chance (0.49)**, even sign-flipped.
- **Topology-resolved scoring**: 7/8 position-aware tools score cytoplasmic
  regions > TM segments (all FDR < 0.05); **SEG** (low-complexity masker) is
  the lone reversal.
- **Cross-predictor consensus tracks the soluble-IDR signature, not curated
  LLPS mechanism**: it is dominated by **hydrophobicity** (ρ = −0.69, strongest
  correlate), rises with **chain length** (ρ = +0.56) and **disorder**
  (ρ = +0.33), and falls with **TMD count** (ρ = −0.41). Large single-pass
  receptors/RTKs with big disordered tails (EGFR, ERBB2, JPH2, CKAP4) rank
  **high**; compact, hydrophobic multi-pass transporters/channels and small
  tail-anchored proteins (MALL, NPY2R, MFSD1) — several with genuine LLPS
  evidence — are systematically **under-ranked**.
- Delineates a **membrane-protein-specific failure mode** — predictors reward
  the generic disorder/length signature and penalise membrane-embedded
  architecture irrespective of experimental mechanism; argues for
  topology-aware training and evaluation.

# Introduction (draft points)

- LLPS organises the cell into membraneless compartments; curated databases
  [@phasepdb2020; @llpsdb2020; @drllps2020; @phasepro2020; @cdcode2023] have
  enabled many sequence- and structure-based predictors.
- Predictor families span: prion-like / low-complexity composition
  [@plaac2014; @seg1993]; π–π and cation–π interaction models [@vernon2018];
  disorder proxies [@espritz2012; @iupred2a2018]; physicochemical scores
  [@catgranule2016]; deep sequence models [@deephase2021]; and ML classifiers
  using AlphaFold features [@alphafold2021; @picnic2024; @pspire2024;
  @psphunter2024; @psap2021; @phasepred2022; @pdl2025].
- **Problem**: almost all trained/validated on soluble proteins. Membrane
  proteins (~25% of the proteome, majority of drug targets) are largely absent
  from training corpora yet routinely scored.
- Membrane proteins break predictor assumptions three ways:
    - TM segments are hydrophobic α-helices whose composition mimics the
      low-complexity regions some tools flag as LLPS-prone;
    - extra-membrane domains are split into cytoplasmic vs extracellular/lumenal
      faces with very different physicochemistry;
    - documented membrane LLPS (TCR/LAT clustering, Wnt signalosome,
      nephrin/NCK/N-WASP) is driven by **cytoplasmic tails**, not the
      membrane-embedded core.
- **This study** — 18 predictors × 60 curated membrane LLPS proteins, with three
  membrane-specific design choices:
    1. same-class (membrane-vs-membrane) background;
    2. protein-by-protein training-set **leakage audit** + leakage-corrected AUROC;
    3. **topology-resolved** scoring (cytoplasmic / TM / extracellular-lumenal).

# Results

## 1. Structure-aware models lead; a soluble-tuned biophysical model fails

![](figures/fig1_benchmark_leakage.png)

- Set: 60 human membrane proteins with curated LLPS evidence (56 PhaSepDB, rest
  LLPSDB); each scored by 18 predictors; discrimination = AUROC vs a background
  of other human membrane proteins.
- **Top tier (structure-aware, AlphaFold features)**: PSPHunter **0.80**, PSPire
  **0.77**, PdPS **0.76** [@psphunter2024; @pspire2024; @phasepred2022].
- **Middle tier**: PICNIC variants + SaPS, AUROC ≈ 0.70–0.71 [@picnic2024].
- **Compositional / disorder proxies** near 0.6: PLAAC 0.62, PSAP 0.59,
  DeePhase 0.58, ESpritz 0.58.
- **At/below chance**: R+Y 0.52, SEG 0.48, **LLPhyScore 0.49**. LLPhyScore is
  biophysically interpretable but tuned on soluble IDRs [@llphyscore2022];
  retains no signal on membrane proteins even sign-flipped. → clearest evidence
  that soluble-IDR inductive biases do not transfer.
- **Training-set leakage (Fig 1, right)** — substantial and *asymmetric*:
    - 25/60 positives (42%) flagged for ≥1 predictor.
    - For the leaders, leaked proteins are mostly **negatives** in the background
      (PSPHunter 22 neg vs 4 pos; PdPS 23 vs 2).
    - Removing leaked proteins therefore **raises** the leaders: PSPHunter →
      **0.90** (95% CI 0.86–0.94), PSPire → 0.81, PdPS → 0.81.
    - PDL (leaked proteins all positives) **drops** 0.63 → 0.58.
    - Leakage-corrected values used throughout.

## 2. Predictor scores follow the membrane topology

![](figures/fig2_topology.png)

- Test: partition per-residue / per-region scores into cytoplasmic, TM, and
  extracellular/lumenal using UniProt topology; compare each region vs
  whole-protein score (paired Wilcoxon, BH-FDR).
- **Structure tools on AlphaFold models**:
    - PICNIC: cytoplasmic **above** whole (median Δ = +0.14, FDR < 10⁻⁷); TM
      **below** (Δ = −0.14, FDR < 10⁻⁷).
    - PSPire: same TM suppression (Δ = −0.11, FDR < 10⁻⁷).
- **Cytoplasmic > TM in 7/8 position-aware tools (all FDR < 0.05, Fig 2b)**:
  PLAAC 98% of proteins, ESpritz 94%, PICNIC 94%, ParSe2 88%, catGRANULE 82%,
  PSPire 80%, PScore 77%.
- Largest median cyto−TM gap: PICNIC Δ = 0.47, PSPire Δ = 0.28.
- **Lone reversal: SEG** (low-complexity masker) — TM scored *higher* (29%
  cyto > TM; median Δ = −0.17, FDR < 10⁻³), as expected for a tool that flags
  hydrophobic α-helices.
- **Read**: given topology, the tools localise their signal to the cytoplasmic
  face — consistent with membrane-condensate biology — even when the
  whole-protein score is uninformative.

## 3. A cross-predictor consensus tracks the soluble-IDR signature, not mechanism

![](figures/fig3_consensus.png)

- Consensus = each tool's per-protein rank (LLPhyScore sign-flipped), averaged
  across the 18 tools; reported as **plain mean rank** (headline) with the
  AUROC-weighted mean as a supplement (the two agree, ρ = 0.996, and give the
  same top/bottom proteins). Between-predictor SD gives a disagreement score.
- **Governed by the classic soluble-IDR signature**, not membrane biology:
    - **Hydrophobicity is the single strongest correlate**: whole-protein
      Kyte–Doolittle hydropathy ρ = **−0.69** — more hydrophobic ⇒ ranked lower.
    - Chain length ρ = **+0.56** — longer ⇒ higher; holds *within* single-pass
      proteins alone.
    - Disordered fraction ρ = **+0.33**; LCR fraction and net charge also
      positive.
    - TMD count ρ = **−0.41** — more membrane passes ⇒ ranked lower.
    - Partial-correlation control confirms **hydrophobicity, length and disorder
      are independent axes** (each survives adjustment for the others), and the
      hydrophobicity effect holds *within* both topology classes — so
      "multi-pass lower" is not merely a length artifact.
- **Consistently high**: large single-pass receptors / RTKs with big disordered
  cytoplasmic tails — JPH2 (0.84), CKAP4 (0.80), EGFR (0.77), ERBB2 (0.75),
  MAVS (0.72), ALK (0.71), LRP6 (0.69), LAT (0.68), ERBB4 (0.65). These are the
  **best-documented** membrane LLPS drivers, and the tools rank them correctly.
- **Consistently low**: compact / hydrophobic proteins — MALL (4 TMD, 153 aa;
  0.19), NPY2R (7-TM GPCR; 0.20), TAZ (0.20), PAR3 (7 TMD; 0.23), MFSD1
  (12 TMD; 0.24), COX7C (63 aa; 0.28), SC6A4 (12 TMD; 0.32). Several
  (MFSD1, PAR3, the multi-pass transporters) carry genuine curated LLPS
  evidence yet fall to the bottom.
- **Why the split** (the "large exposed regions" intuition, confirmed): the
  high-ranked single-pass receptors carry long, disordered cytoplasmic tails
  that match exactly what soluble-IDR predictors learned; the low-ranked set is
  either membrane-dominated (multi-pass, hydrophobic) or too short/compact to
  register. Median length single-pass 770 aa vs multi-pass 503 aa.
- Disagreement (score SD) unrelated to rank (ρ = −0.13, p = 0.32) → both the
  high and low extremes reflect genuine tool **agreement**, not uncertainty; no
  sequence feature explains disagreement (none survive FDR).
- **Post-translational modification burden tracks the same axis**: total
  annotated UniProt PTM sites (phosphorylation, ubiquitination/SUMOylation,
  glycosylation, lipidation) correlate positively with consensus rank
  (ρ = **+0.49**, FDR-p < 0.001), and PTM sites restricted to the cytoplasmic
  region alone give the same signal (ρ = +0.44, FDR-p = 0.001). This tracks
  chain length (longer proteins accumulate more annotated sites) but is not
  fully explained by it: partial correlation controlling for length leaves
  ρ ≈ +0.29–0.33 (still FDR-significant), so PTM density carries information
  beyond a length proxy.
  Mechanistically this is consistent with the same signature as disorder and
  length — heavily modified cytoplasmic tails are also typically long and
  disordered — rather than evidence that any predictor explicitly encodes PTM
  state.

![](figures/fig5_feature_analysis.png)

![](figures/fig2b_predictor_tmd_trend.png)

- The architecture bias is **not** an artifact of the rank-averaging step: it is
  visible in the **raw scores of the individual leading tools**. Across the 60
  proteins, raw score falls monotonically with TMD count for every leader —
  PSPHunter ρ = −0.47, PdPS ρ = −0.44, PSPire ρ = −0.43, FuzDrop ρ = −0.31
  (all p < 0.02) — mirroring the pattern FuzDrop's authors reported for
  multi-pass proteins. Single-pass proteins occupy the high-score band; scores
  collapse as soon as ≥2 TM helices are present.

## 4. Failure modes are architectural, not mechanistic

![](figures/fig4_functional.png)

- Annotated all 60 with UniProt functional class, topology, and curated LLPS
  descriptors (driver/client role; autonomous vs partner-dependent mode);
  tested each vs **plain mean** consensus rank (MWU binary, Spearman continuous,
  Kruskal multi-level; all BH-FDR).
- **Significant — all architectural/functional**:
    - Chain length (Spearman): FDR = 1 × 10⁻⁴ (strongest).
    - Leakage-flagged **higher**: 0.58 vs 0.42 (FDR = 6 × 10⁻⁴).
    - Kinase higher: 0.62 vs 0.46 (FDR = 3 × 10⁻³).
    - RTK higher: 0.62 vs 0.46 (FDR = 4 × 10⁻³).
    - TMD count (Spearman, negative): FDR = 4 × 10⁻³.
    - **Single-pass > multi-pass**: 0.565 vs 0.419 (FDR = 8 × 10⁻³).
    - Transferase higher: 0.61 vs 0.46 (FDR = 0.035).
    - Monotone rise across classes: ion channel 0.36 / transporter 0.39 / GPCR
      0.42 → transferase 0.61 / kinase 0.62 / RTK 0.62.
- **Not significant — the mechanistic descriptors**: driver/client role
  (Kruskal FDR = 0.46); autonomous vs partner-dependent mode (FDR = 0.46).
- **Read**: predictors read architecture (size, membrane passes, hydrophobicity,
  enzyme class), **not** curated LLPS mechanism. They happen to rank the large
  single-pass RTKs correctly — because those look like soluble IDR proteins —
  but for compact or multi-pass proteins that condense via short cytoplasmic
  motifs, the same signature scores them near the floor regardless of the
  experimental evidence.

# Discussion (draft points)

- **Three recurring findings**:
    1. Best generalisers integrate AlphaFold structural features (PSPHunter,
       PSPire, PICNIC); a soluble-IDR biophysical model (LLPhyScore) collapses
       to chance → structural context carries the transferable signal.
    2. Given topology, 7/8 tools correctly localise signal to the cytoplasmic
       face and suppress the membrane core.
    3. Consensus is driven by the generic soluble-IDR signature —
       hydrophobicity (ρ = −0.69), length (+0.56), disorder (+0.33) — and by
       membrane architecture (TMD count −0.41), **not** by curated LLPS
       mechanism (driver/client role and autonomous/partner mode both
       non-significant). The tools rank large single-pass RTKs correctly only
       because those resemble soluble IDR proteins; compact and multi-pass
       proteins with genuine LLPS evidence are penalised for their membrane
       architecture.
- **Two general methodological points**:
    - Leakage direction matters: naïve membrane-background AUROCs *understate*
      the leaders because leaked proteins are mostly easy negatives; a per-tool,
      per-protein audit recovers the honest estimate.
    - Same-class background is essential; scoring vs a soluble proteome inflates
      accuracy by rewarding detection of "membrane-ness."
- **Practical recommendation**: score the **cytoplasmic domain in isolation**,
  not the full-length chain — TM (and, for most tools, extracellular) regions
  add noise, and full-length hydrophobicity dominates the consensus; the
  topology analysis shows the cytoplasmic face is where signal already
  concentrates. This is most consequential for compact and multi-pass proteins,
  whose short disordered motifs are swamped by membrane-embedded sequence.
- **For developers**: topology-aware training/evaluation — topology as an
  explicit feature, membrane-protein positives in training, same-class
  background reporting; do not let whole-chain hydrophobicity stand in as a
  negative LLPS signal for membrane proteins.
- **Limitations**: 60 curated positives; annotation-based topology; consensus is
  in-sample. Coiled-coil and post-translational-modification features
  (phosphorylation of cytoplasmic tails is a known LLPS switch) were not
  computed here and are the clearest follow-up. But the leakage-corrected
  benchmark, topology analysis, and feature-resolved consensus converge on the
  same failure mode: predictors reward the soluble-IDR signature and penalise
  membrane architecture, independent of experimental mechanism.

# Methods

**Protein set.** Sixty human membrane proteins with curated experimental LLPS
evidence were assembled, 56 from PhaSepDB and the remainder from LLPSDB
[@phasepdb2020; @llpsdb2020]. Each was annotated with UniProt membrane topology
(single- vs multi-pass; TM-domain count), chain length, functional class, and
the source database's curated LLPS descriptors (driver/client role; autonomous
vs partner-dependent mode; in-vitro and in-vivo evidence flags)
[@uniprot2023]. The set comprised 39 single-pass and 21 multi-pass proteins.

**Predictors.** Eighteen predictor scores were computed or collated: PICNIC and
PICNIC (GO) [@picnic2024], PSPire [@pspire2024], PSPHunter [@psphunter2024],
PSAP [@psap2021], PhaSePred's SaPS and PdPS sub-scores [@phasepred2022], PDL
[@pdl2025], FuzDrop [@fuzdrop2020], catGRANULE [@catgranule2016], PLAAC
[@plaac2014], PScore and its residue-composition baseline (R+Y) [@vernon2018],
ParSe 2.0 [@parse2_2023], LLPhyScore [@llphyscore2022], PSAP-independent
DeePhase [@deephase2021], the ESpritz disorder proxy [@espritz2012], and the
SEG low-complexity masker [@seg1993]. Structure-based tools (PICNIC, PICNIC
(GO), PSPire) used per-residue AlphaFold2 models [@alphafold2021]; disorder was
computed with IUPred2A where required [@iupred2a2018].

**Benchmarking.** Discrimination was quantified by AUROC against a background of
other human membrane proteins (membrane-vs-membrane), with additional
backgrounds (whole proteome, single-pass-only, multi-pass-only) reported in the
supplement. 95% confidence intervals were obtained by bootstrap. Maximum
Matthews correlation coefficient (MaxMCC) was recorded as a threshold-free
secondary metric. Training-set leakage was audited by matching each of the 60
positives (and each background protein) against every tool's published training
set; per-tool AUROCs were recomputed on the leakage-free subset ("clean").

**Topology-resolved scoring.** For tools producing per-residue or per-region
output, scores were partitioned into cytoplasmic, transmembrane, and
extracellular/lumenal regions using UniProt topology. Each region's score was
compared against the whole-protein score (paired Wilcoxon signed-rank), and
cytoplasmic regions were contrasted directly against TM regions. Because raw
score ranges differ by orders of magnitude across tools, the summary in
**Figure 2b** uses the scale-free fraction of proteins with cytoplasmic > TM
score. All p-values were Benjamini–Hochberg FDR-corrected across tools.

**Consensus.** Each tool's raw score was converted to a normalised per-protein
rank across the 60 proteins (0 = lowest, 1 = highest predicted propensity;
LLPhyScore sign-flipped before ranking so that higher always means more
LLPS-prone). The headline consensus is the **plain mean** normalised rank across
the 18 tools; an AUROC-weighted mean (weights = leakage-corrected
membrane-vs-membrane AUROC) is reported as a supplement and agrees closely
(ρ = 0.996, identical top/bottom proteins). Between-predictor standard deviation
gives a per-protein disagreement score. The full 60 × 18 rank matrix is provided
as Supplementary Table S1 (`output/rank_consensus_matrix.csv`). Associations
between consensus rank and sequence features / annotations were tested with
Spearman correlation (continuous), Mann–Whitney U (binary), and Kruskal–Wallis
(multi-level), all BH-FDR corrected. Sequence features (Kyte–Doolittle
hydropathy, charge descriptors via localCIDER, disordered fraction from ESpritz,
low-complexity fraction from SEG) were computed on whole-protein and
cytoplasmic-only sequences. Post-translational modification burden was pulled
per protein from the UniProtKB REST API [@uniprot2023] (feature types:
Modified residue, Cross-link, Glycosylation, Lipidation), counted both across
the whole sequence
and restricted to residues falling within annotated cytoplasmic topological
domains (same topology boundaries used for the region-based score analysis in
Figure 2).

**Reproducibility.** All scores, statistical tests, and figures are regenerated
by the analysis scripts in the accompanying repository; every quantitative claim
in this manuscript is drawn from a single machine-readable results digest
(`output/results_digest.json`).

# Figure legends

![](figures/fig1_benchmark_leakage.png)

**Figure 1. Structure-aware predictors lead the membrane-vs-membrane benchmark
and training-set leakage is asymmetric.** Left: full-set (grey) and per-tool
leakage-free ("clean", red) AUROC for 18 predictors discriminating 60 membrane
LLPS proteins from a background of other human membrane proteins; error bars are
95% bootstrap CIs, and *n* is the leakage-free positive count. Label colour
encodes predictor mechanism. Right: change in AUROC after leakage removal
(green = rises, red = drops). Removing leaked proteins raises the structure-based
leaders (PSPHunter +0.10, PdPS +0.05, PSPire +0.05) because their leaked
proteins are predominantly background negatives, and lowers PDL (−0.08), whose
leaked proteins are positives.

![](figures/fig2_topology.png)

**Figure 2. Predictor scores follow the membrane topology.** (a) Per-region
minus whole-protein score for the two AlphaFold-based structure tools, PICNIC
and PSPire, across cytoplasmic, transmembrane, and extracellular/lumenal
regions; cytoplasmic regions score above and TM regions below the whole-protein
baseline (paired Wilcoxon, *** FDR < 0.001, n.s. not significant). (b) Fraction
of proteins scoring cytoplasmic > TM for eight position-aware tools; seven of
eight exceed 50% (all FDR < 0.05), and only SEG, a low-complexity masker,
reverses.

![](figures/fig3_consensus.png)

**Figure 3. A cross-predictor consensus for 60 membrane LLPS proteins.** Plain
mean normalised rank (x, higher = consistently predicted phase-separation-prone)
versus between-predictor standard deviation (y, higher = more disagreement).
Point colour/marker encodes single- vs multi-pass topology; shaded bands mark
consistently high (>0.7) and low (<0.3) rank. Large single-pass receptors and
RTKs with long disordered cytoplasmic tails (EGFR, ERBB2, JPH2, CKAP4, MAVS)
populate the high-consensus region; compact and multi-pass transporters,
channels and small tail-anchored proteins (MALL, NPY2R, MFSD1, PAR3, COX7C,
TAZ) populate the low-consensus region. Rank and disagreement are uncorrelated
(ρ = −0.13), so both extremes reflect predictor agreement.

![](figures/fig5_feature_analysis.png)

**Figure 5. The consensus is driven by the generic soluble-IDR sequence
signature.** (a) Spearman correlation of each sequence/architecture/PTM feature
with consensus rank (* FDR < 0.05); whole-protein hydrophobicity is the
strongest (negative) correlate, followed by positive length, total and
cyto-tail PTM-site, LCR and disorder terms, and negative TM-fraction /
TMD-count / net-charge terms. (b) Consensus rank vs Kyte–Doolittle hydropathy
(ρ = −0.69): more hydrophobic proteins rank lower, across both topology
classes. (c) Rank vs disordered fraction (ρ = +0.33), (d) vs log chain length
(ρ = +0.56), and (e) vs total annotated UniProt post-translational
modification sites (ρ = +0.49). Points coloured by topology; extreme proteins
labelled.

![](figures/fig2b_predictor_tmd_trend.png)

**Figure 2b (supporting Figure 3). Leading predictors penalise multi-pass
architecture in their raw scores.** Raw score of each of the four leading tools
(PSPHunter, PSPire, PdPS, FuzDrop) versus transmembrane-domain count for the 60
proteins; black line = binned median. Raw scores fall with increasing TMD count
for all four (Spearman ρ = −0.47, −0.43, −0.44, −0.31; all p < 0.02), confirming
that the architecture bias in the consensus is present in the individual tools
and not introduced by rank-averaging.

![](figures/fig4_functional.png)

**Figure 4. The consensus tracks architecture and function, not curated LLPS
mechanism.** (a) Mean consensus rank by UniProt functional class (median bars);
rank rises from compact multi-pass ion channels, transporters and GPCRs at the
bottom to large single-pass transferases, kinases and receptor tyrosine kinases
at the top. (b) Consensus rank for single- vs multi-pass proteins, coloured by
leakage status; single-pass proteins rank higher (MWU p = 3.2 × 10⁻³), and the
high-ranked single-pass set is enriched for leakage-flagged receptors.

![](figures/figS1_roc_scenarios.png)

**Supplementary Figure S1. ROC curves across background scenarios.** Full ROC
curves for the leading predictors under the four background definitions
(membrane-vs-membrane, vs multi-pass, vs single-pass, vs whole proteome).

![](figures/figS2_predictor_correlation.png)

**Supplementary Figure S2. Cross-predictor score correlation.** Spearman
correlation matrix of the 18 predictor scores across the 60 proteins,
illustrating which tools share underlying signal (e.g. the PhaSePred-derived
sub-scores and the IDR proxies).

![](figures/fig6_kde_score_vs_tmd.png)

**Supplementary Figure S3. The membrane-architecture bias is a proteome-wide
property of the leading predictors, not an artifact of the curated 60-protein
set.** Score density (histogram + KDE, left column), membrane-vs-non-membrane
score distributions (centre column), and score distributions stratified by
TMD count (right column, viridis) for the three leading predictors
(PSPHunter, PSPire, PdPS) and the two worst-performing predictors (SEG,
LLPhyScore) by leakage-corrected membrane-vs-membrane AUROC, computed across
the full human proteome background (`output/background_scored.csv`,
n ≈ 20,400). LLPhyScore's unbounded linear-sum score is shown clipped to its
0.5–99.5th percentile range (see Methods) because its raw range spans
approximately −19,600 to +10,300 and is otherwise dominated by extreme
outliers. Proteome-wide, PSPHunter, PSPire, PdPS and LLPhyScore all show
significant negative score-vs-TMD-count correlation (Spearman ρ = −0.35 to
−0.41, all p < 10⁻⁵ ) and membrane proteins score significantly lower than
non-membrane proteins (Mann–Whitney p < 10⁻⁷ for all four); SEG alone shows
negligible architecture dependence (ρ = +0.04) — consistent with its lone
reversal in the cytoplasmic-vs-TM topology test (Figure 2b) and its
chance-level discrimination in the ROC benchmark (Figure 1).

# References
