# Predictor Categories — Short Report

18 predictors are scored against the membrane LLPS benchmark. Each is classified
on two independent axes used to group/colour the plots in `plot_roc.py`.

## Axis 1 — Target type

What the tool was actually designed or trained to predict.

| Target type | Predictors | Meaning |
|---|---|---|
| **LLPS-specific** | PICNIC, PICNIC (GO), PSPire, PSPHunter, SaPS, PdPS, PDL, LLPhyScore, PScore, R+Y, ParSe2, FuzDrop, PSAP, DeepPhase, catGRANULE | Trained on, or theoretically derived for, phase separation specifically. |
| **IDR proxy** | PLAAC, ESpritz, SEG | General-purpose disorder/low-complexity/prion-domain tools, never trained on LLPS labels. Used here only as weak proxies — expect lower AUROC, and read their results with that caveat. |

On the grids, the IDR-proxy group is forced onto its own row, and its titles
are italicised with a "(IDR proxy)" suffix.

## Axis 2 — Mechanism

The dominant feature or theory driving the method, per
`output/predictor_classification.md`.

| Mechanism | Colour | Predictors |
|---|---|---|
| **AF2 structure + sequence** | blue `#3366CC` | PICNIC, PICNIC (GO), PSPire — the only three that actually take AlphaFold2-predicted structure as input. |
| **Multi-feature ML ensemble** | purple `#7755AA` | PSPHunter (ensemble ML + annotations, not structure), SaPS, PdPS (PhaSePred 10-feature gradient boost), PDL (ProtT5 + KmerConv), LLPhyScore (8 biophysical feature sum). None of these use a structural model. |
| **Pi-pi / cation-pi** | red `#CC3333` | PScore (π–π contact frequency), R+Y (Arg+Tyr composition). |
| **Sticker-spacer / polymer theory** | amber `#CC8800` | ParSe2 (polymer scaling exponent vs. helix propensity), FuzDrop (context-dependent stickers within IDRs). |
| **IDR / disorder** | green `#228855` | PSAP, DeepPhase, catGRANULE (LLPS-trained), plus PLAAC, ESpritz, SEG (general-purpose, IDR-proxy). |

## Caveats baked into every figure's caption

- **SaPS and PdPS are PhaSePred's own scores** (Chen et al., PNAS 2022;
  SaPS-10fea/PdPS-10fea are its two prediction tasks — scaffold-forming vs.
  partner-dependent/client), sourced here from PhaSePred's raw JSON. Both
  take catGRANULE/PLAAC/PScore/ESpritz/SEG/DeepPhase scores as *input
  features* — they are **not independent** of those six tools also in this
  comparison; treat their agreement with each other, or with that group, as
  partly circular rather than corroborating evidence.
- **PDL is not from PhaSePred.** It's Wenbin Li's own tool ("Protein
  Dual-model Language", LLM + KmerConv, published as PSPsPredict on GitHub).
  His benchmark compilation spreadsheet (`PSPspredict_full_proteome.xlsx`)
  re-publishes PhaSePred's SaPS/PdPS (percentile-ranked) alongside his own
  PDL score for comparison in his own paper — but our `SaPS_score`/`PdPS_score`
  come from PhaSePred's raw JSON, not that spreadsheet, so PDL has no actual
  numeric link to SaPS/PdPS here, despite the superficial filename overlap.
- **LLPhyScore's true-positive set is keyed by construct/gene name, not
  UniProt accession** (e.g. "GFP_CSTF2", "APVGVG_35", "NICD_WT_GFP"), and
  it's a heterogeneous mix: real human proteins, domain fragments,
  GFP-fusion constructs, synthetic elastin-like repeat peptides, and
  non-human orthologs (Arabidopsis, C. elegans, etc.). Exact gene-symbol
  matching against our 60 benchmark proteins found zero hits, but **exact/
  fragment sequence matching** against the human proteome (either direction
  — training sequence found inside a protein, or a protein found embedded
  inside a longer fusion construct) resolved **66 of the 305 training
  sequences (22%) to real UniProt IDs**, with zero ambiguous (>1 ID)
  matches. Two of our 60 membrane benchmark proteins are confirmed
  **positive** training leaks this way: **O60500/Nephrin/NPHS1** (matched
  via "NICD_WT_GFP"/"NICD_delNTD_GFP" — Nephrin's intracellular domain is a
  classic LLPS clustering study system) and **P08908/HTR1A** (matched
  exactly via "5HT1A"). The remaining 239 unresolved sequences
  (fusions/non-human orthologs/synthetic repeats with no real human-protein
  counterpart) stay genuinely `unknown` for every other benchmark protein —
  not falsely "confirmed clean".
- LLPhyScore's **negative** training set
  (`LLPhyScore_human_negatives_training_tags.pkl`, fetched from
  `julie-forman-kay-lab/LLPhyScore`'s `data/processed/training/`) is a
  2000-protein generic-human-proteome sample keyed
  `ACCESSION_GENE_ORGANISM` — directly UniProt-mappable, no sequence
  matching needed. **7 of our 60 membrane proteins are confirmed LLPhyScore
  training negatives**: LAT (O43561), ERN1/IRE1 (O75460), APP (P05067),
  NOTCH1 (P46531), CKAP4 (Q07065), FZD8 (Q9H461), F11R/JAM1 (Q9Y624) — their
  LLPhyScore scores are likely deflated.
- Net effect: `clean_LLPhyScore_score` is `False` for all 60 (2 confirmed
  positives, 7 confirmed negatives, 51 still unknown — none confirmed
  *both* non-positive and non-negative), and the global-clean set excludes
  exactly these 9 proteins as confirmed contaminated rather than unverified
  unknowns.
