# Phase 0b — Predictor Training Set Provenance
Date: 2026-06-17

---

## Critical finding

56/60 of our membrane LLPS benchmark proteins originate from PhasePDB. Five of the
eleven ML predictors (PICNIC, PSPire, PSPHunter, PDL, SaPS, PdPS) trained their
positive sets from PhasePDB or databases that aggregate it (CD-CODE, DrLLPS,
PhaSePro, LLPSDB). Structural leakage is near-certain for these tools. The degree
depends on snapshot dates and protein-level filtering, which requires the actual
training protein lists to resolve.

**Safe for evaluation without downloads:** PSAP, catGRANULE, PLAAC, PScore, ESpritz, SEG
**Require training lists before clean evaluation:** PICNIC, PSPire, PSPHunter, PDL, SaPS, PdPS, DeePhase

---

## Tool-by-tool table

| Tool | Classifier type | N_pos_train | N_neg_train | Training databases | Membrane exclusion | Training list status | Leakage risk (our 60) |
|---|---|---|---|---|---|---|---|
| **PSAP** | ML (Random Forest, AA composition) | 90 | None explicit (proteome BG) | Literature-curated, pre-PhaSepDB | Not mentioned | GitHub `/data/pps_proteins.txt` | **Negligible** — pre-database era, canonical nuclear PSPs |
| **PSPHunter** | ML ensemble (SVM/NB/NN/RF/LGB/XGB) | 135 (human model) | 5,754 | PhaSepDB + LLPSDB + DrLLPS + PhaSePro | Not excluded | Paper Supp Table 1 (needs download) | **HIGH** — 60/60 from same source DBs |
| **PSPire** | ML (XGBoost, AF2 features) | 259 train | 8,323 train | PhaSePred seed + LLPSDB + PhaSePro + PhaSepDB + DrLLPS | Not excluded | Paper Supp Data 4+5 (needs download) | **HIGH** — 60/60 from same source DBs; 64 noID-PSPs from DrLLPS are highest risk |
| **PICNIC** | ML (Random Forest, AF2 + GO) | 2,142 | 1,709 | CD-CODE v1.00 (aggregates PhaSepDB, LLPSDB, DrLLPS, PhaSePro) | Not excluded | Edmond doi:10.17617/3.0Y9Q8N (needs download) | **VERY HIGH** — 58/60 from PhasePDB which CD-CODE directly aggregated |
| **PDL** | ML (ProtT5 + KmerConv) | 640 | 717 | PhaSepDB + DrLLPS + LLPSDB + PhaSePro | Not excluded | GitHub xmuzhanglab/PSPsPredict Datasets S1–S4 | **HIGH** — 60/60 from same source DBs |
| **catGRANULE** | Scoring function (linear) | 120 (yeast granule proteins) | ~4,145 (yeast proteome) | Mitchell et al. 2013 (yeast foci) | N/A — yeast, 2016 | No deposit (cited paper) | **Negligible** — yeast training, human membrane proteins not in set |
| **PLAAC** | HMM scorer (NOT classifier on LLPS labels) | 4→28 yeast prion domains | None | Yeast prion domain sequences | N/A | github.com/whitehead/plaac (HMM params) | **None** — no LLPS labels used |
| **PScore** | Biophysical scorer (NOT classifier) | 11 (pi-contact weight optimisation) | PDB structures | PDB + 11 known LLPS proteins | N/A | eLife source data files | **None** — 11 proteins are canonical IDR PSPs (FUS, hnRNPA1 etc.) |
| **DeePhase** | ML (EF + LM ensemble) | 77 UniProt IDs (137 constructs) | 84 (LLPS−) + 1,563 (PDB*) | LLPSDB (May 2020 snapshot, in vitro homotypic <100 µM) | Not excluded | github.com/kadiliissaar/DeePhase Dataset S1 | **LOW** — in vitro homotypic filter makes membrane protein inclusion unlikely; 2/60 from LLPSDB |
| **SaPS (PhaSePred)** | ML (gradient boosted, 10 features) | 128 | 48,187 (multi-species NoPS) | PhaSepDB + LLPSDB + PhaSePro | Membrane-BOUND organelle proteins excluded from positive set | PNAS Supp Dataset S2 | **MODERATE** — membrane-bound proteins excluded, but "membrane-associated" LLPS drivers may remain |
| **PdPS (PhaSePred)** | ML (gradient boosted, 10 features) | 214 | 48,187 | PhaSepDB + LLPSDB + PhaSePro | Same as SaPS | PNAS Supp Dataset S2 | **MODERATE** |
| **ESpritz-DisProt** | Disorder predictor (BRNN) | Not an LLPS classifier | — | DisProt (disorder, not LLPS) | N/A | Not applicable | **None** — disorder proxy, not LLPS-trained |
| **SEG** | Low-complexity scorer (entropy) | Not a classifier | — | None (algorithmic) | N/A | Not applicable | **None** — algorithmic, no training data |

PSPHunter* and PSPire* are PSPspredict re-runs of PSPHunter and PSPire — same training sets,
same leakage risk, but percentile-normalised to proteome. They inherit the same leakage flags.

---

## Files needed for exact leakage resolution (Phase 0c)

| Priority | Tool | Where to get | File to look for |
|---|---|---|---|
| 1 (highest) | PICNIC | Edmond doi:10.17617/3.0Y9Q8N | Dataset S1 — positive training proteins |
| 2 | PSPire | Nat Comms 15:2147 supplementary | Supplementary Data 4 (train positives) |
| 3 | PSPHunter | Nat Comms 15:2662 supplementary | Supplementary Table 1 (hPS167 list) |
| 4 | PDL | github.com/xmuzhanglab/PSPsPredict | Datasets S1–S4 |
| 5 | PSAP | github.com/vanheeringen-lab/psap | data/pps_proteins.txt |
| 6 | DeePhase | github.com/kadiliissaar/DeePhase | Paper_code/Dataset_S1* |
| 7 | SaPS/PdPS | PNAS 119:e2115369119 supplementary | Dataset S2 (training proteins) |

For each file: extract UniProt IDs and cross-reference against our 60 membrane proteins.
Proteins that match should be excluded from that tool's evaluation (per-tool masking).

---

## Sequence similarity check (Phase 0c extended)

After ID-based cross-reference, run CD-HIT at 30% identity across:
  - All 60 membrane proteins (query)
  - Union of all training sets (reference)

Flag any membrane protein in a cluster with a training protein as a homolog hit.
Install: `sudo apt install cd-hit`

---

## Leakage policy for this study

**Recommended: per-tool masking** (not global clean set)

Global clean set would eliminate all 60 proteins for PSPHunter/PSPire/PICNIC/PDL,
leaving N=0 to evaluate — unusable.

Per-tool masking: each tool evaluated only on proteins it did NOT see in training.
- Requires training lists (see priority table above)
- Tools cannot be directly compared on identical test sets
- Report results with explicit per-tool N and a leakage-flag column
- PSAP, catGRANULE, PLAAC, PScore, ESpritz, SEG can be evaluated on all 60

**Secondary benchmark (cleanest possible):** Restrict to proteins NOT in any of
PhaSepDB / LLPSDB / DrLLPS / PhaSePro at any version. This is a future dataset
(new experimental LLPS membrane proteins), not available now.
