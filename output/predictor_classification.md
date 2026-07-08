# Predictor Classification for Membrane LLPS Benchmark
Date: 2026-06-17

Classification based on feature type / methodological approach (not algorithm type),
following the framework proposed by J. Hill's collaborator.

---

## Full classification + status

| Tool | Class | Feature basis | LLPS-DB trained? | Leakage risk (our 60) | Proteome scores | In current analysis |
|---|---|---|---|---|---|---|
| **PSAP** | ML | AA composition (RF) | No (literature, pre-DB) | Negligible | XLSX in repo | Yes |
| **DeePhase** | ML | Biophysical + LM embeddings | Yes (LLPSDB) | Low | No (GitHub) | Yes |
| **MolPhase** | ML | 39 physicochemical features (RF) | Yes (multi-DB) | HIGH | No (web submit) | No |
| **PredLLPS_PSSM** | ML | PSSM + Doc2Vec + composition | Yes | HIGH | No (slow) | No |
| **Opt_PredLLPS** | ML | PSSM + HMM + XGBoost | Yes | HIGH | No | No |
| **MambaPhase** | ML | ESM2 + Mamba SSM | Yes | HIGH | No (GPU) | No |
| **PSPredictor** | ML | Word2vec k-mers + GBDT | Yes (LLPSDB v1) | HIGH | No | No |
| **FuzDrop** | Disorder/biophysical | FuzPred + Espritz (theory) | No | None | Partial (PNAS suppl.) | In dataset (p(LLPS)) |
| **ParSe 2.0** | Disorder/biophysical | Polymer scaling + β-turn | No | None | Run-it-yourself | No |
| **ESpritz-DisProt** | Disorder proxy | BRNN disorder predictor | No (disorder, not LLPS) | None | Via PSPspredict | Yes |
| **SEG** | Disorder proxy (LCR) | Sequence entropy | No | None | Via PSPspredict | Yes |
| **PScore** | Physicochemical | π–π contact frequency | No (biophysical scorer) | None | eLife suppl. | Yes |
| **catGRANULE** | Physicochemical | Disorder + RNA-binding + RGG | No (yeast training) | None | Via PSPspredict | Yes |
| **R+Y model** | Physicochemical (baseline) | Arg + Tyr composition | No (rule) | None | Trivially computable | No |
| **LLPhyScore** | Physicochemical/ML | 8 biophysical feature classes | Yes (~565 IDR proteins) | MODERATE | No (GitHub) | No |
| **PLAAC** | Domain/motif | Q/N prion-like HMM | No (yeast PLD HMM) | None | Via PSPspredict | Yes |
| **PDL** | Domain/motif | ProtT5 + KmerConv (prion-like) | Yes (multi-DB) | HIGH | Via PSPspredict | Yes |
| **PSPer** | Domain/motif | PLD + RBD co-occurrence | No (rule-based) | None | Bio2Byte Atlas | No |
| **PICNIC** | Hybrid/whole-protein | RF + AF2 + sequence | Yes (CD-CODE) | VERY HIGH | XLSX in repo | Yes |
| **PICNIC-GO** | Hybrid/whole-protein | RF + AF2 + GO annotations | Yes (CD-CODE) | VERY HIGH | XLSX in repo | Yes (flag for exclusion) |
| **Seq2Phase** | Hybrid/whole-protein (client) | ProtT5 embeddings + stacking | Yes (DrLLPS clients) | HIGH | GitHub | No |
| **PSPHunter** | Hybrid/region-level | Ensemble ML + annotations | Yes (multi-DB) | HIGH | XLSX in repo | Yes |
| **PSPire** | Hybrid/region-level | XGBoost + AF2 structure | Yes (multi-DB) | HIGH | CSV in repo | Yes |
| **dSCOPE** | Hybrid/region-level | RF on SCOPE windows | Yes (experimental SCOPEs) | MODERATE | dscope.omicsbio.info | No |
| **SaPS (PhaSePred)** | Physicochemical/ML | 10-feature gradient boost | Yes (DrLLPS etc.) | MODERATE | predict.phasep.pro | Yes (via PSPspredict) |
| **PdPS (PhaSePred)** | Physicochemical/ML | 10-feature gradient boost | Yes (DrLLPS etc.) | MODERATE | predict.phasep.pro | Yes (via PSPspredict) |
| **PSPspredict (meta)** | Meta-predictor | Aggregates above sub-scores | Inherited | Inherited | XLSX in repo | Yes (sub-scores) |

---

## 2-per-class recommended selection

Rationale: maximise methodological diversity, minimise redundancy, favour tools with
no or low leakage risk where possible, and prefer tools with available proteome scores.

| Class | Recommended pick 1 | Recommended pick 2 | Why |
|---|---|---|---|
| **ML** | PSAP | DeePhase | No leakage (PSAP) + LM embeddings (DeePhase); complementary |
| **Disorder/biophysical** | FuzDrop | ParSe 2.0 | Theory-based (no leakage); FuzDrop whole-protein, ParSe region-level |
| **Physicochemical** | PScore | catGRANULE | π–π vs RNA-binding/RGG; mechanistically complementary |
| **Domain/motif** | PLAAC | PSPer | PLAAC = Q/N prion-like; PSPer = PLD+RBD co-occurrence; different specificity |
| **Hybrid/whole-protein** | PICNIC | Seq2Phase | PICNIC = scaffold bias; Seq2Phase = client-specific (unique in this list) |
| **Hybrid/region-level** | PSPire | PSPHunter | Both in current analysis; complementary AF2 vs annotation features |

**PICNIC-GO**: exclude from primary comparison; show in supplementary as GO-bias reference.
**PSPspredict meta-scores**: exclude from primary comparison; individual sub-scores captured above.
**SaPS/PdPS**: captured under PhaSePred; if included, note they use catGRANULE/PLAAC/PScore as input
  features — not independent.
**MolPhase**: add if scores for 60 proteins can be batch-submitted (adds a 2024 RF model as third ML option).

---

## What needs to be obtained for the recommended selection

| Tool | Action needed |
|---|---|
| FuzDrop | Already in dataset — add to ROC script |
| ParSe 2.0 | Upload membrane protein FASTA to stevewhitten.github.io/Parse_v2_FASTA, download CSV |
| PSPer | Download from bio2byte.be/proteome (all reviewed human proteins) |
| Seq2Phase | Download from github.com/IwasakiLab/Seq2Phase (human proteome predictions) |
| PLAAC, PScore, catGRANULE, PSAP, DeePhase, PICNIC, PSPire, PSPHunter | Already in analysis |
