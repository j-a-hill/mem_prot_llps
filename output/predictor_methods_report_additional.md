# LLPS Predictor Methods Report — Additional Tools
Date: 2026-06-17

Companion to predictor_methods_report.md (which covers the 16 predictors already in the analysis).
This file covers tools proposed by colleagues for potential inclusion.

---

## Summary Table

| Tool | Year | Output type | Scale | Trained on LLPS DBs? | Whole-proteome download | Maintained |
|---|---|---|---|---|---|---|
| MolPhase | 2024 | Probability | 0–1 | Yes (multi-DB) | No (web submit) | Yes |
| PredLLPS_PSSM | 2023 | Probability 3-class | 0–1 | Yes | No (PSI-BLAST, slow) | Yes |
| Opt_PredLLPS | 2024 | Binary + subtype | 0/1 | Yes | No (GitHub only) | No server |
| MambaPhase | 2025 | Probability + multi-class | 0–1 | Yes | No (GPU required) | No server |
| PSPredictor | 2022 | Probability | 0–1 | Yes (LLPSDB v1) | No | Yes (server) |
| FuzDrop | 2022 | p(LLPS) + DPR region | 0–1 | No (biophysical) | Partial (PNAS 2020 suppl.) | Yes |
| ParSe 2.0 | 2023 | Region classification + score | Physicochemical | No | Run-it-yourself (minutes) | Yes |
| PhaSePred | 2022 | SaPS + PdPS scores | Continuous | Yes (DrLLPS, etc.) | **YES** predict.phasep.pro/download | Yes |
| R+Y model | 2018 | Composition ratio | 0–1 | No (rule-based) | Trivially computable | Not a tool |
| PSPer | 2019 | Composite score + regions | Unbounded | No (rule-based) | **YES** bio2byte.be/proteome | Yes |
| LLPhyScore | 2022 | Linear sum | Unbounded | Yes (~565 IDR proteins) | No (GitHub) | No server |
| Seq2Phase | 2023 | Client probability | 0–1 | Yes (DrLLPS clients) | **YES** GitHub IwasakiLab/Seq2Phase | Yes |
| dSCOPE | 2023 | Region binary + score | Per-window | Yes (experimental SCOPEs) | **YES** dscope.omicsbio.info | Yes |

---

## Name clarifications

- **PhaSePred = SaPS/PdPS.** PhaSePred (Chen et al. PNAS 2022) is the publication; SaPS-10fea and PdPS-10fea
  are the two prediction tasks within it. Already in the analysis as PSPspredict_SaPS_10fea and
  PSPspredict_PdPS_10fea (percentile-normalised versions). Whole-proteome raw scores downloadable
  from predict.phasep.pro/download/ — these would be the un-normalised versions.

- **Opt_PredLLPS is a direct successor to PredLLPS_PSSM** (same lab, same framework, larger dataset).
  Only one of the two needs to be considered.

- **R+Y is not a standalone tool.** It is a compositional rule (fraction of Arg + Tyr residues)
  used as a benchmark baseline in the LLPS literature. Computable in one line of Python from any
  FASTA file. Not an independent predictor but a useful feature baseline.

- **FuzDrop is already in the dataset** as the p(LLPS) column in membrane_exp_db_matches.csv
  (60/60 coverage, mean 0.498). It uses biophysical theory (FuzPred + Espritz), not LLPS-DB
  training, making it independent of the database contamination affecting other tools.

---

## Tool-by-tool notes

### MolPhase (EMBO Journal 2024)
Random forest on 39 physicochemical features (IDR, LCR, PLAAC, pi-stacking, charge, composition).
Trained on 606 multi-species LLPS proteins from LLPSDB v2.1, PhaSePro, DrLLPS, PhaSepDB 2.1,
CD-CODE. Negatives: 1,362 PDB globular proteins. AUC 0.972 (test). No bulk download; web server
accepts up to 200 sequences (batch possible for 60 proteins).
Membrane bias: TM helices score low on disorder/composition features. No membrane evaluation.

### PredLLPS_PSSM (Briefings in Bioinformatics 2023)
CNN + BiLSTM on PSSM evolutionary features + Doc2Vec + composition. 3-class output:
LLPS/non-LLPS → PS-Self/PS-Part. Trained on 578 proteins from LLPSDB 2.0 + PhaSepDB 2.1
+ PhaSePro. Slow (PSI-BLAST per protein). Successor = Opt_PredLLPS.

### Opt_PredLLPS (Briefings in Bioinformatics 2024)
Optimised successor to PredLLPS_PSSM. Adds HMM + PSSM features, XGBoost for subtype.
798 training positives from LLPSDB 2.0 + PhaSePro + PhaSepDB 2.1 + DrLLPS. GitHub only,
no server. If only one of {PredLLPS_PSSM, Opt_PredLLPS} is used, prefer this one.

### MambaPhase (Briefings in Bioinformatics 2025)
ESM2 + LoRA fine-tuning + Mamba SSM. Multi-label: LLPS/non-LLPS + scaffold/client +
pH class + salt class. 798 LLPS positives from DrLLPS + LLPSDB 2.0. Requires GPU.
Most architecturally modern tool in this list.

### PSPredictor (BMC Bioinformatics 2022)
Word2vec k-mer embeddings + GBDT. 353 proteins from LLPSDB v1 — smallest and oldest
training set here. Web server at pkumdl.cn/PSPredictor; no bulk download.

### FuzDrop (PNAS 2020 + NAR 2022)
Theory-based (FuzPred + Espritz), NOT trained on LLPS databases. Outputs:
  - p(LLPS) [0,1]: whole-protein droplet probability (threshold 0.60 = driver)
  - DPR: per-residue droplet-promoting region score
**Already in analysis**: p(LLPS) in membrane_exp_db_matches.csv, 60/60 coverage.
Training-set independent — no leakage concern. Most appropriate tool for membrane proteins
among disorder-based predictors because it can output per-region DPR scores separately
from TM helices. Nature Protocols tutorial published 2025.

### ParSe 2.0 (Protein Science 2023)
Rule-based (model polymer scaling exponent + beta-turn propensity). Classifies residues
as folded / IDR / PS-IDR. NOT ML-trained on LLPS labels. Run from FASTA file in minutes,
CSV export, no training data concerns. Best for: identifying which regions of a membrane
protein's cytoplasmic tail are PS-IDR candidates.
No whole-protein score analogous to p(LLPS) — region-level only.

### R+Y model (Vernon et al. eLife 2018)
Composition rule: (count_R + count_Y) / protein_length. No training. Proxy for
cation-π capacity. TM helices dilute R/Y content → systematic underestimation for
membrane proteins. Computable trivially. Use as a feature baseline, not primary predictor.

### PSPer (Bioinformatics 2019)
Detects prion-like domain (PLD) + RNA-binding domain (RBD) co-occurrence. Rule-based
feature combination. Specific to the FUS/TDP-43/hnRNP class. Will miss most membrane
LLPS mechanisms. **Whole-proteome scores available from Bio2Byte Proteome Atlas**
(bio2byte.be/proteome), updated August 2023. Narrow specificity but directly downloadable.

### LLPhyScore (Biomolecules 2022)
Interpretable linear model on 8 biophysical feature classes (PDB-derived statistics).
~565 curated IDR proteins for positive training. Linear decomposition allows feature-level
interpretation. "Disordered protein phase separation" in title — explicitly IDR-centric.
GitHub only; no server. Worth considering for interpretability but limited membrane applicability.

### Seq2Phase (Bioinformatics Advances 2023)
ProtT5-XL-U50 embeddings + stacking ML ensemble. Unique: the ONLY tool here specifically
designed to predict **client** proteins (recruited to condensates, not scaffold/drivers).
DrLLPS client annotations for training. Human proteome predictions (20,398 proteins) available
on GitHub (IwasakiLab/Seq2Phase). Language model embeddings include TM context; the classifier
above may partially generalise to membrane client proteins.

### dSCOPE (Briefings in Bioinformatics 2023)
Region-level Random Forest using experimental SCOPE (Sequence Critical for Phase sEparation)
windows from mutagenesis literature. Training data = experimentally validated critical segments,
not whole-protein labels. Per-window output. Precomputed for all reviewed human proteins at
dscope.omicsbio.info — searchable by UniProt ID. The only region-level predictor with a
precomputed human proteome searchable database.

---

## Immediately addable to the analysis (scores obtainable now)

| Tool | What to do | Class |
|---|---|---|
| FuzDrop | Already in dataset — add to ROC | Disorder/biophysical |
| PSPer | Download from bio2byte.be/proteome | Domain/motif |
| Seq2Phase | Download from GitHub IwasakiLab/Seq2Phase | Hybrid/ML-client |
| dSCOPE | Query 60 proteins at dscope.omicsbio.info | Region-level |
| R+Y | Compute from UniProt FASTA (one script) | Physicochemical baseline |
| MolPhase | Batch-submit 60 proteins to molphase.sbs.ntu.edu.sg | ML |

## Requires running locally (non-trivial)

| Tool | Blocker |
|---|---|
| ParSe 2.0 | Need UniProt FASTA file; run takes minutes |
| LLPhyScore | GitHub clone + Python environment |
| Opt_PredLLPS | GitHub + HMM/PSSM setup |
| MambaPhase | Requires GPU |
| PredLLPS_PSSM | PSI-BLAST, slow (~hours for 60 proteins) |
| PSPredictor | Web server, no bulk; tedious for 60 proteins |
