# Topology Scoring Pipeline

**Purpose:** Score LLPS predictors on topology-defined sequence regions (Cytoplasmic, Transmembrane, Extracellular/Lumenal) for 60 human membrane proteins, to test whether predicted LLPS propensity differs by topology region and whether the tools agree on region ranking.

**Output:** `output/topology_scores_master.csv` — 7,794 rows, columns: `UniProt_ID`, `tool`, `approach`, `region`, `seg_idx`, `score`.

---

## Scoring approaches

| Approach | Description | Tools |
|---|---|---|
| `region_mean` | Per-residue mean within each topology region, computed from the full-protein per-residue score array | catGRANULE, PLAAC, PScore, ESpritz, SEG, ParSe2 |
| `whole` / `whole_db` | Whole-protein score (fresh local run or precomputed database) | All tools |
| `concatenated` | All residues in a topology class joined into one sequence and scored as a single query | ParSe2, LLPhyScore, PSPsPredict, PSAP, FuzDrop, DeePhase, PSPHunter |
| `segment` | Each contiguous topological span scored individually (one row per segment) | ParSe2, LLPhyScore, PSPsPredict, PSAP, FuzDrop, DeePhase, PSPHunter |

For the per-residue tools (catGRANULE, PLAAC, PScore, ESpritz, SEG), the `concatenated` and `segment` rows are also derived from the full-protein per-residue arrays (slicing, not re-running on isolated sequences). These reflect full-protein context, not isolated-region context.

---

## Tools and coverage

| Tool | whole / whole_db | concatenated | region_mean | segment | Runner |
|---|---|---|---|---|---|
| catGRANULE | 60 (db) | 142* | 142 | 350* | PhaSePred JSON |
| PLAAC | 60 (db) | 159* | 159 | 366* | PhaSePred JSON |
| PScore | 58 (db) | 154* | 154 | 360* | PhaSePred JSON |
| ESpritz | 60 (db) | 159* | 159 | 366* | PhaSePred JSON |
| SEG | 60 (db) | 159* | 159 | 366* | PhaSePred JSON |
| ParSe2 | 60 | 159 | 155 | 366 | `score_topology_parse2.py` |
| LLPhyScore | 60 | 159 | — | 366 | `run_llphyscore_topology.py` |
| PSPsPredict | 60 | 159 | — | 366 | `run_pspspredict_segments.py` |
| PSAP | 60 | 159 | — | 366 | `run_pspspredict_segments.py` |
| FuzDrop | 59 | 107 | — | 101 | `run_fuzdrop_topology.py`, `run_fuzdrop_whole.py` |
| DeePhase | 60 | 154 | — | 343 | `run_deephase_topology.py` |
| PSPHunter | 60 | 159 | — | 366 | `run_psphunter.py` |
| PSPire | 60 (db) | — | — | — | `PSPire_Homo_sapiens_phos_scores.csv` |
| PICNIC | 60 (db) + 58 (fresh) | — | — | — | `picnic auto`; `PICNIC-9606-data.csv` |

\* Derived from full-protein per-residue arrays (same context as `region_mean`), not isolated-sequence re-runs.

Coverage gaps: FuzDrop TM segment/concat is near-zero by design (hydrophobic sequences fail the disorder window). Extracellular concatenated and segment counts are lower than Cytoplasmic/TM because many extracellular regions are very short.

---

## Validation (whole-protein, Spearman ρ vs reference)

| Tool | ρ | n | Reference | Notes |
|---|---|---|---|---|
| PSPHunter | 0.79 | 60 | PSPsPredict master sheet (whole-proteome) | Pure Python reimplementation (Perl script truncates headers) |
| DeePhase | — | 60 | Expected FUS=0.967; got 0.9669 ✓ | No full-proteome reference available |
| PICNIC (fresh vs db) | 0.928 | 58 | PICNIC-9606-data.csv | 2/60 failed: O75445, Q8WXG9 — no AlphaFold model |
| FuzDrop | 0.51 (all) / 0.71 (D_frac>0) | 59/31 | PhaSePred JSON D_frac | 28/59 proteins have ESpritz D_frac=0 in reference; removing those gives ρ=0.71 |
| PSAP | not applicable | 60 | — | MinMaxScaler normalises per-batch; cross-batch Pearson r is not meaningful; rank ordering within approach valid |
| LLPhyScore | ρ≈0.7 | 60 | (from previous session) | |

---

## Key caveats

1. **Per-residue tool concat/segment rows use full-protein context.** The PhaSePred JSON per-residue scores for catGRANULE, PLAAC, PScore, ESpritz, SEG were computed on the whole protein. Slicing those arrays to topology regions preserves full-protein context — this is the correct thing to do for `region_mean` (you want the tool's assessment of residue X in context), but for `concatenated` and `segment` the tool has not seen the isolated region. This is explicitly what we are testing: whether the full-protein score differs between regions.

2. **FuzDrop TM failure is expected.** Hydrophobic TM sequences fail FuzDrop's disorder-window algorithm. Near-zero TM coverage is a property of the tool, not a pipeline bug.

3. **PSAP batch normalisation.** PSAP applies MinMaxScaler across the batch of input sequences, so running 60 proteins gives different raw scores than running 17,800. Rank-based comparisons (Spearman ρ) within each approach are valid; absolute score comparisons are not.

4. **DeePhase ProtVec pickle stub.** The DeePhase word2vec model was pickled with the `ProtVec` class defined in `__main__` scope. The runner script (`run_deephase_topology.py`) defines a minimal `ProtVec` stub class before loading the model so pickle can reconstruct it. The algorithm is unchanged — the stub only provides the attributes gensim needs to unpickle the object.

5. **PICNIC O75445 and Q8WXG9.** These two proteins have no AlphaFold model in the EBI API (HTTP 404). They are absent from the `whole` (fresh) PICNIC scores but present in the `whole_db` scores from the pre-computed PICNIC-9606-data.csv.

---

## TODO

### Structure-based tools: region-level scoring via PDB slicing

PSPire and PICNIC both use AlphaFold2 PDB structures as input. The current master table has only `whole` and `whole_db` scores for these tools — no `concatenated` or `segment` approaches.

To extend these to topology regions:

1. **Download AlphaFold PDBs** for all 60 proteins (EBI API: `https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb`). For proteins >1400 aa, multi-fragment models exist (F1, F2, …).

2. **Slice PDB to topology regions.** Using UniProt topology annotations (already in `data/uniprot_topology_cache.csv`), extract the residue ranges for each region and write region-specific PDB files. BioPython's `PDBIO` with a `Select` subclass can write subsets by residue number.

3. **Run PICNIC `manual` mode** on each region PDB:
   ```
   picnic manual AF-{uid}-{region}-F1.pdb -o output/picnic_region/
   ```
   The filename must follow the `AF-<uid>-F<i>-v<j>.pdb` convention (or `extract_uniprot_id_from_pdb_file` must be patched). Note: scoring a TM-only PDB fragment through PICNIC may not be physically meaningful since STRIDE assigns secondary structure per-residue in isolation, and IUPred2A does not use structure at all.

4. **Run PSPire** on region PDB files. PSPire uses surface-accessible structured patches (RSA > 25%) — this metric is meaningful for soluble regions but undefined for membrane-embedded helices in isolation. Scoring TM fragments is interpretable only as a negative control.

5. **Interpret with caution.** Structure-based tools are validated on soluble proteins; applying them to membrane-extracted fragments is outside their training distribution. Results for Cytoplasmic and Extracellular/Lumenal regions are interpretable; Transmembrane results are not.

### Per-residue tools: isolated-sequence runs (optional)

If true isolated-context scores are needed (rather than full-protein-context slices), standalone runs would require:

- **catGRANULE**: no public binary; would need to use the catRAPID API or reverse-engineer the linear model
- **ESpritz**: no public binary; web server only
- **SEG**: binary available (`segmasker` from BLAST+); could run on individual region sequences
- **PLAAC**: Go binary; `plaac run --input region.fasta` after `go install`
- **PScore**: `pscore` Python package has a broken dependency (`Appium-Python-Client`); not installable in current form
