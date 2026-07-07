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
| PSPire | 60 (db) | — (pdb_region: 146) | — | — | `PSPire_Homo_sapiens_phos_scores.csv`; `run_pspire_pdb_region.py` |
| PICNIC | 60 (db) + 58 (fresh) | — (pdb_region: 146) | — | — | `picnic auto`; `PICNIC-9606-data.csv`; `run_picnic_pdb_region.py` |

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

## `pdb_region` approach (PICNIC — complete; PSPire — complete)

PICNIC and PSPire both use AlphaFold2 PDB structures as input. The `pdb_region` approach is the structural analogue of `concatenated` for sequence-based tools: each topology region is sliced from the full AF2 PDB and scored as if it were a standalone protein.

### What `pdb_region` does and does not represent

| Feature | Source | Effect of region slice |
|---|---|---|
| pLDDT | B-factor column (AF2 full-protein prediction) | ✓ Preserved exactly — full-protein AF2 confidence values |
| STRIDE secondary structure / RSA | 3D coordinates | ✗ Changes at TM-boundary residues; interior of large domains unaffected |
| IUPred2A disorder | Sequence extracted from ATOM records | ✗ Boundary effects at TM junctions; minor for large regions |
| Sequence composition / complexity | ATOM-derived sequence | ✓ Correct for the region |

SEQRES header records are deliberately stripped from region PDB files so PICNIC falls back to ATOM-based sequence extraction, ensuring IUPred sees only the region sequence (consistent with the `concatenated` approach for sequence tools).

**TM region results are a negative control** — isolated hydrophobic helices appear maximally surface-exposed; STRIDE secondary structure assignments change at cut-points. Analogous to FuzDrop's TM failure.

### PICNIC `pdb_region` — ✅ complete

**Script:** `run_picnic_pdb_region.py`  
**Integration:** `integrate_picnic_pdb_region.py`  
**Output:** `output/topology_scores_raw/PICNIC_pdb_region_raw.csv` (146 rows, 58 proteins)

Coverage: 58/60 proteins (O75445 and Q8WXG9 have no AF2 model). Regions skipped if annotated residue count < 10.

| Region | n | Median | Mean |
|---|---|---|---|
| Cytoplasmic | 50 | 0.63 | 0.57 |
| Extracellular/Lumenal | 40 | 0.37 | 0.43 |
| Transmembrane | 56 | 0.14 | 0.16 |

41/50 proteins score Cytoplasmic > whole-protein (mean Δ = +0.20), consistent with cytoplasmic IDRs driving LLPS propensity. TM scores are low throughout, as expected.

**Run command:**
```bash
PYTHONPATH=/home/jake/iupred2a_lib \
PATH="/home/jake/miniconda3/envs/PSPHunter/bin:$PATH" \
/home/jake/miniconda3/envs/PSPHunter/bin/python3 run_picnic_pdb_region.py
```

AlphaFold2 PDB files (v6, downloaded from EBI API) are cached to `data/alphafold_pdbs/` (gitignored — large binary files, reproducible by re-running the script). Region slices are saved to `data/alphafold_pdbs/regions/{region}/`.

### PSPire `pdb_region` — ✅ complete

PSPire (Hou et al., *Nature Communications* 15, 2147 (2024); DOI: 10.1038/s41467-024-46445-y) predicts LLPS via surface-exposed structured patches (SSUPs: non-IDR residues with RSA > 25%) computed from AlphaFold2 structures. It uses DSSP for RSA and an XGBoost model.

**Script:** `run_pspire_pdb_region.py`  
**Integration:** `integrate_pspire_pdb_region.py`  
**Output:** `output/topology_scores_raw/PSPire_pdb_region_raw.csv` (146 rows, 58 proteins)

Coverage: 58/60 proteins (O75445 and Q8WXG9 have no AF2 model). PSPire code cloned from `github.com/TongjiZhanglab/PSPire` into `external_tools/PSPire/` (gitignored).

| Region | n | Median | Mean |
|---|---|---|---|
| Cytoplasmic | 50 | 0.305 | 0.366 |
| Extracellular/Lumenal | 40 | 0.117 | 0.254 |
| Transmembrane | 56 | 0.041 | 0.059 |

Cytoplasmic > Extracellular/Lumenal > Transmembrane, mirroring the PICNIC pattern. PSPire SSUP scores are generally higher in cytoplasmic regions, consistent with disordered-but-structured cytoplasmic regions forming more surface patches.

**Run command:**
```bash
python3 run_pspire_pdb_region.py
python3 integrate_pspire_pdb_region.py
```

**Technical notes:**

- PSPire was cloned from GitHub (not available as a pip package). The `lib/dssp_run.py` in the clone has two compatibility patches applied for the system mkdssp v4.2.2 (vs the expected v2.x):
  1. Command syntax: `mkdssp --output-format dssp input output` (v4) instead of `mkdssp -i input -o output` (v2)
  2. Mapping dict: added `'P': 0` to handle the polyproline II helix class (`P`) new in DSSP v4; without this, unmapped residues cause float promotion and `ValueError: invalid literal for int() with base 10: '.'`

- Region PDB files (from the PICNIC pipeline) have HEADER/TITLE/COMPND/SOURCE lines prepended from the full AF2 PDB before passing to PSPire; DSSP v4 requires a HEADER record to recognise files as PDB format.

- PSPire uses pLDDT B-factors to identify IDRs (threshold < 50). For region slices, these B-factors are from the full-protein AF2 prediction (correctly preserved), so IDR identification reflects the original structural confidence.

**Interpretation caveat:** PSPire's SSUP metric (RSA > 25%) requires a meaningful solvent-accessible surface. For TM region slices, residues that are normally buried in the lipid bilayer appear maximally surface-exposed in isolation. TM scores are uninformative and serve as a negative control (same caveat as PICNIC TM results).

## TODO (remaining)

The `pdb_region` approach for both structure-based tools (PICNIC, PSPire) is now complete. The master table (`output/topology_scores_master.csv`) has 8,086 rows covering all 14 tools and all implemented approaches.

### Per-residue tools: isolated-sequence runs (optional)

If true isolated-context scores are needed (rather than full-protein-context slices), standalone runs would require:

- **catGRANULE**: no public binary; would need to use the catRAPID API or reverse-engineer the linear model
- **ESpritz**: no public binary; web server only
- **SEG**: binary available (`segmasker` from BLAST+); could run on individual region sequences
- **PLAAC**: Go binary; `plaac run --input region.fasta` after `go install`
- **PScore**: `pscore` Python package has a broken dependency (`Appium-Python-Client`); not installable in current form
