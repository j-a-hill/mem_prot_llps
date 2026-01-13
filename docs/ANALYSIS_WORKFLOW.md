# LLPS Protein Analysis Workflow

## Overview

This analysis pipeline identifies functional groups of membrane proteins enriched in phase-separation propensity (pLLPS), fetches their protein-protein interactions from STRING, and visualizes the resulting networks with nodes colored by pLLPS scores.

## Analysis Notebooks

Run notebooks in sequence for comprehensive analysis:

### 01: Data Loading and Classification
**File:** `01_data_loading_and_classification.ipynb`

**Purpose:** Load and classify the LLPS dataset

**What it does:**
- Loads the Human Phase Separation dataset (20,366 proteins)
- Classifies proteins into High/Medium/Low pLLPS categories
- Identifies 6,463 membrane proteins
- Generates initial visualizations

**Outputs:**
```
results/full_dataset.csv              # Complete classified dataset
results/high_pllps_proteins.csv       # High pLLPS proteins (≥0.7)
results/membrane_proteins.csv         # Membrane proteins only
results/classification_summary.json   # Summary statistics
```

**Runtime:** ~30 seconds

---

### 02: pLLPS-Enriched Functional Groups
**File:** `02_pllps_enriched_functional_groups.ipynb`

**Purpose:** Identify functional groups enriched in pLLPS scores

**What it does:**
1. Loads pLLPS scores and membrane protein data
2. Classifies membrane proteins by functional group:
   - Structural proteins
   - Ion Channels
   - GPCRs/Receptors  
   - Receptor Tyrosine Kinases
   - Transporters
   - Enzymes
   - Other membrane proteins
3. Analyzes pLLPS distribution within each group
4. Performs chi-squared tests to identify significantly enriched groups (p < 0.05)
5. Calculates enrichment factors
6. Visualizes enrichment patterns

**Key Finding:** Structural proteins show ~1.21x enrichment (significantly enriched in high pLLPS)

**Outputs:**
- `functional_groups_with_pllps.csv` - Membrane proteins with functional classification
- `functional_enrichment_stats.csv` - Statistical analysis per group
- `pllps_enriched_groups.json` - List of significantly enriched groups
- `functional_groups_pllps_enrichment.png` - Enrichment visualization

**Runtime:** ~5 minutes

---

### 03: STRING Network Analysis
**File:** `03_string_networks_pllps_enriched.ipynb`

**Purpose:** Fetch STRING interactions for enriched functional groups

**What it does:**
- Loads enriched groups from notebook 02
- Fetches high-confidence STRING interactions (score ≥ 700/1000)
- Matches interactors with pLLPS scores
- Calculates interaction statistics
- Processes groups: Structural, Ion Channel, Transporter, GPCR

**Outputs:**
- `results/string_networks_by_group/{group}_interactions.csv`
- `results/enriched_group_interactions_summary.json`

**Runtime:** ~10-15 minutes (depends on group sizes)

---

### 04: Network Visualization
**File:** `04_visualize_pllps_networks.ipynb`

**Purpose:** Visualize pLLPS-colored interaction networks

**What it does:**
- Generates network graphs for each functional group
- Colors nodes by pLLPS score (blue=low, white=medium, red=high)
- Sizes nodes by connectivity degree
- Weights edges by STRING confidence score

**Outputs:**
- Publication-quality PNG networks (300 DPI)
- Interactive HTML visualizations (PyVis)
- Network statistics (CSV + JSON)

**Runtime:** ~2-5 minutes

---

### 05: Interactive Functional Group Networks
**File:** `05_interactive_functional_group_networks.ipynb`

**Purpose:** Create detailed interactive network visualizations

**Features:**
- ✅ Node Coloring by pLLPS Score (blue=low, white=medium, red=high)
- ✅ Edge Weighting by STRING Score
- ✅ Node Size by Connectivity
- ✅ Multiple Output Formats (PNG, HTML, CSV)

**How to Use:**
1. Run cells 1-3 to load data
2. In Cell 5, change `TARGET_GROUP` to desired functional group
3. Run remaining cells to generate networks

**Runtime:** ~1-2 minutes per group

---

## Quick Start

### Option 1: Full Analysis (Recommended)
```bash
# Run all notebooks in order
jupyter notebook
# Open and run: 01 → 02 → 03 → 04 → 05
```

### Option 2: Quick Exploration
```bash
# Use the interactive dashboard
streamlit run shiny_app.py
```

## Directory Structure

```
/workspaces/mem_prot_llps/
├── docs/                              # Documentation
│   ├── ANALYSIS_WORKFLOW.md          # This file
│   ├── notebooks/                    # Notebook-specific guides
│   └── guides/                       # Additional guides
├── scripts/                          # Analysis and utility scripts
│   ├── analysis/                     # Analysis scripts
│   └── utils/                        # Utility functions
├── data/                             # Input data and cache
│   └── string_cache/                 # STRING database cache
├── results/                          # Analysis results
│   ├── functional_group_networks/    # Network visualizations
│   └── string_networks_by_group/     # STRING interaction data
├── *.ipynb                           # Analysis notebooks
└── requirements.txt                  # Python dependencies
```

## Key Findings

### Significantly Enriched Functional Groups (p < 0.05)
- **Structural**: 1.21x enrichment - proteins with structural roles are enriched in pLLPS
- **Other**: 1.04x enrichment - remaining proteins show slight enrichment
- **Depletion**: Ion Channels, Transporters, GPCRs, Enzymes show depletion in pLLPS

### Network Insights
- High pLLPS proteins preferentially interact with other high pLLPS proteins
- Structural proteins form densely connected clusters
- Ion channel networks show hub-and-spoke topology

## Performance

- Data loading: ~30 seconds
- Enrichment analysis: ~5 minutes
- STRING fetching: ~10-15 minutes (network dependent)
- Visualization: ~2-5 minutes
- **Total end-to-end runtime:** ~30 minutes

## Tools Used

- **Data Processing**: pandas, numpy, scipy
- **Network Analysis**: networkx, igraph
- **Visualization**: matplotlib, plotly, pyvis
- **Statistical Testing**: scipy.stats
- **Data Source**: STRING database, UniProt pLLPS scores

## References

- STRING Database: https://string-db.org/
- pLLPS Scores: https://uniprot.org/ (Phase Separation field)
- Publication: pLLPS: Prediction of Phase Separation Propensity
