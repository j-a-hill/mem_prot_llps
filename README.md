# LLPS Protein Data Explorer

An interactive dashboard and analysis workflow for exploring Liquid-Liquid Phase Separation (LLPS) protein data, with integrated STRING network visualization. Built with Jupyter notebooks and a Shiny for Python dashboard.

## Quick Start

### Analysis Notebooks (Recommended)

Run the Jupyter notebooks in sequence for comprehensive analysis:

1. **[01_data_loading_and_classification.ipynb](01_data_loading_and_classification.ipynb)** – Load and classify pLLPS data
2. **[02_pllps_enriched_functional_groups.ipynb](02_pllps_enriched_functional_groups.ipynb)** – Identify functionally enriched groups
3. **[03_string_networks_pllps_enriched.ipynb](03_string_networks_pllps_enriched.ipynb)** – Fetch STRING interactions for enriched groups
4. **[04_visualize_pllps_networks.ipynb](04_visualize_pllps_networks.ipynb)** – Visualise pLLPS-coloured networks
5. **[05_interactive_functional_group_networks.ipynb](05_interactive_functional_group_networks.ipynb)** – Create detailed interactive networks
6. **[06_pllps_scores_analysis.ipynb](06_pllps_scores_analysis.ipynb)** – Deep-dive pLLPS score analysis

📖 See [docs/ANALYSIS_WORKFLOW.md](docs/ANALYSIS_WORKFLOW.md) for detailed workflow documentation.

### Interactive Dashboard (Alternative)

For quick data exploration, launch the Shiny for Python dashboard:
```bash
shiny run scripts/shiny_app.py --reload --port 8000
```

### Shinylive Dashboard (Browser-Based, No Server Required)

The `dashboard/` directory contains a self-contained [Shinylive](https://shiny.posit.co/py/docs/shinylive.html)
app that runs entirely in the browser via WebAssembly. No Python server is needed.

**Run locally:**
```bash
pip install shinylive pandas plotly numpy openpyxl
shinylive export dashboard/ /tmp/llps_dashboard
python3 -m http.server --directory /tmp/llps_dashboard 8008
# open http://localhost:8008
```

**Features:** p(LLPS) / length / location / function filters, distribution histograms,
scatter plots, location and function bar charts, and CSV export.

The dashboard is also automatically deployed to GitHub Pages via the
`.github/workflows/deploy-dashboard.yml` workflow on every push to `main`.

## Features

### 📊 Data Explorer
- **Data Upload**: Upload your XLSX files containing protein LLPS data
- **Search & Filter**: Search proteins by name, entry ID, or keywords with customisable filters
- **Interactive Visualisations**: distribution plots, scatter plots, subcellular location analysis, functional category analysis
- **Functional Category Classification**: Automatically categorise proteins by function (ion channels, GPCRs, kinases, transporters, etc.) using the YAML-based rules in `data/functional_classification_terms.yaml`

### 🔗 Protein Interaction Analysis
- **STRING Integration**: Fetch protein-protein interactions directly from the STRING database
- **Enrichment Analysis**: Test whether high pLLPS proteins preferentially interact with each other via chi-squared tests
- **Customisable Parameters**: adjustable pLLPS threshold, STRING confidence score, and sample size
- **Data Export**: Download interaction data with pLLPS annotations

## Installation

1. Clone the repository:
```bash
git clone https://github.com/j-a-hill/mem_prot_llps.git
cd mem_prot_llps
```

2. Create a virtual environment (recommended):
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Project Structure

```
mem_prot_llps/
├── 01_data_loading_and_classification.ipynb       # Step 1: load & classify proteins
├── 02_pllps_enriched_functional_groups.ipynb      # Step 2: identify enriched groups
├── 03_string_networks_pllps_enriched.ipynb        # Step 3: fetch STRING interactions
├── 04_visualize_pllps_networks.ipynb              # Step 4: visualise networks
├── 05_interactive_functional_group_networks.ipynb # Step 5: interactive network explorer
├── 06_pllps_scores_analysis.ipynb                 # Step 6: pLLPS score deep-dive
├── protein_pllps_lookup.ipynb                     # Utility: look up individual proteins
├── llps_functions.py                              # Backward-compatible re-export shim
├── llps/                                          # Core analysis library
│   ├── constants.py                               # STRING API constants, StringQueryConfig
│   ├── data.py                                    # Data loading & pLLPS classification
│   ├── location.py                                # Subcellular location parsing
│   ├── string_api.py                              # STRING database API queries
│   ├── network.py                                 # NetworkX network analysis
│   ├── enrichment.py                              # Interaction enrichment analysis
│   ├── visualization.py                           # Heatmap plotting & reports
│   ├── io.py                                      # Caching, export, result save/load
│   └── functional.py                              # Functional category classification
├── tests/                                         # Unit tests (pytest)
├── scripts/
│   ├── shiny_app.py                               # Interactive Shiny dashboard
│   └── analysis/
│       ├── generate_all_networks.py               # Batch network generation
│       └── generate_string_cache.py               # Pre-cache STRING data offline
├── data/
│   ├── Human Phase separation data.xlsx           # Source dataset
│   ├── sample_data.xlsx                           # Small sample for testing
│   └── functional_classification_terms.yaml       # Functional category rules (YAML)
├── results/                                       # Generated analysis outputs (CSV/JSON/PNG)
├── docs/                                          # Documentation
│   ├── ANALYSIS_WORKFLOW.md                       # Full workflow walkthrough
│   └── guides/                                    # Additional reference guides
├── deprecated/                                    # Archived old notebooks
├── requirements.txt                               # Python dependencies
└── pyproject.toml                                 # Package config & tool settings
```

## Programmatic Usage

All analysis functions are available via the `llps` package (or the backward-compatible `llps_functions` module):

```python
from llps_functions import (
    load_llps_data,
    fetch_string_interactions,
    match_interactions_to_pllps,
    analyze_interaction_enrichment,
    StringQueryConfig,
)

# Load and classify data
df = load_llps_data('data/Human Phase separation data.xlsx')

# Fetch STRING interactions
protein_ids = ['P04637', 'P38398', 'P51587']
cfg = StringQueryConfig(score_threshold=700)
interactions_df, errors = fetch_string_interactions(protein_ids, config=cfg)

# Match interactions to pLLPS scores
matched_df = match_interactions_to_pllps(interactions_df, df)

# Analyse enrichment
results = analyze_interaction_enrichment(matched_df, threshold=0.7)
```

## Data Format

The notebooks and dashboard expect an XLSX file with the following columns:

| Column | Description |
|--------|-------------|
| `Entry` | UniProt entry ID |
| `Entry name` | UniProt entry name |
| `Protein names` | Full protein names |
| `p(LLPS)` | Probability of LLPS (0–1) |
| `Length` | Protein sequence length |
| `Function [CC]` | Function annotation |
| `Subcellular location [CC]` | Subcellular location |

A sample dataset is included at `data/sample_data.xlsx`.

## STRING Interaction Caching

For environments without network access to `string-db.org`, pre-generate a cache file:

```bash
python scripts/analysis/generate_string_cache.py --threshold 0.7 --score 700 --max-proteins 500
# Saves to: data/string_cache_700.json
```

The dashboard and notebooks will automatically use the cached file if it is present.

## Running Tests

```bash
pip install pytest
pytest tests/
```

## Contributing

Contributions are welcome. Please open a pull request.

## License

This project is open source and available under the MIT License.
