# LLPS Protein Data Explorer

An interactive dashboard for exploring Liquid-Liquid Phase Separation (LLPS) protein data. Built with Streamlit for high traceability and easy sharing.

## Features

### 📊 Data Explorer
- 📁 **Data Upload**: Upload your XLSX files containing protein LLPS data
- 🔍 **Search & Filter**: Search proteins by name, entry ID, or keywords with customizable filters
- 📊 **Interactive Visualizations**: 
  - Distribution plots for p(LLPS) and protein length
  - Scatter plots with customizable axes
  - Subcellular location analysis
  - Functional category analysis
  - Correlation analysis
- 🏷️ **Functional Category Classification**: Automatically categorize proteins by function using pattern matching on Function [CC] and Protein names columns. Supported categories include:
  - Membrane proteins (Multi-pass, Single-pass)
  - Ion channels (Voltage-gated, Ligand-gated)
  - Receptors (G protein-coupled, general)
  - Enzymes (Kinase, Phosphatase, Protease, etc.)
  - Transporters (ABC transporter, general)
  - DNA/RNA binding proteins
  - And many more...

### 🔗 Protein Interaction Analysis (NEW!)
- **STRING Database Integration**: Fetch protein-protein interactions directly from STRING database
- **Enrichment Analysis**: Test whether high pLLPS proteins preferentially interact with each other
- **Statistical Testing**: Chi-squared tests to determine significance of interaction patterns
- **Interactive Visualizations**: 
  - Observed vs Expected interaction distributions
  - Interaction network scatter plots
  - Enrichment factor calculations
- **Customizable Parameters**:
  - Adjustable pLLPS threshold for classification
  - STRING confidence score filtering (medium, high, highest)
  - Sample size control to limit API calls
- **Data Export**: Download interaction data for further analysis

### 💾 Data Export
- Download filtered protein data as CSV
- Export interaction data with pLLPS annotations
- Column information and metadata

## Installation

1. Clone this repository:
```bash
git clone https://github.com/j-a-hill/mem_prot_llps.git
cd mem_prot_llps
```

2. Create a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Interactive Web Application

**Shiny for Python Dashboard:**

```bash
shiny run shiny_app.py --reload --port 8000
```

The dashboard will open in your default browser at `http://localhost:8000`.

### Programmatic Usage

**Using LLPS Functions in Python:**

```python
from llps_functions import (
    fetch_string_interactions,
    match_interactions_to_pllps,
    analyze_interaction_enrichment
)

# Fetch interactions
protein_ids = ['P04637', 'P38398', 'P51587']
interactions_df, errors = fetch_string_interactions(protein_ids, score_threshold=700)

# Match to your dataset
matched_df = match_interactions_to_pllps(interactions_df, pllps_df)

# Analyze enrichment
results = analyze_interaction_enrichment(matched_df, threshold=0.7)
```

See `example_string_usage.py` for complete examples.

### Data Format

The dashboard expects an XLSX file with the following columns:

| Column | Description |
|--------|-------------|
| `Entry` | UniProt entry ID |
| `Entry name` | UniProt entry name |
| `Protein names` | Full protein names |
| `p(LLPS)` | Probability of LLPS (0-1) |
| `n(DPR=> 25)` | Number of droplet promoting regions (DPR) |
| `Length` | Protein sequence length |
| `Function [CC]` | Function annotation |
| `Subcellular location [CC]` | Subcellular location |
| `Involvement in disease` | Disease associations |
| `Cross-reference (PDB)` | PDB structure references |

### Sample Data

A sample dataset is included in `data/sample_data.xlsx` for testing purposes.

## Deployment Options

### Streamlit Community Cloud (Free)

1. Push your code to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect your GitHub repository
4. Deploy with one click

### Other Deployment Options

- **Docker**: Build and deploy as a container
- **Heroku**: Deploy using the Streamlit buildpack
- **AWS/GCP/Azure**: Deploy on cloud VMs

## Protein Interaction Analysis

The interactive webapp now includes integrated protein interaction analysis! Access it through the **🔗 Protein Interactions** tab in the main dashboard.

### Web Interface Features:
- **Interactive Analysis**: Point-and-click interface for fetching interactions from STRING
- **Real-time Enrichment**: Immediate statistical analysis of interaction patterns
- **Visual Results**: Interactive charts showing observed vs expected interactions
- **Export Capabilities**: Download interaction data for external analysis
- **Caching Support**: Pre-cache STRING data for offline/restricted environments

### STRING Interaction Caching

For deployments without network access to string-db.org (e.g., restricted environments), you can pre-generate a cache:

```bash
# Generate cache file (requires network access)
python generate_string_cache.py --threshold 0.7 --score 700 --max-proteins 500

# Cache will be saved to: data/string_cache_700.json
```

The Shiny app will automatically use cached data if:
1. A cache file exists in `data/string_cache_{score}.json`
2. Network access to STRING is unavailable

**Note**: Always generate the cache file locally with network access before deploying to restricted environments.

### Command-Line Tools (Alternative):

For programmatic access or batch processing, you can also use the standalone modules:

```bash
# Simple interaction analysis (using master functions)
python pllps_interaction_simple.py

# Full network analysis with NetworkX (using master functions)
python string_interaction_analysis.py

# Location-based interaction analysis (using master functions)
python interaction_analysis.py
```

**Note:** All standalone scripts now import from the consolidated `llps_functions.py` module.

These modules provide:
- **STRING API integration**: Fetch protein-protein interactions
- **Network analysis**: Analyze interaction patterns using NetworkX
- **Hub detection**: Identify if high pLLPS proteins are network hubs
- **Cluster analysis**: Detect communities of interacting proteins
- **Enrichment analysis**: Test if high pLLPS proteins preferentially interact

See `docs/protein_interaction_analysis_exploration.md` for detailed documentation on:
- Available APIs (STRING, BioGRID, IntAct, Reactome)
- Analysis strategies and hypotheses
- Implementation approaches
- Expected outcomes and interpretations

## Project Structure

```
mem_prot_llps/
├── shiny_app.py                                 # Main Shiny application with integrated interaction analysis
├── llps_functions.py                            # 🆕 MASTER functions module (all functions consolidated here)
├── example_string_usage.py                      # Example usage of llps_functions
├── generate_string_cache.py                     # Cache generator for offline use
├── pllps_interaction_analysis.ipynb             # Jupyter notebook for interaction analysis
├── exploration_notebook.ipynb                   # Jupyter notebook for step-by-step analysis
├── kegg_pathway_analysis.ipynb                  # KEGG pathway analysis notebook
├── requirements.txt                             # Python dependencies
├── interaction_analysis.py                      # ⚠️  DEPRECATED - Use llps_functions.py
├── string_functions.py                          # ⚠️  DEPRECATED - Use llps_functions.py
├── pllps_interaction_simple.py                  # ⚠️  DEPRECATED - Use llps_functions.py
├── string_interaction_analysis.py               # ⚠️  DEPRECATED - Use llps_functions.py
├── data/
│   ├── sample_data.xlsx                         # Sample dataset
│   └── Human Phase separation data.xlsx         # Full dataset (if available)
├── docs/
│   └── protein_interaction_analysis_exploration.md  # Analysis exploration document
└── README.md                                    # This file
```

**Important Note:** All analysis functions have been consolidated into `llps_functions.py` for easier maintenance and reusability. The old module files (`interaction_analysis.py`, `string_functions.py`, `pllps_interaction_simple.py`, `string_interaction_analysis.py`) are kept for backward compatibility but are deprecated and will be removed in future versions.

**Important Note:** All analysis functions have been consolidated into `llps_functions.py` for easier maintenance and reusability. The old module files (`interaction_analysis.py`, `string_functions.py`, `pllps_interaction_simple.py`, `string_interaction_analysis.py`) are kept for backward compatibility but are deprecated and will be removed in future versions.

## Jupyter Notebooks

### General Data Exploration

For a more detailed, step-by-step view of the data transformations, use the Jupyter notebook:

```bash
jupyter notebook exploration_notebook.ipynb
```

The notebook provides:
- **Transparent data transformations**: See exactly what happens at each stage
- **Interactive exploration**: Modify parameters and re-run cells
- **Detailed statistics**: Extended analysis beyond the dashboard
- **Export capabilities**: Save filtered data for further analysis

### KEGG Pathway Analysis

Analyze pLLPS scores in the context of KEGG biological pathways:

```bash
jupyter notebook kegg_pathway_analysis.ipynb
```

This notebook provides:
- **Membrane Protein Classification**: Identify high/low pLLPS membrane proteins
- **Pathway Enrichment Analysis**: Discover pathways enriched/depleted in high/low pLLPS proteins
- **Score Similarity Analysis**: Test whether similar pLLPS scores co-occur in pathways
- **Pathway Visualization**: Generate KEGG pathway diagrams annotated with pLLPS scores
- **Statistical Validation**: Statistical tests to validate observed patterns
- **Export Results**: Save pathway analysis results for further investigation

**Key Features:**
- Integrates with KEGG REST API (respects 3 requests/second rate limit)
- Reuses function parser to identify membrane proteins
- Calculates pathway enrichment ratios
- Visualizes pLLPS score distributions within pathways
- Provides direct links to KEGG pathway diagrams

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is open source and available under the MIT License