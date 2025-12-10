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

### Running Locally

```bash
streamlit run app.py
```

The dashboard will open in your default browser at `http://localhost:8501`.

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

### Command-Line Tools (Alternative):

For programmatic access or batch processing, you can also use the standalone modules:

```bash
# Simple interaction analysis
python pllps_interaction_simple.py

# Full network analysis with NetworkX
python string_interaction_analysis.py

# Location-based interaction analysis
python interaction_analysis.py
```

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
├── app.py                                       # Main Streamlit application with integrated interaction analysis
├── pllps_interaction_simple.py                  # Standalone simple interaction analysis
├── string_interaction_analysis.py               # Standalone full network analysis module
├── interaction_analysis.py                      # Location-based interaction analysis
├── pllps_interaction_analysis.ipynb             # Jupyter notebook for interaction analysis
├── exploration_notebook.ipynb                   # Jupyter notebook for step-by-step analysis
├── requirements.txt                             # Python dependencies
├── data/
│   ├── sample_data.xlsx                         # Sample dataset
│   └── Human Phase separation data.xlsx         # Full dataset (if available)
├── docs/
│   └── protein_interaction_analysis_exploration.md  # Analysis exploration document
└── README.md                                    # This file
```

## Jupyter Notebook

For a more detailed, step-by-step view of the data transformations, use the Jupyter notebook:

```bash
jupyter notebook exploration_notebook.ipynb
```

The notebook provides:
- **Transparent data transformations**: See exactly what happens at each stage
- **Interactive exploration**: Modify parameters and re-run cells
- **Detailed statistics**: Extended analysis beyond the dashboard
- **Export capabilities**: Save filtered data for further analysis

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is open source and available under the MIT License