# LLPS Protein Data Explorer

An interactive dashboard for exploring Liquid-Liquid Phase Separation (LLPS) protein data. Built with Streamlit for high traceability and easy sharing.

## Features

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
- 💾 **Data Export**: Download filtered data as CSV
- 📱 **Responsive Design**: Works on desktop and mobile browsers

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
| `n(DPR=> 25)` | Number of dipeptide repeats |
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

## Project Structure

```
mem_prot_llps/
├── app.py                    # Main Streamlit application
├── exploration_notebook.ipynb # Jupyter notebook for step-by-step analysis
├── requirements.txt          # Python dependencies
├── data/
│   └── sample_data.xlsx      # Sample dataset
└── README.md                 # This file
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