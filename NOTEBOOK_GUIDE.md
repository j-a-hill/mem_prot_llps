# Notebook Analysis Guide

## Overview

The LLPS protein interaction analysis has been split into focused, manageable notebooks to avoid memory crashes and improve usability. Each notebook handles a specific stage of the analysis and saves results for downstream use.

## Benefits of the New Structure

✅ **Smaller file sizes** - Each notebook is <100KB vs 18MB for the old notebook  
✅ **Faster loading** - Notebooks load instantly without embedded outputs  
✅ **No crashes** - Results are saved to disk instead of embedded in notebooks  
✅ **Reusable results** - Load saved results without recomputing  
✅ **Modular workflow** - Run only the analyses you need  
✅ **Easier debugging** - Isolate issues to specific analysis stages  

## Notebook Workflow

### 01: Data Loading and Classification
**File:** `01_data_loading_and_classification.ipynb`

**Purpose:** Load and classify the LLPS dataset

**What it does:**
- Loads the Human Phase Separation dataset
- Classifies proteins into High/Medium/Low pLLPS categories
- Identifies membrane proteins
- Generates initial visualizations

**Outputs:**
```
results/full_dataset.csv              # Complete classified dataset
results/high_pllps_proteins.csv       # High pLLPS proteins (≥0.7)
results/membrane_proteins.csv         # Membrane proteins only
results/classification_summary.json   # Summary statistics
results/pllps_distribution.png        # Visualization
```

**Runtime:** ~30 seconds

---

### 02: STRING Interactions
**File:** `02_string_interactions.ipynb`

**Purpose:** Fetch protein-protein interactions from STRING database

**What it does:**
- Loads high pLLPS proteins from step 01
- Fetches interactions from STRING API (or loads cached data)
- Matches interactions to the pLLPS dataset
- Visualizes interaction patterns

**Inputs:**
```
results/high_pllps_proteins.csv
results/full_dataset.csv
```

**Outputs:**
```
results/string_interactions_raw.csv      # Raw STRING interactions
results/string_interactions_matched.csv  # Matched to pLLPS data
results/interaction_summary.json         # Summary statistics
results/string_interaction_stats.png     # Visualization
```

**Runtime:** ~5-10 minutes (first run with API calls), ~30 seconds (cached)

**Note:** The first run makes API calls to STRING. Results are cached to avoid re-fetching.

---

### 03: Enrichment Analysis
**File:** `03_enrichment_analysis.ipynb`

**Purpose:** Test if high pLLPS proteins preferentially interact

**What it does:**
- Tests for high-high pLLPS interaction enrichment
- Performs chi-squared statistical tests
- Creates 3x3 interaction matrices (High/Med/Low)
- Generates heatmap visualizations

**Inputs:**
```
results/string_interactions_matched.csv
results/full_dataset.csv
```

**Outputs:**
```
results/enrichment_results.json           # Statistical test results
results/enrichment_matrix_observed.csv    # Observed interactions
results/enrichment_matrix_expected.csv    # Expected interactions
results/enrichment_matrix_ratios.csv      # Enrichment ratios
results/enrichment_analysis.png           # Visualization
results/interaction_matrix_heatmaps.png   # Heatmaps
```

**Runtime:** ~30 seconds

---

### 04: Network Analysis
**File:** `04_network_analysis.ipynb`

**Purpose:** Analyze protein interaction networks

**What it does:**
- Builds network graph using NetworkX
- Identifies network hubs (highly connected proteins)
- Tests if high pLLPS proteins are hubs
- Performs community detection
- Analyzes community characteristics

**Inputs:**
```
results/string_interactions_matched.csv
results/membrane_proteins.csv
```

**Outputs:**
```
results/hub_analysis_results.csv       # Hub proteins and metrics
results/community_analysis_results.csv # Community structure
results/network_stats.json             # Network statistics
results/hub_analysis.png               # Visualization
results/community_analysis.png         # Visualization
```

**Runtime:** ~1-2 minutes

---

### 05: Functional Groups
**File:** `05_functional_groups.ipynb`

**Purpose:** Analyze patterns by functional category

**What it does:**
- Classifies proteins into functional categories (Ion Channels, GPCRs, etc.)
- Analyzes pLLPS patterns within categories
- Examines interaction patterns between categories
- Visualizes functional group characteristics

**Inputs:**
```
results/full_dataset.csv
results/string_interactions_matched.csv
```

**Outputs:**
```
results/functional_categories.csv              # Proteins with annotations
results/functional_group_stats.json            # Statistics per group
results/functional_group_interaction_matrix.csv # Interaction matrix
results/functional_category_distribution.png    # Visualization
results/pllps_by_functional_group.png          # Visualization
```

**Runtime:** ~1 minute

---

### 06: Visualization Summary
**File:** `06_visualization_summary.ipynb`

**Purpose:** Generate comprehensive summary visualizations

**What it does:**
- Loads all previous analysis results
- Creates publication-ready multi-panel figures
- Generates text summary report
- Lists all generated visualizations

**Inputs:** All result files from previous notebooks

**Outputs:**
```
results/comprehensive_summary.png      # Multi-panel summary figure
results/analysis_summary_report.txt    # Text report
```

**Runtime:** ~30 seconds

---

## How to Use

### Quick Start

Run notebooks in order:

```bash
jupyter notebook 01_data_loading_and_classification.ipynb
# Run all cells, wait for completion

jupyter notebook 02_string_interactions.ipynb
# Run all cells, wait for completion

jupyter notebook 03_enrichment_analysis.ipynb
# Run all cells

jupyter notebook 04_network_analysis.ipynb
# Run all cells

jupyter notebook 05_functional_groups.ipynb
# Run all cells

jupyter notebook 06_visualization_summary.ipynb
# Run all cells for final summary
```

### Re-running Analysis

If you want to re-run a specific stage:

1. **Starting fresh:** Delete the `results/` directory and run all notebooks from 01
2. **Re-run from step N:** Delete only the outputs from step N onwards, then run from that notebook
3. **Update only visualizations:** Just re-run notebook 06

### Loading Saved Results

To load and explore saved results in a new notebook:

```python
import llps_functions as lf

# Load any saved result
df = lf.load_analysis_result('full_dataset', format='csv')
stats = lf.load_analysis_result('enrichment_results', format='json')

# List all available results
lf.list_saved_results()
```

## Result File Reference

### CSV Files (DataFrames)
- `full_dataset.csv` - Complete protein dataset with classifications
- `high_pllps_proteins.csv` - Subset of high pLLPS proteins
- `membrane_proteins.csv` - Subset of membrane proteins
- `string_interactions_raw.csv` - Raw STRING interactions
- `string_interactions_matched.csv` - Interactions with pLLPS annotations
- `hub_analysis_results.csv` - Network hub analysis
- `community_analysis_results.csv` - Community detection results
- `functional_categories.csv` - Proteins with functional annotations
- `enrichment_matrix_*.csv` - Interaction matrices

### JSON Files (Summaries/Stats)
- `classification_summary.json` - Dataset classification statistics
- `interaction_summary.json` - STRING interaction statistics
- `enrichment_results.json` - Enrichment analysis results
- `network_stats.json` - Network topology statistics
- `functional_group_stats.json` - Functional category statistics

### PNG Files (Visualizations)
- All visualization plots saved with descriptive names

## Comparison with Old Notebooks

### Old Structure (DEPRECATED)
```
pllps_interaction_analysis.ipynb    # 18MB - Everything in one file
exploration_notebook.ipynb           # 1.1MB - General exploration
```

**Problems:**
- ❌ Huge file sizes cause crashes
- ❌ All outputs embedded in notebook
- ❌ Difficult to navigate
- ❌ Must re-run everything to update one section
- ❌ Hard to find specific analyses

### New Structure
```
01_data_loading_and_classification.ipynb    # <50KB each
02_string_interactions.ipynb
03_enrichment_analysis.ipynb
04_network_analysis.ipynb
05_functional_groups.ipynb
06_visualization_summary.ipynb
```

**Benefits:**
- ✅ Small, fast-loading notebooks
- ✅ No embedded outputs (saved to results/)
- ✅ Clear, focused analyses
- ✅ Modular - run only what you need
- ✅ Easy to find and understand

## Tips and Best Practices

1. **Always run in order first time** - Each notebook depends on outputs from previous ones

2. **Check results directory** - Use `lf.list_saved_results()` to see what's available

3. **Don't embed large outputs** - The new notebooks are designed to save results to disk

4. **Clear outputs before committing** - Use "Cell → All Output → Clear" to keep notebooks small

5. **Use caching** - STRING interactions are cached automatically to avoid re-fetching

6. **Customize parameters** - Each notebook has configuration cells at the top (thresholds, etc.)

7. **Save frequently** - Each notebook saves its results, so you won't lose progress

## Migration from Old Notebooks

If you have existing code in the old notebooks:

1. **Find the analysis stage** - Identify which new notebook covers that analysis
2. **Extract the code** - Copy relevant cells to the appropriate new notebook
3. **Use llps_functions** - Replace inline functions with calls to `llps_functions.py`
4. **Save results** - Use `lf.save_analysis_result()` to persist outputs
5. **Update downstream** - Ensure subsequent notebooks load your new results

## Support

For questions or issues:
1. Check function docstrings: `help(lf.function_name)`
2. See `example_string_usage.py` for code examples
3. Refer to `MIGRATION_GUIDE.md` for function references
4. Open an issue on GitHub

## Version History

- **2025-01**: Split large notebook into 6 focused notebooks
- **2024-12**: Consolidated functions into `llps_functions.py`
- **2024**: Original analysis in single large notebook
