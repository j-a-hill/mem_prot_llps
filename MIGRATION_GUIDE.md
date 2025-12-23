# Migration Guide: Using llps_functions.py

## Overview

All analysis functions have been consolidated into a single master file: `llps_functions.py`

This consolidation provides:
- ✅ **Easier maintenance**: All functions in one place
- ✅ **Better documentation**: Comprehensive docstrings for all functions
- ✅ **No duplication**: Removed duplicate code across modules
- ✅ **Backward compatibility**: Old modules still work with deprecation warnings

## Quick Start

### In Python Scripts

```python
# New way (recommended)
from llps_functions import (
    load_llps_data,
    fetch_string_interactions,
    match_interactions_to_pllps,
    analyze_interaction_enrichment
)

# Or import all
import llps_functions as lf
df = lf.load_llps_data('data.xlsx')
```

### In Jupyter Notebooks

```python
# Import master functions module
import llps_functions as lf

# Use with lf. prefix
df = lf.load_llps_data('Human Phase separation data.xlsx')
df_classified, high_ids = lf.get_high_pllps_proteins(df, threshold=0.7)
interactions_df, errors = lf.fetch_string_interactions(high_ids[:50])
```

## Migration from Old Modules

### From string_functions.py

```python
# Old
from string_functions import fetch_string_interactions

# New
from llps_functions import fetch_string_interactions
```

### From pllps_interaction_simple.py

```python
# Old
from pllps_interaction_simple import load_and_filter_high_pllps

# New
from llps_functions import get_high_pllps_proteins
```

### From string_interaction_analysis.py

```python
# Old
from string_interaction_analysis import analyze_network

# New
from llps_functions import analyze_network
```

### From interaction_analysis.py

```python
# Old
from interaction_analysis import parse_location, analyze_interaction_matrix

# New
from llps_functions import parse_location, analyze_interaction_matrix
```

## Available Functions

### Data Loading and Classification
- `load_llps_data()` - Load protein data from Excel
- `load_and_classify_data()` - Load and classify proteins into High/Medium/Low
- `get_high_pllps_proteins()` - Identify high pLLPS proteins

### Location Analysis
- `parse_location()` - Parse subcellular location strings
- `add_location_columns()` - Add parsed location columns to DataFrame
- `analyze_interactions_by_location()` - Analyze interactions by location

### STRING Database
- `get_string_mapping()` - Map UniProt IDs to STRING IDs
- `fetch_string_interactions()` - Fetch interactions from STRING API
- `fetch_interaction_partners()` - Fetch interaction partners
- `load_string_network_file()` - Load STRING network from TSV file

### Network Analysis
- `match_interactions_to_pllps()` - Match interactions to pLLPS data
- `match_interactors_to_pllps()` - Alternative matching function
- `analyze_network()` - Full network analysis with NetworkX

### Enrichment Analysis
- `analyze_interaction_enrichment()` - Chi-squared enrichment test
- `analyze_interaction_matrix()` - 3x3 enrichment matrix analysis

### Visualization
- `plot_interaction_heatmap()` - Plot enrichment heatmap
- `plot_location_heatmaps()` - Plot location-specific heatmaps
- `print_analysis_report()` - Print formatted analysis report

### Export and Caching
- `save_interactions_to_cache()` - Save interactions for offline use
- `export_protein_list()` - Export protein list for STRING web interface

## Deprecated Modules

The following modules are **deprecated** and will be removed in future versions:

- ⚠️ `string_functions.py` → Use `llps_functions.py`
- ⚠️ `pllps_interaction_simple.py` → Use `llps_functions.py`
- ⚠️ `string_interaction_analysis.py` → Use `llps_functions.py`
- ⚠️ `interaction_analysis.py` → Use `llps_functions.py`

These files still work but will show deprecation warnings. Please update your imports.

## Examples

### Example 1: Load Data and Classify

```python
import llps_functions as lf

# Load data
df = lf.load_llps_data('Human Phase separation data.xlsx')

# Get high pLLPS proteins
df_classified, high_ids = lf.get_high_pllps_proteins(df, threshold=0.7)
print(f"Found {len(high_ids)} high pLLPS proteins")
```

### Example 2: Fetch STRING Interactions

```python
import llps_functions as lf

# Fetch interactions
protein_ids = ['P04637', 'P38398', 'P51587']
interactions_df, errors = lf.fetch_string_interactions(
    protein_ids,
    score_threshold=700,
    progress_callback=lambda msg: print(msg)
)
```

### Example 3: Enrichment Analysis

```python
import llps_functions as lf
import pandas as pd

# Load data
df = lf.load_llps_data('data.xlsx')
interactions_df, _ = lf.fetch_string_interactions(df['Entry'].head(100))

# Match and analyze
matched_df = lf.match_interactions_to_pllps(interactions_df, df)
results = lf.analyze_interaction_enrichment(matched_df, threshold=0.7)

if results and results['p_value'] < 0.05:
    print(f"Significant enrichment: {results['enrichment']:.2f}x")
```

## Support

For issues or questions:
1. Check function docstrings: `help(lf.function_name)`
2. See `example_string_usage.py` for complete examples
3. Refer to notebooks for workflow examples
4. Open an issue on GitHub

## Version History

- **v2.0** (2025-12): Consolidated all functions into `llps_functions.py`
- **v1.0** (2024): Original modular structure
