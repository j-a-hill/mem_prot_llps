# Output Management Guide

## Problem
The notebook was generating outputs too large to display, causing "413 Request Entity Too Large" errors.

## Solution
Modified the notebook to:
1. **Save large datasets to CSV files** instead of displaying them
2. **Show only summaries** in the notebook output
3. **Use `.to_string()` for tables** instead of `display()` which creates HTML
4. **Limit displayed rows** to top 10-20 results

## Output Files Location
All detailed results are saved to the `results/` folder:

- `membrane_proteins_full.csv` - All membrane proteins with pLLPS scores
- `communities_analysis.csv` - All community detection results
- `high_pllps_proteins_connectivity.csv` - High pLLPS proteins with degree info
- `hub_analysis_full.csv` - Complete hub analysis (all proteins)
- `high_pllps_hubs.csv` - High pLLPS hub proteins (filtered)

## How to View Detailed Results

### Option 1: Load CSV files
```python
import pandas as pd

# Load any result file
df = pd.read_csv('results/high_pllps_hubs.csv')
df.head(20)  # View top 20
```

### Option 2: Query specific data
```python
# Example: Find specific protein
hubs = pd.read_csv('results/hub_analysis_full.csv')
hubs[hubs['Protein'] == 'DLG4']
```

### Option 3: Filter and analyze
```python
# Example: Get top hubs with pLLPS > 0.8
hubs = pd.read_csv('results/hub_analysis_full.csv')
top_hubs = hubs[(hubs['pLLPS'] > 0.8) & (hubs['Degree'] > 50)]
print(top_hubs.sort_values('Degree', ascending=False))
```

## Best Practices

1. **Always save large dataframes to CSV before displaying**
   ```python
   # Good ✓
   large_df.to_csv('results/output.csv', index=False)
   print(f"Saved {len(large_df)} rows")
   display(large_df.head(10))
   
   # Bad ✗
   display(large_df)  # Could be thousands of rows!
   ```

2. **Use print() instead of display() for tables**
   ```python
   # Good ✓ (plain text)
   print(df.head(10).to_string())
   
   # Bad ✗ (creates large HTML)
   display(df.head(10))
   ```

3. **Limit visualizations**
   - Don't create 100+ subplots
   - Show only top N most interesting cases
   - Save plots to files for later review

4. **Provide summaries**
   ```python
   print(f"Found {len(results)} proteins")
   print(f"Top 5: {results['protein'].head().tolist()}")
   ```

## Re-running the Notebook

If you need to re-run cells that generate CSV files, they will be automatically created in the `results/` folder. The summary cell at the bottom will show what files were created and their sizes.
