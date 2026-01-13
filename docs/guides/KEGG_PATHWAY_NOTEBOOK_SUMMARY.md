# KEGG Pathway Analysis Notebook - Summary

## Overview
The `kegg_pathway_analysis.ipynb` notebook provides comprehensive analysis of pLLPS (Liquid-Liquid Phase Separation) scores in the context of KEGG biological pathways.

## What This Notebook Does

### 1. Data Loading and Classification
- Loads the LLPS protein dataset from `Human Phase separation data.xlsx`
- Classifies proteins as membrane proteins using pattern matching on:
  - Function [CC] annotations
  - Protein names
  - Subcellular location [CC]
- Categorizes proteins by pLLPS score:
  - **High**: p(LLPS) ≥ 0.7
  - **Medium**: 0.3 < p(LLPS) < 0.7
  - **Low**: p(LLPS) ≤ 0.3

### 2. KEGG Pathway Data Retrieval
- **Focus**: Fetches pathway data for **HIGH and LOW pLLPS membrane proteins only** (excludes medium)
  - Reduces API calls from ~6,300 to ~4,891 proteins by default
  - **NEW: Optional filtering** by function (e.g., ion channels) and location keywords
  - Can reduce to 50-500 proteins for targeted analyses (90%+ reduction)
  - Significantly reduces runtime (from hours to minutes)
- Fetches KEGG pathway annotations using the bioservices library
- Converts UniProt IDs to KEGG gene IDs using 3 fallback methods
- Maps proteins to all KEGG pathways they participate in
- **Rate Limiting**: Respects KEGG API limit of 3 requests/second (0.34s delay)
- Caches results to avoid redundant API calls

### 3. Pathway Enrichment Analysis
Answers the question: **Are particular pathways enriched or depleted in high/low pLLPS proteins?**

- Calculates enrichment ratios for each pathway:
  - High pLLPS enrichment ratio
  - Low pLLPS enrichment ratio
- Identifies pathways with:
  - Highest/lowest mean pLLPS scores
  - Most high pLLPS proteins
  - Most low pLLPS proteins
- Visualizes enrichment patterns with interactive plots

### 4. Score Similarity Analysis
Answers the question: **Do similar pLLPS scores interact in the same pathway?**

- Calculates within-pathway variance of pLLPS scores
- Compares to overall variance across all proteins
- Identifies pathways where proteins have:
  - Very similar pLLPS scores (low variance)
  - Diverse pLLPS scores (high variance)
- Visualizes score distributions with box plots

### 5. Pathway Visualization
Provides tools to **draw KEGG diagrams annotated with pLLPS scores**:

- Generates URLs to KEGG pathway diagrams
- Creates custom visualizations showing:
  - pLLPS scores for all proteins in a pathway
  - Color-coded bars (red=high, blue=low, gray=medium)
  - Direct links to KEGG web interface
- Allows selection of specific pathways for detailed analysis

### 6. Statistical Validation
Performs rigorous statistical tests:

- **Within-pathway similarity test**: Tests if proteins in the same pathway have more similar pLLPS scores than expected by chance
- **Correlation analysis**: Tests for correlation between pathway size and mean pLLPS
- **Mann-Whitney U test**: Compares pathway participation between high and low pLLPS proteins
- All tests include p-values and clear interpretation

### 7. Results Export
- Saves analysis results to Excel file (`kegg_pathway_analysis_results.xlsx`) with sheets:
  - Pathway statistics (enrichment, variance, counts)
  - High pLLPS membrane proteins
  - Low pLLPS membrane proteins

## Key Questions Answered

1. ✅ **Do similar pLLPS scores interact in the same pathway?**
   - Analyzed through within-pathway variance
   - Statistical test compares to null distribution

2. ✅ **Are particular pathways enriched or depleted in high/low pLLPS proteins?**
   - Calculates enrichment ratios vs. expected
   - Identifies top enriched/depleted pathways
   - Visualizes patterns with interactive plots

3. ✅ **Which pathways contain high or low pLLPS membrane proteins?**
   - Maps all membrane proteins to their pathways
   - Identifies pathway participation patterns
   - Exports lists for further investigation

4. ✅ **Can we visualize KEGG diagrams with pLLPS annotations?**
   - Generates pathway URLs for browser viewing
   - Creates custom pLLPS score visualizations
   - Provides protein-level details for each pathway

## Usage

### Prerequisites
- Python 3.9+
- Jupyter Notebook
- Internet access to KEGG REST API (rest.kegg.jp)

### Running the Notebook

```bash
# Install dependencies
pip install -r requirements.txt

# Launch Jupyter
jupyter notebook kegg_pathway_analysis.ipynb

# Run all cells in order
```

### Expected Runtime
- Pathway name retrieval: ~1 second
- **Default (no filtering)**: ~4,891 proteins, ~55 minutes
- **With function filtering** (e.g., ion channels): ~100-500 proteins, ~5-25 minutes
- **With combined filtering** (function + location): ~50-300 proteins, ~3-15 minutes
- All analysis and visualization: ~1-2 minutes

**Filtering dramatically reduces runtime** - configure filters based on your analysis needs.

### Customization Options
- **Filter by function category**: Ion channels, receptors, transporters, kinases, etc.
- **Filter by subcellular location**: Plasma membrane, ER, mitochondrion, etc.
- Adjust pLLPS thresholds (HIGH_PLLPS_THRESHOLD, LOW_PLLPS_THRESHOLD)
- Include medium pLLPS proteins (currently excluded to minimize API calls)
- Analyze all proteins instead of just membrane proteins
- Focus on specific pathways of interest
- Modify statistical tests and visualizations

**Recommended**: Use function/location filtering for exploratory analysis to minimize API calls.

## Output Files
- `kegg_pathway_analysis_results.xlsx`: Comprehensive results
- Interactive Plotly visualizations (displayed in notebook)
- KEGG pathway URLs for browser viewing

## Network Requirements
This notebook requires internet access to:
- KEGG REST API at rest.kegg.jp
- KEGG web interface for pathway diagrams

If running in a restricted environment, consider:
- Running locally with network access
- Using pre-cached pathway data
- Exporting protein lists for offline KEGG queries

## Troubleshooting UniProt to KEGG Mapping

### Expected Behavior
Not all UniProt IDs will successfully map to KEGG entries. This is **normal and expected**:

- **Typical mapping rate**: 30-60% of proteins successfully map to KEGG
- **Why proteins don't map**: KEGG focuses on metabolic and signaling pathways. Many proteins (structural, some regulatory, etc.) are not included in KEGG.

### Improved ID Conversion (v2)
The notebook now uses multiple conversion methods to maximize mapping success:

1. **Method 1**: Direct conversion with `uniprot:` prefix
   ```python
   kegg.conv('hsa', 'uniprot:P12345')
   ```

2. **Method 2**: Alternative with `up:` prefix (fallback)
   ```python
   kegg.conv('hsa', 'up:P12345')
   ```

3. **Method 3**: KEGG find service (second fallback)
   ```python
   kegg.find('genes', 'P12345')
   ```

### Diagnostic Output
The notebook now provides detailed mapping statistics:
- Number and percentage of proteins successfully mapped
- Number and percentage of proteins not found in KEGG
- Examples of unmapped proteins for verification
- Mean/median pathway counts for mapped proteins

### Validation
Check the diagnostic output after running the pathway retrieval section:
- ✅ If 30-60% map successfully: **Normal operation**
- ⚠️ If <10% map successfully: Check internet connectivity to rest.kegg.jp
- ❌ If 0% map successfully: Network access to KEGG API is blocked

## Integration with Existing Code
This notebook reuses:
- Membrane protein classification patterns from `app.py`
- pLLPS threshold conventions from `string_interaction_analysis.py`
- Data loading and parsing approaches from `exploration_notebook.ipynb`

## Future Enhancements
Potential additions:
- Integration with protein interaction networks
- Pathway-pathway similarity analysis
- Time series or condition-specific pathway analysis
- Machine learning for pathway-pLLPS predictions
- Automated pathway diagram annotation
