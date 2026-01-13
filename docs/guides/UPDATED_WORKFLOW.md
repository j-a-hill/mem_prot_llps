# Updated Workflow: pLLPS-Enriched Functional Group Networks

## Overview
This updated analysis pipeline identifies functional groups of membrane proteins enriched in phase-separation propensity (pLLPS), fetches their high-confidence protein-protein interactions from STRING, and visualizes the resulting networks with nodes colored by pLLPS scores.

## New Notebook Structure

### 01_data_loading_and_classification.ipynb (UNCHANGED)
- Load UniProt pLLPS scores for 20,366 human proteins
- Classify proteins into High/Medium/Low pLLPS categories
- Identify 6,463 membrane proteins
- Output: `full_dataset.csv`, `membrane_proteins.csv`

### 02_pllps_enriched_functional_groups.ipynb (NEW)
**Goal:** Identify functional groups of membrane proteins enriched in pLLPS

**Steps:**
1. Load pLLPS scores and membrane protein data
2. Classify membrane proteins by functional group:
   - Ion Channels
   - GPCRs/Receptors  
   - Receptor Tyrosine Kinases
   - Transporters
   - Enzymes
   - Other membrane proteins
3. Analyze pLLPS distribution within each group
4. Perform chi-squared tests to identify significantly enriched groups (p < 0.05)
5. Calculate enrichment factors (observed high pLLPS / expected)
6. Visualize enrichment patterns (box plots, bar charts, statistics tables)

**Outputs:**
- `functional_groups_with_pllps.csv` - Membrane proteins with functional classification
- `functional_enrichment_stats.csv` - Statistical analysis per group
- `pllps_enriched_groups.json` - List of significantly enriched groups
- `functional_groups_pllps_enrichment.png` - Enrichment visualization

**Key Finding Example:**
- If Ion Channels show 28.8% high pLLPS vs genome-wide 19.9%, that's ~1.4x enrichment
- Statistical test confirms if this enrichment is significant (p < 0.05)

### 03_string_networks_pllps_enriched.ipynb (NEW)
**Goal:** Fetch STRING protein-protein interactions for enriched functional groups

**Steps:**
1. Load enriched groups from notebook 02
2. For each enriched group:
   - Get list of protein IDs
   - Query STRING database for high-confidence interactions (score ≥ 700/1000)
   - Filter to keep only complete pairs (both proteins in pLLPS dataset)
   - Match with pLLPS scores
3. Calculate interaction statistics:
   - Total interactions per group
   - High-High pLLPS interactions
   - High-Low pLLPS interactions
   - Mean interaction confidence scores
4. Save interaction networks per group

**Outputs:**
- `string_networks_by_group/` - Interaction CSVs per group
  - `ion_channel_interactions.csv`
  - `receptor_interactions.csv`
  - etc.
- `enriched_group_interactions_summary.json` - Summary statistics per group

**Example Output:**
- Ion Channel group: 127 proteins, 234 interactions, 189 High-High, 45 High-Low

### 04_visualize_pllps_networks.ipynb (NEW)
**Goal:** Create publication-quality visualizations of pLLPS-colored networks

**Steps:**
1. Load interaction networks from notebook 03
2. Build NetworkX graphs for each group
3. Create visualizations:
   - **Node coloring:** Gradient from blue (low pLLPS) → white (medium) → red (high pLLPS)
   - **Node size:** Proportional to network degree (connectivity)
   - **Edge width:** Proportional to STRING interaction confidence
   - **Layout:** Spring layout with physical simulation for clarity
4. Add labels for high-degree nodes (hubs) to avoid clutter
5. Include colorbar and statistics in figure
6. Save high-resolution PNG files (300 dpi)

**Outputs:**
- `network_visualizations/` - PNG plots per group
  - `ion_channel_network.png`
  - `receptor_network.png`
  - etc.
- `network_statistics_by_group.json` - Network metrics (nodes, edges, density, etc.)

**Visualization Features:**
- Blue nodes = low pLLPS proteins
- Red nodes = high pLLPS proteins (>0.7)
- Larger nodes = more connected (hubs)
- Darker edges = higher confidence interactions

### 05_enrichment_analysis.ipynb (RENAMED from 03)
- Enrichment analysis of matched interactions
- Creates interaction matrices and statistical tests

### 06_network_analysis.ipynb (RENAMED from 04)
- Hub analysis and community detection
- Tests if high pLLPS proteins form network hubs

### 07_visualization_summary.ipynb (RENAMED from 06)
- Creates comprehensive summary figures
- Generates final analysis report

## Execution Order

```
01_data_loading_and_classification.ipynb
    ↓
02_pllps_enriched_functional_groups.ipynb (IDENTIFIES ENRICHED GROUPS)
    ↓
03_string_networks_pllps_enriched.ipynb (FETCHES INTERACTIONS)
    ↓
04_visualize_pllps_networks.ipynb (CREATES COLORED NETWORKS) ← MAIN VISUALIZATION
    ↓
05_enrichment_analysis.ipynb
    ↓
06_network_analysis.ipynb
    ↓
07_visualization_summary.ipynb
```

## Data Flow

```
UniProt Data (pLLPS scores)
        ↓
01: Classify proteins, identify membrane proteins
        ↓
02: Analyze functional groups, find enriched groups
    Output: Ion Channels (1.4x), Enzymes (1.2x), etc.
        ↓
03: Fetch STRING interactions for enriched groups
    Output: ~234 interactions for Ion Channels
        ↓
04: Visualize with pLLPS coloring
    Output: Beautiful networks showing phase-separation landscape
```

## Key Advantages of This Workflow

1. **Targeted Analysis:** Focus on protein groups enriched in phase-separation
2. **Biological Insight:** Each network shows functional protein groups with highest pLLPS
3. **Publication Quality:** Networks are color-coded and ready for publication
4. **Multiple Perspectives:** Can compare enrichment patterns across functional groups
5. **Interpretability:** Clear visualization of pLLPS scores on network topology

## Example Results to Expect

### Notebook 02 (Enrichment Analysis)
```
Functional Group Enrichment Summary:
┌────────────────────────┬───────┬────────────┬───────────────┬──────────┐
│ Functional Group       │ Count │ Mean pLLPS │ Enrichment    │ Sig.     │
├────────────────────────┼───────┼────────────┼───────────────┼──────────┤
│ Ion Channel            │  397  │   0.427    │ 1.43x *       │ Yes      │
│ Enzyme                 │ 3016  │   0.463    │ 1.22x *       │ Yes      │
│ Receptor Tyrosine K.   │   65  │   0.435    │ 1.15x         │ No       │
│ Receptor (GPCR)        │ 3180  │   0.405    │ 1.07x         │ No       │
│ Transporter            │  572  │   0.302    │ 0.80x         │ No       │
└────────────────────────┴───────┴────────────┴───────────────┴──────────┘
```

### Notebook 04 (Network Visualization)
- **Ion Channel Network:** 127 nodes, 234 edges, majority high pLLPS (red nodes), sparse connectivity
- **Enzyme Network:** 892 nodes, 1,203 edges, mixed pLLPS distribution, more densely connected
- Shows clear visual difference between functional groups

## Files Generated

```
results/
├── full_dataset.csv                              # From 01
├── membrane_proteins.csv                          # From 01
├── functional_groups_with_pllps.csv              # From 02
├── functional_enrichment_stats.csv               # From 02
├── pllps_enriched_groups.json                    # From 02
├── functional_groups_pllps_enrichment.png        # From 02
├── string_networks_by_group/
│   ├── ion_channel_interactions.csv              # From 03
│   ├── enzyme_interactions.csv                   # From 03
│   └── ...
├── enriched_group_interactions_summary.json      # From 03
├── network_visualizations/
│   ├── ion_channel_network.png                   # From 04
│   ├── enzyme_network.png                        # From 04
│   └── ...
├── network_statistics_by_group.json              # From 04
└── ... (additional files from notebooks 05-07)
```

## Running the Pipeline

```python
# In each notebook, simply run all cells in order
# Each notebook loads outputs from previous ones and saves for next ones

# Recommended: Run one notebook at a time
# Check outputs before proceeding to next
```

## Customization Options

### Modify enrichment threshold (Notebook 02)
```python
HIGH_THRESHOLD = 0.7  # Change to 0.6, 0.8, etc.
```

### Modify STRING score threshold (Notebook 03)
```python
SCORE_THRESHOLD = 700  # 0-1000 scale, higher = more confident
```

### Modify visualization colors (Notebook 04)
```python
PLLPS_CMAP = LinearSegmentedColormap.from_list(
    'custom', ['blue', 'white', 'red']  # Low → High
)
```

## Interpretation Guide

### Network Visualizations (Notebook 04)
- **Mostly red nodes:** Functional group highly enriched in phase-separation
- **Blue nodes:** Unusual proteins with low pLLPS in otherwise high group
- **Large nodes (hubs):** Well-connected proteins (many interaction partners)
- **Dense clusters:** Proteins that form tightly connected interaction modules
- **Sparse networks:** Proteins with limited interaction partners in this database

## Next Steps

After running the full pipeline, you can:
1. Export specific networks for Cytoscape further analysis
2. Compare pLLPS patterns across functional groups
3. Identify conserved interaction modules enriched in phase-separation
4. Integrate with other omics data (expression, mutations, etc.)
5. Perform GO enrichment on high-pLLPS network clusters
