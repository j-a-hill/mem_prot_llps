# Notebook 05: Interactive STRING Interaction Networks Colored by pLLPS Score

## Overview

This new notebook creates exactly what you requested: **STRING interaction networks for functional groups (e.g., Ion channels) where each protein node is colored by its pLLPS score and edges represent STRING database interactions weighted by confidence score.**

## Key Features

### ✅ Node Coloring by pLLPS Score
- **Blue nodes** = Low pLLPS score (< 0.4) - proteins unlikely to phase separate
- **White nodes** = Medium pLLPS score (≈ 0.5) - proteins with moderate phase separation propensity
- **Red nodes** = High pLLPS score (> 0.7) - proteins likely to phase separate
- **Color scale** spans continuously from 0.0 to 1.0

### ✅ Edge Weighting by STRING Score
- **Edge thickness** represents STRING interaction confidence (0-1000 scale)
- **Edge opacity** also encodes confidence - thick, opaque edges = high confidence interactions
- Edges are drawn as gray lines connecting interacting proteins

### ✅ Node Size by Connectivity
- **Larger nodes** = more protein-protein interactions (higher degree)
- **Smaller nodes** = fewer interactions
- Helps identify hub proteins and central network components

### ✅ Multiple Output Formats

1. **Publication-Quality PNG (300 DPI)**
   - High-resolution image suitable for papers and presentations
   - Shows full network with pLLPS coloring, edge weights, and hub labels
   - File: `[group_name]_network.png`

2. **Interactive HTML Visualization (PyVis)**
   - Zoom, pan, and drag nodes to explore
   - Physics simulation for dynamic layout
   - Hover over nodes/edges for detailed information
   - File: `[group_name]_network_interactive.html`
   - Open in any web browser

3. **Network Data (CSV)**
   - **Nodes file**: Protein IDs, names, pLLPS scores, connectivity degree
   - **Edges file**: Protein pairs, STRING scores, pLLPS values for both proteins
   - **Summary stats**: Network metrics in JSON format

## How to Use

### Basic Usage
1. Open `05_interactive_functional_group_networks.ipynb`
2. Run cells 1-3 to install packages and load data
3. In **Cell 5**, change `TARGET_GROUP = 'Ion Channel'` to one of:
   - `'Ion Channel'`
   - `'GPCR'`
   - `'Transporter'`
   - `'Structural'`
   - `'Receptor Tyrosine Kinase'`
   - `'Enzyme'`
   - `'Other'`
4. Run remaining cells to generate visualizations

### Cell-by-Cell Breakdown

| Cell | Purpose | Output |
|------|---------|--------|
| 1-2 | Install packages | Library readiness |
| 3 | Import and configure | Ready environment |
| 4 | Load data | Display available groups |
| 5 | Select group | Statistics for chosen group |
| 6 | Load interactions | Check if pre-computed STRING data exists |
| 7 | Build network | Graph with nodes and edges |
| 8 | Prepare viz attributes | Node colors and sizes |
| 9 | Static PNG visualization | Publication-ready image |
| 10 | Interactive HTML visualization | Browser-explorable network |
| 11 | Export CSV data | Node, edge, and summary files |

## Current Status

The notebook has been tested with Ion Channels. The visualization shows:
- **334 Ion Channel proteins** (from the functional classification)
- **pLLPS scores** ranging from 0.090 to 1.000
- **81 high pLLPS proteins** (> 0.7) - shown in red
- Mean pLLPS of 0.428 for the group

Currently, no pre-computed STRING interactions exist for the Ion Channel group, so the network shows nodes without edges. To add edges, you need to either:

1. Run notebook **03_string_networks_pllps_enriched.ipynb** to fetch STRING interactions for the group, OR
2. Provide a pre-computed interactions CSV file at: `results/string_networks_by_group/ion_channel_interactions.csv`

## Visualization Interpretations

### Publication Figure Usage
The PNG outputs are publication-ready (300 DPI) and show:
- Clear pLLPS color gradient
- Node size as degree indicator
- Black edges for interactions
- Colorbar legend for pLLPS scores

### Interactive Exploration
The HTML files allow you to:
- **Zoom**: Scroll wheel or pinch gestures
- **Pan**: Click and drag the background
- **Drag nodes**: Click and drag individual proteins to explore connections
- **View details**: Hover over nodes to see pLLPS scores and protein names
- **Physics simulation**: Nodes repel/attract to find optimal layout

## Generated Files

For Ion Channels, the following files are created in `results/functional_group_networks/`:

```
ion_channel_network.png                      # Publication-quality visualization (1.8 MB)
ion_channel_network_interactive.html         # Interactive PyVis network (87 KB)
ion_channel_nodes.csv                        # Node attributes (13 KB)
ion_channel_edges.csv                        # Edge list (empty if no interactions)
ion_channel_summary_stats.json               # Network statistics
```

## Advantages Over Previous Notebooks

✅ **Focused visualization**: Each functional group gets its own dedicated network
✅ **Publication-ready output**: High-resolution PNG with proper color schemes
✅ **Interactive exploration**: PyVis HTML allows detailed exploration
✅ **Data export**: CSV files for further analysis in R, Cytoscape, etc.
✅ **pLLPS-centric coloring**: Primary focus on phase separation propensity
✅ **Customizable**: Easy to switch between functional groups

## Next Steps

1. **Generate STRING interactions** - Run notebook 03 to fetch interactions for your desired functional groups
2. **Create networks with edges** - Re-run notebook 05 after interaction data is available
3. **Analyze hub proteins** - Use the node CSV to identify high-degree, high-pLLPS proteins
4. **Export for further analysis** - Use the CSV files in Cytoscape, R, or other network tools

## Notes

- The colormap (blue→white→red) matches the pLLPS probability scale (0.0→0.5→1.0)
- Edge weights are normalized to the range [0.5, 4.5] pixels for visibility
- The spring layout algorithm uses seed=42 for reproducibility
- Hub labels (degree ≥ 3) are shown to highlight key proteins
- All data is normalized to the 0-1 pLLPS scale
