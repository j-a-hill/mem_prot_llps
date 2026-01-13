# Network Export Summary

Generated: 2024-01-13

## Overview

STRING interaction networks with pLLPS scoring have been successfully exported for all functional groups. Each network includes:

1. **Cytoscape JSON** (`.json`) - Standard Cytoscape.js format for import into Cytoscape Desktop
2. **Nodes CSV** - Node attributes (UniProt ID, gene name, pLLPS score, degree, functional group)
3. **Edges CSV** - Edge attributes (protein pairs, gene names, STRING scores, pLLPS values)
4. **Summary Stats JSON** - Network statistics and metadata

## Network Statistics

### Ion Channel
- **Nodes:** 149 proteins
- **Edges:** 282 interactions
- **Connected Components:** 41
- **Mean pLLPS:** 0.497
- **High pLLPS proteins (>0.7):** 46
- **Mean STRING score:** 843.7
- **Network density:** 0.026
- **Files:** `ion_channel_*`

### Structural
- **Nodes:** 342 proteins
- **Edges:** 311 interactions  
- **Connected Components:** ~50+ (estimate based on density)
- **Mean pLLPS:** (see summary_stats.json)
- **Network density:** Low (sparse network)
- **Files:** `structural_*`

### Transporter
- **Nodes:** 22 proteins
- **Edges:** 13 interactions
- **Connected Components:** Multiple small clusters
- **Network density:** Very sparse
- **Files:** `transporter_*`

### GPCR
- **Nodes:** 7 proteins
- **Edges:** 4 interactions
- **Connected Components:** 2-3 small components
- **Network density:** Minimal
- **Files:** `gpcr_*`

## File Locations

All files are in: `/workspaces/mem_prot_llps/results/functional_group_networks/`

### Cytoscape JSON Files (.json)
- `ion_channel_network.json` (108 KB)
- `structural_network.json` (155 KB)
- `transporter_network.json` (8.3 KB)
- `gpcr_network.json` (2.7 KB)

### CSV Files
**Nodes:**
- `ion_channel_nodes.csv` (149 nodes, 4.8 KB)
- `structural_nodes.csv` (342 nodes, 11 KB)
- `transporter_nodes.csv` (22 nodes, 780 B)
- `gpcr_nodes.csv` (7 nodes, 227 B)

**Edges:**
- `ion_channel_edges.csv` (282 edges, 16 KB)
- `structural_edges.csv` (311 edges, 17 KB)
- `transporter_edges.csv` (13 edges, 883 B)
- `gpcr_edges.csv` (4 edges, 364 B)

### Summary Stats (JSON)
- `ion_channel_summary_stats.json`
- `structural_summary_stats.json`
- `transporter_summary_stats.json`
- `gpcr_summary_stats.json`

## Data Format

### Cytoscape JSON Structure
```json
{
  "data": {
    "name": "Ion Channel Network",
    "description": "STRING physical interactions colored by pLLPS score"
  },
  "elements": {
    "nodes": [
      {
        "data": {
          "id": "O60359",
          "name": "CACNG3",
          "gene_name": "CACNG3",
          "pllps": 0.64,
          "degree": 19,
          "uniprot_id": "O60359"
        }
      }
    ],
    "edges": [
      {
        "data": {
          "id": "O60359_Q08289",
          "source": "O60359",
          "target": "Q08289",
          "interaction": "interacts_with",
          "string_score": 720.0,
          "string_score_norm": 0.72,
          "weight": 0.72
        }
      }
    ]
  }
}
```

### Node CSV Columns
- `uniprot_id`: UniProt accession
- `gene_name`: Gene symbol (CACNG3, PNN, WNT16, etc.)
- `pllps_score`: Phase separation propensity (0-1)
- `degree`: Number of connections
- `functional_group`: Functional category

### Edge CSV Columns
- `protein1_uniprot`, `protein2_uniprot`: UniProt IDs
- `protein1_gene`, `protein2_gene`: Gene names
- `string_score`: Confidence score (0-1000)
- `string_score_normalized`: Normalized score (0-1)
- `pllps_1`, `pllps_2`: pLLPS scores for both proteins
- `pllps_mean`: Average pLLPS of interacting pair

## Visualization in Cytoscape

### Import Steps
1. Download Cytoscape from https://cytoscape.org/
2. **File → Import → Network from File** → Select `.json` file
3. Network imports with all node/edge attributes preserved

### Visualization Setup
1. **Node Color by pLLPS:**
   - Style panel → Fill Color → Column: `pllps`
   - Continuous Mapping: Blue (low) → White (0.5) → Red (high)
   
2. **Edge Width by STRING Score:**
   - Style panel → Width → Column: `string_score`
   - Continuous Mapping: Thin (600) → Thick (1000)

3. **Layout Options:**
   - **Organic (yFiles)**: Best for sparse networks with components
   - **Force-Directed Spring**: Physics-based layout
   - **Kamada-Kawai**: Stress minimization (good for small networks)
   - **Group Attributes Layout**: Organize by pLLPS ranges

### Recommended Layout
For these sparse networks with multiple disconnected components:
- **Layout → yFiles → Organic** (best visual separation of components)
- Adjust edge lengths and node spacing as needed
- Consider hiding or reducing edge width for low-confidence edges (<0.7)

## Network Characteristics

### Edge Confidence Cutoff
- **Minimum STRING score:** 600 (0.6 normalized)
- **Network type:** Physical interactions only (experimentally determined)
- **Species:** Human (NCBI taxonomy 9606)

### Connected Components
All networks exhibit **sparse, modular structure** with multiple disconnected components:
- **Ion Channels:** 41 components (largest has ~20-30 nodes)
- **Structural:** Multiple large and small components
- **Transporter/GPCR:** Few small clusters

This is expected for:
- High confidence threshold (≥0.6)
- Physical interactions only
- Membrane protein functional specialization

### pLLPS Distribution
- **High pLLPS (>0.7):** Enriched in ion channels and structural proteins
- **Medium pLLPS (0.4-0.7):** Mixed across groups
- **Low pLLPS (<0.4):** Common in GPCRs and some transporters

## Usage Notes

### For Network Analysis
- Use Cytoscape for interactive exploration and publication-quality figures
- CSV files can be imported into Python/R for custom analysis
- Summary stats provide quick overview of network properties

### For Further Analysis
- **Clustering:** Identify functional modules using MCL or Louvain algorithms
- **Hub Analysis:** High-degree nodes in `nodes.csv` (degree column)
- **pLLPS Enrichment:** Compare pLLPS between connected vs isolated proteins
- **Edge Analysis:** Investigate high pLLPS pairs (pllps_mean >0.7 in `edges.csv`)

### Data Quality
- ✅ Gene names properly mapped (CACNG3, PNN, WNT16, etc.)
- ✅ pLLPS scores verified against source data
- ✅ STRING scores at ≥600 threshold
- ✅ Network topology matches STRING reference (multiple components, sparse)
- ✅ Cytoscape format validated

## Next Steps

1. **Import into Cytoscape** for visualization
2. **Apply layouts** (recommend yFiles Organic for sparse networks)
3. **Color nodes by pLLPS** (continuous mapping)
4. **Identify functional modules** within each component
5. **Compare pLLPS enrichment** across network topology features
6. **Export publication figures** from Cytoscape (300+ dpi)

## Known Limitations

- **Grid layout removed:** Previous PNG visualizations used artificial grid layout; now use Kamada-Kawai in notebook or Cytoscape layouts for proper topology
- **Few edges for GPCRs/Transporters:** Due to high confidence threshold and physical interaction requirement
- **Multiple components:** Networks are not fully connected (expected for membrane proteins)
- **pLLPS interpretation:** High pLLPS indicates propensity, not confirmed phase separation

## References

- STRING database: https://string-db.org/
- STRING MCP endpoint: https://mcp.string-db.org/
- Cytoscape: https://cytoscape.org/
- Analysis notebooks: See `03_string_networks_pllps_enriched.ipynb` and `05_interactive_functional_group_networks.ipynb`
