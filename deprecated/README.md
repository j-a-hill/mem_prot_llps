# Deprecated Notebooks and Documentation

This folder contains notebooks and documentation from earlier versions of the analysis workflow. These files have been superseded by the current workflow but are retained for reference.

## Deprecated Notebooks

### Original Analysis Workflow (Superseded 2025-01-14)

- `02_string_interactions.ipynb` - Original STRING interaction fetching (replaced by `03_string_networks_pllps_enriched.ipynb`)
- `03_enrichment_analysis.ipynb` - Original enrichment analysis (replaced by `02_pllps_enriched_functional_groups.ipynb`)
- `04_network_analysis.ipynb` - Original network analysis (replaced by `04_visualize_pllps_networks.ipynb`)
- `05_functional_groups.ipynb` - Original functional grouping (integrated into new workflow)
- `06_visualization_summary.ipynb` - Original visualization summary (replaced by new notebooks)

### Exploration and Other Deprecated Notebooks

- `pllps_interaction_analysis_DEPRECATED.ipynb` - Early exploratory analysis
- `exploration_notebook.ipynb` - General exploration notebook
- `kegg_pathway_analysis.ipynb` - KEGG pathway analysis (see `KEGG_PATHWAY_NOTEBOOK_SUMMARY.md` in root)

## Deprecated Documentation

- `DEPRECATED_NOTEBOOK_README.md` - Original deprecation notice
- `MIGRATION_GUIDE.md` - Migration guide for old workflow
- `REFACTORING_SUMMARY.md` - Summary of refactoring changes
- `TIDY_UP_SUMMARY.txt` - Previous tidy-up notes

## Current Workflow

Please refer to the root directory for the current analysis workflow:

1. `01_data_loading_and_classification.ipynb` - Load and classify pLLPS data
2. `02_pllps_enriched_functional_groups.ipynb` - Identify functionally enriched groups
3. `03_string_networks_pllps_enriched.ipynb` - Fetch STRING interactions for enriched groups
4. `04_visualize_pllps_networks.ipynb` - Visualize pLLPS-colored networks

See `README.md` and `NEW_WORKFLOW_COMPLETION.md` in the root directory for full documentation.

---

**Date Deprecated**: January 14, 2025  
**Reason**: Workflow redesigned to focus on functionally enriched groups with pLLPS-colored STRING network visualizations
