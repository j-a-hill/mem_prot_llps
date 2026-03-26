"""
LLPS Analysis Master Functions Module

This module consolidates all functions for analyzing Liquid-Liquid Phase Separation (LLPS) 
protein data and their interaction networks. All analysis notebooks and scripts should 
import from this master file.

USAGE:
------
In Python scripts:
    from llps_functions import load_llps_data, fetch_string_interactions
    
In Jupyter notebooks:
    import llps_functions as lf
    df = lf.load_llps_data('path/to/data.xlsx')

Organized sections:
1. Data Loading and Classification
2. Location Parsing and Analysis
3. STRING Database Interactions
4. Network Analysis
5. Enrichment Analysis
6. Visualization Functions
7. Export and Caching Functions
8. Result Saving and Loading Utilities
9. Functional Classification Utilities

Author: LLPS Analysis Team
Date: 2025
"""

# Re-export everything from the llps sub-package for backward compatibility.
from llps import (
    STRING_API_BASE,
    STRING_API_GET_IDS,
    STRING_API_NETWORK,
    STRING_API_PARTNERS,
    StringQueryConfig,
    load_llps_data,
    load_and_classify_data,
    get_high_pllps_proteins,
    parse_location,
    add_location_columns,
    categorize_location_to_compartment,
    analyze_interactions_by_location,
    get_string_mapping,
    fetch_string_interactions,
    fetch_interaction_partners,
    load_string_network_file,
    match_interactions_to_pllps,
    match_interactors_to_pllps,
    analyze_network,
    analyze_interaction_enrichment,
    analyze_interaction_matrix,
    plot_interaction_heatmap,
    plot_location_heatmaps,
    print_analysis_report,
    save_interactions_to_cache,
    export_protein_list,
    get_string_interactions,
    save_analysis_result,
    load_analysis_result,
    list_saved_results,
    parse_function_categories,
    is_membrane_protein,
    classify_protein_function,
    filter_membrane_proteins,
    add_functional_categories,
)  # noqa: F401
