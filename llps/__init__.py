"""
llps package - LLPS Analysis Functions

This package provides all functions for analysing Liquid-Liquid Phase Separation
(LLPS) protein data and their interaction networks.

All public symbols are re-exported here so that existing code using
``from llps_functions import X`` continues to work after the module is
refactored to import from this package.
"""

from llps.constants import (
    STRING_API_BASE,
    STRING_API_GET_IDS,
    STRING_API_NETWORK,
    STRING_API_PARTNERS,
    StringQueryConfig,
)

from llps.data import (
    load_llps_data,
    load_and_classify_data,
    get_high_pllps_proteins,
    fetch_uniprot_tm_annotations,
    add_tmd_count,
    fetch_uniprot_go_annotations,
    add_go_annotations,
    fetch_uniprot_location_sl_ids,
)

from llps.location import (
    parse_location,
    parse_sl_ids,
    compartment_from_sl_ids,
    add_location_columns,
    categorize_location_to_compartment,
    load_subcell_ontology,
    download_subcell_terms,
    is_membrane_localized,
)

from llps.string_api import (
    get_string_mapping,
    fetch_string_interactions,
    fetch_interaction_partners,
    load_string_network_file,
)

from llps.network import (
    match_interactions_to_pllps,
    match_interactors_to_pllps,
    analyze_network,
)

from llps.enrichment import (
    analyze_interaction_enrichment,
    analyze_interaction_matrix,
    analyze_interactions_by_location,
)

from llps.visualization import (
    plot_interaction_heatmap,
    plot_location_heatmaps,
    print_analysis_report,
)

from llps.io import (
    save_interactions_to_cache,
    export_protein_list,
    get_string_interactions,
    save_analysis_result,
    load_analysis_result,
    list_saved_results,
)

from llps.functional import (
    parse_function_categories,
    is_membrane_protein,
    classify_protein_function,
    filter_membrane_proteins,
    add_functional_categories,
    count_tm_domains,
    parse_go_ids,
    add_go_slim_categories,
    map_go_ids_to_slim,
)

__all__ = [
    # constants
    "STRING_API_BASE",
    "STRING_API_GET_IDS",
    "STRING_API_NETWORK",
    "STRING_API_PARTNERS",
    "StringQueryConfig",
    # data
    "load_llps_data",
    "load_and_classify_data",
    "get_high_pllps_proteins",
    "fetch_uniprot_tm_annotations",
    "add_tmd_count",
    "fetch_uniprot_go_annotations",
    "add_go_annotations",
    # location
    "parse_location",
    "add_location_columns",
    "categorize_location_to_compartment",
    "load_subcell_ontology",
    "download_subcell_terms",
    "is_membrane_localized",
    "analyze_interactions_by_location",
    # string_api
    "get_string_mapping",
    "fetch_string_interactions",
    "fetch_interaction_partners",
    "load_string_network_file",
    # network
    "match_interactions_to_pllps",
    "match_interactors_to_pllps",
    "analyze_network",
    # enrichment
    "analyze_interaction_enrichment",
    "analyze_interaction_matrix",
    # visualization
    "plot_interaction_heatmap",
    "plot_location_heatmaps",
    "print_analysis_report",
    # io
    "save_interactions_to_cache",
    "export_protein_list",
    "get_string_interactions",
    "save_analysis_result",
    "load_analysis_result",
    "list_saved_results",
    # functional
    "parse_function_categories",
    "is_membrane_protein",
    "classify_protein_function",
    "filter_membrane_proteins",
    "add_functional_categories",
    "count_tm_domains",
    "parse_go_ids",
    "add_go_slim_categories",
    "map_go_ids_to_slim",
]
