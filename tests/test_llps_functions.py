"""
Backward-compatibility import test for llps_functions module.

Verifies that all public symbols are accessible via the legacy import path.
"""

import pytest


def test_backward_compat_imports() -> None:
    """All public symbols remain importable from the top-level module."""
    from llps_functions import (
        STRING_API_BASE,
        STRING_API_NETWORK,
        STRING_API_GET_IDS,
        STRING_API_PARTNERS,
        StringQueryConfig,
        load_llps_data,
        load_and_classify_data,
        get_high_pllps_proteins,
        parse_location,
        add_location_columns,
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
        save_analysis_result,
        load_analysis_result,
        list_saved_results,
        parse_function_categories,
        is_membrane_protein,
        classify_protein_function,
        filter_membrane_proteins,
        add_functional_categories,
    )
    assert STRING_API_BASE == "https://string-db.org/api/json"
    assert callable(load_llps_data)
    assert callable(fetch_string_interactions)
    assert callable(analyze_network)
    assert callable(analyze_interaction_enrichment)
