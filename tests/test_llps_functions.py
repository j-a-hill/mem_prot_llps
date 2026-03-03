"""
Tests for llps_functions.py

Tests focus on pure functions that do not require network access or large data files.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock

from llps_functions import (
    parse_location,
    add_location_columns,
    get_high_pllps_proteins,
    load_and_classify_data,
    match_interactions_to_pllps,
    analyze_interaction_enrichment,
    StringQueryConfig,
    STRING_API_BASE,
    STRING_API_NETWORK,
    STRING_API_GET_IDS,
    STRING_API_PARTNERS,
    _classify_network_nodes,
    _compute_network_metrics,
)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture()
def sample_pllps_df() -> pd.DataFrame:
    """Small synthetic pLLPS dataset."""
    return pd.DataFrame({
        'Entry': ['P001', 'P002', 'P003', 'P004', 'P005'],
        'Entry name': ['GENE1_HUMAN', 'GENE2_HUMAN', 'GENE3_HUMAN', 'GENE4_HUMAN', 'GENE5_HUMAN'],
        'Protein names': ['Protein 1', 'Protein 2', 'Protein 3', 'Protein 4', 'Protein 5'],
        'p(LLPS)': [0.9, 0.8, 0.5, 0.3, 0.1],
        'Subcellular location [CC]': [
            'SUBCELLULAR LOCATION: Nucleus. Cytoplasm.',
            'Nucleus; Cytosol.',
            '',
            'Mitochondrion.',
            None,
        ],
        'Length': [500, 400, 300, 200, 100],
    })


@pytest.fixture()
def sample_interactions_df() -> pd.DataFrame:
    """Small synthetic interaction dataset (STRING preferredName format)."""
    return pd.DataFrame({
        'preferredName_A': ['GENE1', 'GENE1', 'GENE2', 'GENE3'],
        'preferredName_B': ['GENE2', 'GENE3', 'GENE4', 'GENE4'],
        'score': [900, 800, 750, 700],
    })


# =============================================================================
# STRING API Constants
# =============================================================================

def test_string_api_constants() -> None:
    assert STRING_API_BASE == "https://string-db.org/api/json"
    assert STRING_API_NETWORK == f"{STRING_API_BASE}/network"
    assert STRING_API_GET_IDS == f"{STRING_API_BASE}/get_string_ids"
    assert STRING_API_PARTNERS == f"{STRING_API_BASE}/interaction_partners"


# =============================================================================
# StringQueryConfig
# =============================================================================

def test_string_query_config_defaults() -> None:
    cfg = StringQueryConfig()
    assert cfg.species == 9606
    assert cfg.score_threshold == 700
    assert cfg.batch_size == 100
    assert cfg.use_cache is True
    assert cfg.network_type == "physical"


def test_string_query_config_custom() -> None:
    cfg = StringQueryConfig(species=10090, score_threshold=400, use_cache=False)
    assert cfg.species == 10090
    assert cfg.score_threshold == 400
    assert cfg.use_cache is False
    assert cfg.batch_size == 100  # default unchanged


# =============================================================================
# parse_location
# =============================================================================

def test_parse_location_basic() -> None:
    result = parse_location("SUBCELLULAR LOCATION: Nucleus. Cytoplasm.")
    assert "Nucleus" in result
    assert "Cytoplasm" in result


def test_parse_location_semicolons() -> None:
    result = parse_location("Nucleus; Cytosol")
    assert "Nucleus" in result
    assert "Cytosol" in result


def test_parse_location_empty_string() -> None:
    assert parse_location('') == []


def test_parse_location_nan() -> None:
    assert parse_location(float('nan')) == []


def test_parse_location_none() -> None:
    assert parse_location(None) == []


def test_parse_location_removes_evidence_codes() -> None:
    result = parse_location("Nucleus {ECO:0000269|PubMed:12345}.")
    assert all("{" not in r for r in result)
    assert "Nucleus" in result


def test_parse_location_removes_isoform_tags() -> None:
    result = parse_location("Isoform 1: Nucleus. Cytoplasm.")
    assert "Nucleus" in result or "Cytoplasm" in result


def test_parse_location_case_normalised() -> None:
    result = parse_location("cytosol")
    assert "Cytosol" in result


def test_parse_location_no_duplicates() -> None:
    result = parse_location("Nucleus; Nucleus; Cytoplasm")
    assert result.count("Nucleus") == 1


# =============================================================================
# add_location_columns
# =============================================================================

def test_add_location_columns_adds_column(sample_pllps_df: pd.DataFrame) -> None:
    result = add_location_columns(sample_pllps_df)
    assert 'Location Categories' in result.columns


def test_add_location_columns_is_list(sample_pllps_df: pd.DataFrame) -> None:
    result = add_location_columns(sample_pllps_df)
    assert all(isinstance(v, list) for v in result['Location Categories'])


def test_add_location_columns_no_mutation(sample_pllps_df: pd.DataFrame) -> None:
    """Original dataframe should not be modified."""
    original_cols = set(sample_pllps_df.columns)
    add_location_columns(sample_pllps_df)
    assert set(sample_pllps_df.columns) == original_cols


def test_add_location_columns_missing_col() -> None:
    """DataFrame without the location column should be returned unchanged."""
    df = pd.DataFrame({'Entry': ['P001'], 'p(LLPS)': [0.9]})
    result = add_location_columns(df)
    assert 'Location Categories' not in result.columns


# =============================================================================
# get_high_pllps_proteins
# =============================================================================

def test_get_high_pllps_proteins_absolute(sample_pllps_df: pd.DataFrame) -> None:
    df_out, high_ids = get_high_pllps_proteins(sample_pllps_df, threshold=0.7)
    assert 'P001' in high_ids
    assert 'P002' in high_ids
    assert 'P005' not in high_ids


def test_get_high_pllps_proteins_adds_class_column(sample_pllps_df: pd.DataFrame) -> None:
    df_out, _ = get_high_pllps_proteins(sample_pllps_df, threshold=0.7)
    assert 'pLLPS_class' in df_out.columns


def test_get_high_pllps_proteins_percentile(sample_pllps_df: pd.DataFrame) -> None:
    df_out, high_ids = get_high_pllps_proteins(sample_pllps_df, threshold=80, method='percentile')
    assert len(high_ids) > 0


def test_get_high_pllps_proteins_missing_column() -> None:
    df = pd.DataFrame({'x': [1, 2]})
    with pytest.raises(ValueError):
        get_high_pllps_proteins(df)


# =============================================================================
# load_and_classify_data
# =============================================================================

def test_load_and_classify_data(sample_pllps_df: pd.DataFrame, tmp_path: Path) -> None:
    filepath = tmp_path / "test_data.xlsx"
    sample_pllps_df.to_excel(filepath, index=False)
    result = load_and_classify_data(str(filepath))
    assert 'pLLPS_class' in result.columns
    assert set(result['pLLPS_class'].unique()).issubset({'High', 'Medium', 'Low'})


def test_load_and_classify_data_thresholds(sample_pllps_df: pd.DataFrame, tmp_path: Path) -> None:
    filepath = tmp_path / "test_data.xlsx"
    sample_pllps_df.to_excel(filepath, index=False)
    result = load_and_classify_data(str(filepath), high_threshold=0.8, low_threshold=0.3)
    high_count = (result['pLLPS_class'] == 'High').sum()
    assert high_count == 2  # P001 (0.9) and P002 (0.8) both >= 0.8


# =============================================================================
# match_interactions_to_pllps
# =============================================================================

def test_match_interactions_empty(sample_pllps_df: pd.DataFrame) -> None:
    empty = pd.DataFrame()
    result = match_interactions_to_pllps(empty, sample_pllps_df)
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 0


def test_match_interactions_basic(
    sample_interactions_df: pd.DataFrame, sample_pllps_df: pd.DataFrame
) -> None:
    result = match_interactions_to_pllps(sample_interactions_df, sample_pllps_df)
    assert isinstance(result, pd.DataFrame)


# =============================================================================
# analyze_interaction_enrichment
# =============================================================================

def test_analyze_enrichment_empty() -> None:
    assert analyze_interaction_enrichment(pd.DataFrame()) is None


def test_analyze_enrichment_basic(sample_pllps_df: pd.DataFrame) -> None:
    matched_df = pd.DataFrame({
        'pllps_1': [0.9, 0.8, 0.5, 0.3, 0.1, 0.9],
        'pllps_2': [0.8, 0.5, 0.3, 0.1, 0.9, 0.3],
        'protein1': ['P001', 'P001', 'P002', 'P003', 'P005', 'P001'],
        'protein2': ['P002', 'P003', 'P004', 'P004', 'P001', 'P004'],
    })
    result = analyze_interaction_enrichment(matched_df, threshold=0.7)
    if result is not None:
        assert 'enrichment' in result
        assert 'total' in result
        assert result['total'] > 0


# =============================================================================
# _classify_network_nodes (requires networkx)
# =============================================================================

def test_classify_network_nodes() -> None:
    pytest.importorskip("networkx")
    import networkx as nx
    G = nx.Graph()
    G.add_node('A', pLLPS=0.9)
    G.add_node('B', pLLPS=0.8)
    G.add_node('C', pLLPS=0.3)
    G.add_node('D', pLLPS=None)
    high, low, unknown = _classify_network_nodes(G, high_threshold=0.7)
    assert 'A' in high
    assert 'B' in high
    assert 'C' in low
    assert 'D' in unknown


# =============================================================================
# _compute_network_metrics (requires networkx)
# =============================================================================

def test_compute_network_metrics_empty_graph() -> None:
    pytest.importorskip("networkx")
    import networkx as nx
    G = nx.Graph()
    results = _compute_network_metrics(G, [], [], [])
    assert results['total_nodes'] == 0
    assert results['total_edges'] == 0
    assert results['enrichment_ratio'] == 0.0


def test_compute_network_metrics_simple() -> None:
    pytest.importorskip("networkx")
    import networkx as nx
    G = nx.Graph()
    G.add_node('A', pLLPS=0.9)
    G.add_node('B', pLLPS=0.8)
    G.add_node('C', pLLPS=0.3)
    G.add_edge('A', 'B')
    G.add_edge('A', 'C')
    results = _compute_network_metrics(G, ['A', 'B'], ['C'], [])
    assert results['total_nodes'] == 3
    assert results['total_edges'] == 2
    assert results['high_high_interactions'] == 1
    assert results['high_low_interactions'] == 1
    assert results['low_low_interactions'] == 0
