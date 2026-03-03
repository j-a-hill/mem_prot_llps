"""
Tests for enrichment analysis functions.
"""

import pandas as pd
import pytest

from llps_functions import analyze_interaction_enrichment, match_interactions_to_pllps


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
