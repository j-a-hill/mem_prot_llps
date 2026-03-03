"""
Tests for location parsing functions.
"""

import pytest
import pandas as pd

from llps_functions import parse_location, add_location_columns


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
