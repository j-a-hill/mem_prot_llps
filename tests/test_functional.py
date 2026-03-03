"""
Tests for functional classification utilities.
"""

import pytest
import pandas as pd

from llps.functional import (
    parse_function_categories,
    is_membrane_protein,
    classify_protein_function,
    add_functional_categories,
)


def test_is_membrane_protein_positive() -> None:
    names = "Membrane glycoprotein; Transmembrane domain"
    assert is_membrane_protein(names) is True


def test_is_membrane_protein_negative() -> None:
    names = "Cytoplasmic kinase"
    assert is_membrane_protein(names) is False


def test_is_membrane_protein_empty() -> None:
    assert is_membrane_protein("") is False


def test_is_membrane_protein_nan() -> None:
    assert is_membrane_protein(float("nan")) is False


def test_parse_function_categories_returns_list() -> None:
    result = parse_function_categories("DNA binding; transcription factor; Zinc finger")
    assert isinstance(result, list)


def test_parse_function_categories_empty() -> None:
    assert parse_function_categories("") == []


def test_classify_protein_function_returns_list() -> None:
    result = classify_protein_function("DNA-binding transcription factor")
    assert isinstance(result, list)


def test_add_functional_categories_adds_column(sample_pllps_df: pd.DataFrame) -> None:
    df_with_func = sample_pllps_df.copy()
    df_with_func['Function [CC]'] = [
        'FUNCTION: DNA binding.',
        'FUNCTION: Kinase activity.',
        '',
        None,
        'FUNCTION: Membrane protein.',
    ]
    result = add_functional_categories(df_with_func)
    assert 'Functional_Categories' in result.columns


def test_add_functional_categories_list_values(sample_pllps_df: pd.DataFrame) -> None:
    df_with_func = sample_pllps_df.copy()
    df_with_func['Function [CC]'] = [
        'FUNCTION: DNA binding.',
        'FUNCTION: Kinase activity.',
        '',
        None,
        'FUNCTION: Membrane protein.',
    ]
    result = add_functional_categories(df_with_func)
    assert all(isinstance(v, list) for v in result['Functional_Categories'])
