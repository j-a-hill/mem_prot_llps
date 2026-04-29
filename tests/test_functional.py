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
    # tmd_count is the authoritative signal after the GO/ontology migration
    assert is_membrane_protein("", tmd_count=2) is True


def test_is_membrane_protein_negative() -> None:
    names = "Cytoplasmic kinase"
    assert is_membrane_protein(names) is False


def test_is_membrane_protein_empty() -> None:
    assert is_membrane_protein("") is False


def test_is_membrane_protein_nan() -> None:
    assert is_membrane_protein(float("nan")) is False


def test_is_membrane_protein_ignore_function_by_default() -> None:
    # By default function text is ignored: a function-only signal should not mark membrane
    names = "Membrane glycoprotein; Transmembrane domain"
    assert is_membrane_protein(names) is False


def test_is_membrane_protein_weak_membrane_word_without_tmd() -> None:
    # Example like PCLO: function mentions 'membrane fusion' but there is no TM domain
    func = "Scaffold protein involved in membrane fusion and vesicle transport"
    loc = "Cell junction, synapse, presynaptic active zone"
    # No TMD info provided -> should not be classified as membrane
    assert is_membrane_protein(func, protein_name_str="PCLO", location_str=loc, tmd_count=0) is False


def test_parse_function_categories_returns_list() -> None:
    result = parse_function_categories("DNA binding; transcription factor; Zinc finger")
    assert isinstance(result, list)


def test_parse_function_categories_empty() -> None:
    assert parse_function_categories("") == []


def test_classify_protein_function_returns_list() -> None:
    result = classify_protein_function("DNA-binding transcription factor")
    assert isinstance(result, list)


def test_add_functional_categories_adds_column(sample_pllps_df: pd.DataFrame) -> None:
    df = sample_pllps_df.copy()
    # GO:0003700 = transcription factor activity; GO:0016301 = kinase activity
    df["GO_IDs"] = [
        "GO:0003700",
        "GO:0016301",
        "",
        None,
        "GO:0005216",
    ]
    result = add_functional_categories(df)
    assert "Functional_Categories" in result.columns


def test_add_functional_categories_list_values(sample_pllps_df: pd.DataFrame) -> None:
    df = sample_pllps_df.copy()
    df["GO_IDs"] = [
        "GO:0003700",
        "GO:0016301",
        "",
        None,
        "GO:0005216",
    ]
    result = add_functional_categories(df)
    assert all(isinstance(v, list) for v in result["Functional_Categories"])
