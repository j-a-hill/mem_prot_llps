"""
Tests for the dashboard helper functions (dashboard/app.py).

These tests cover the pure-Python logic that can be exercised without a
running Shiny server: data enrichment, location parsing, and function
classification.
"""

import sys
from pathlib import Path
import pandas as pd
import pytest

# Make dashboard/ importable without installing it as a package.
sys.path.insert(0, str(Path(__file__).parent.parent / "dashboard"))
from app import _parse_location, _parse_functions, _enrich_df, _all_values  # noqa: E402


# ---------------------------------------------------------------------------
# _parse_location
# ---------------------------------------------------------------------------

def test_parse_location_nucleus() -> None:
    result = _parse_location("SUBCELLULAR LOCATION: Nucleus. Cytoplasm.")
    assert "Nucleus" in result
    assert "Cytoplasm" in result


def test_parse_location_removes_evidence_codes() -> None:
    result = _parse_location("SUBCELLULAR LOCATION: Nucleus {ECO:0000255}.")
    assert any("Nucleus" in v for v in result)
    assert not any("ECO" in v for v in result)


def test_parse_location_empty() -> None:
    assert _parse_location("") == []
    assert _parse_location(None) == []


def test_parse_location_deduplicates() -> None:
    result = _parse_location("Nucleus. Nucleus.")
    assert result.count("Nucleus") == 1


# ---------------------------------------------------------------------------
# _parse_functions
# ---------------------------------------------------------------------------

def test_parse_functions_kinase() -> None:
    result = _parse_functions("Serine/threonine-protein kinase activity")
    assert "Kinase" in result


def test_parse_functions_receptor_from_name() -> None:
    result = _parse_functions(None, "G-protein coupled receptor 5")
    assert "Receptor" in result


def test_parse_functions_no_match() -> None:
    result = _parse_functions("Unknown function", "Unknown protein")
    assert result == []


def test_parse_functions_multiple_categories() -> None:
    result = _parse_functions("This protein is a kinase and a transporter")
    assert "Kinase" in result
    assert "Transporter" in result


# ---------------------------------------------------------------------------
# _enrich_df
# ---------------------------------------------------------------------------

def test_enrich_df_adds_location_categories() -> None:
    df = pd.DataFrame({
        "Entry": ["P001"],
        "p(LLPS)": [0.8],
        "Subcellular location [CC]": ["Nucleus."],
        "Length": [400],
    })
    result = _enrich_df(df)
    assert "Location Categories" in result.columns
    assert isinstance(result["Location Categories"].iloc[0], list)


def test_enrich_df_adds_function_categories() -> None:
    df = pd.DataFrame({
        "Entry": ["P001"],
        "p(LLPS)": [0.8],
        "Function [CC]": ["kinase activity"],
        "Protein names": ["Kinase 1"],
    })
    result = _enrich_df(df)
    assert "Function Categories" in result.columns


def test_enrich_df_pllps_class_high() -> None:
    df = pd.DataFrame({"Entry": ["A"], "p(LLPS)": [0.9]})
    result = _enrich_df(df)
    assert str(result["pLLPS_class"].iloc[0]) == "High"


def test_enrich_df_pllps_class_low() -> None:
    df = pd.DataFrame({"Entry": ["A"], "p(LLPS)": [0.2]})
    result = _enrich_df(df)
    assert str(result["pLLPS_class"].iloc[0]) == "Low"


def test_enrich_df_no_location_column() -> None:
    df = pd.DataFrame({"Entry": ["A"], "p(LLPS)": [0.5]})
    result = _enrich_df(df)
    assert "Location Categories" in result.columns
    assert result["Location Categories"].iloc[0] == []


# ---------------------------------------------------------------------------
# _all_values
# ---------------------------------------------------------------------------

def test_all_values_returns_sorted_unique() -> None:
    df = pd.DataFrame({
        "cats": [["Kinase", "Receptor"], ["Kinase"], ["Transporter"]],
    })
    result = _all_values(df, "cats")
    assert result == sorted({"Kinase", "Receptor", "Transporter"})


def test_all_values_skips_empty_strings() -> None:
    df = pd.DataFrame({"cats": [["Kinase", ""], [""]]})
    result = _all_values(df, "cats")
    assert "" not in result
