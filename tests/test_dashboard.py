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
from app import (  # noqa: E402
    _parse_location,
    _parse_functions,
    _enrich_df,
    _all_values,
    _chart_to_div,
    _df_to_csv_url,
)


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


# ---------------------------------------------------------------------------
# _chart_to_div
# ---------------------------------------------------------------------------

def test_chart_to_div_structure() -> None:
    import altair as alt
    df = pd.DataFrame({"x": [1, 2, 3], "y": [4, 5, 6]})
    chart = alt.Chart(df).mark_bar().encode(x="x:O", y="y:Q")
    html = _chart_to_div(chart, "test_div")
    assert 'id="test_div"' in html
    assert "vegaEmbed" in html
    assert "<script>" in html
    assert "cdn.plot.ly" not in html


def test_chart_to_div_unique_ids() -> None:
    import altair as alt
    df = pd.DataFrame({"x": [1], "y": [2]})
    chart = alt.Chart(df).mark_point().encode(x="x:Q", y="y:Q")
    html_a = _chart_to_div(chart, "chart_a")
    html_b = _chart_to_div(chart, "chart_b")
    assert "#chart_a" in html_a
    assert "#chart_b" in html_b
    assert "#chart_a" not in html_b


# ---------------------------------------------------------------------------
# _df_to_csv_url
# ---------------------------------------------------------------------------

def test_df_to_csv_url_format() -> None:
    df = pd.DataFrame({"a": [1.0, 2.0], "b": ["X", "Y"]})
    source = _df_to_csv_url(df)
    assert source.url.startswith("data:text/csv;base64,")


def test_df_to_csv_url_roundtrip() -> None:
    import base64
    df = pd.DataFrame({"a": [1.0, 2.0], "b": ["X", "Y"]})
    source = _df_to_csv_url(df)
    decoded = base64.b64decode(source.url.split(",", 1)[1]).decode()
    assert "a,b" in decoded
    assert "1.0" in decoded
    assert "X" in decoded
