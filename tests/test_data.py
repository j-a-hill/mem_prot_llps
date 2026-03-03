"""
Tests for data loading and classification functions.
"""

import pytest
import pandas as pd
from pathlib import Path

from llps_functions import get_high_pllps_proteins, load_and_classify_data


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
