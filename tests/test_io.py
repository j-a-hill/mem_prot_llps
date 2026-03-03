"""
Tests for I/O and caching functions.
"""

import json
import pytest
import pandas as pd
from pathlib import Path

from llps.io import (
    save_interactions_to_cache,
    save_analysis_result,
    load_analysis_result,
    list_saved_results,
)


def test_save_interactions_to_cache(tmp_path: Path, sample_interactions_df: pd.DataFrame) -> None:
    cache_path = save_interactions_to_cache(
        sample_interactions_df,
        score_threshold=700,
        output_dir=str(tmp_path),
    )
    assert Path(cache_path).exists()
    with open(cache_path) as f:
        data = json.load(f)
    assert len(data) == len(sample_interactions_df)


def test_save_and_load_csv(tmp_path: Path, sample_pllps_df: pd.DataFrame) -> None:
    filepath = save_analysis_result(
        sample_pllps_df, "proteins", results_dir=str(tmp_path), format="csv"
    )
    loaded = load_analysis_result("proteins.csv", results_dir=str(tmp_path), format="csv")
    assert len(loaded) == len(sample_pllps_df)


def test_save_and_load_json(tmp_path: Path) -> None:
    data = {"key": "value", "count": 42}
    filepath = save_analysis_result(
        data, "test_data", results_dir=str(tmp_path), format="json"
    )
    loaded = load_analysis_result("test_data.json", results_dir=str(tmp_path), format="json")
    assert loaded["count"] == 42


def test_list_saved_results(tmp_path: Path, sample_pllps_df: pd.DataFrame) -> None:
    save_analysis_result(sample_pllps_df, "run1", results_dir=str(tmp_path), format="csv")
    save_analysis_result(sample_pllps_df, "run2", results_dir=str(tmp_path), format="csv")
    results = list_saved_results(results_dir=str(tmp_path))
    assert any("run1" in r for r in results)
    assert any("run2" in r for r in results)
