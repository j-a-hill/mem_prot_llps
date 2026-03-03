"""
Tests for STRING API constants and StringQueryConfig.
"""

from llps_functions import (
    STRING_API_BASE,
    STRING_API_NETWORK,
    STRING_API_GET_IDS,
    STRING_API_PARTNERS,
    StringQueryConfig,
)


def test_string_api_constants() -> None:
    assert STRING_API_BASE == "https://string-db.org/api/json"
    assert STRING_API_NETWORK == f"{STRING_API_BASE}/network"
    assert STRING_API_GET_IDS == f"{STRING_API_BASE}/get_string_ids"
    assert STRING_API_PARTNERS == f"{STRING_API_BASE}/interaction_partners"


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
