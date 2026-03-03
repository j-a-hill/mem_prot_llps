"""
STRING API constants and query configuration dataclass.
"""

from dataclasses import dataclass


# =============================================================================
# STRING API CONSTANTS
# =============================================================================

STRING_API_BASE = "https://string-db.org/api/json"
STRING_API_GET_IDS = f"{STRING_API_BASE}/get_string_ids"
STRING_API_NETWORK = f"{STRING_API_BASE}/network"
STRING_API_PARTNERS = f"{STRING_API_BASE}/interaction_partners"


@dataclass
class StringQueryConfig:
    """Configuration for STRING database API queries."""

    species: int = 9606
    score_threshold: int = 700
    batch_size: int = 100
    use_cache: bool = True
    network_type: str = "physical"
