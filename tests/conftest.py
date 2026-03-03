"""
Shared fixtures for llps test suite.
"""

import pytest
import pandas as pd
from pathlib import Path


@pytest.fixture()
def sample_pllps_df() -> pd.DataFrame:
    """Small synthetic pLLPS dataset."""
    return pd.DataFrame({
        'Entry': ['P001', 'P002', 'P003', 'P004', 'P005'],
        'Entry name': ['GENE1_HUMAN', 'GENE2_HUMAN', 'GENE3_HUMAN', 'GENE4_HUMAN', 'GENE5_HUMAN'],
        'Protein names': ['Protein 1', 'Protein 2', 'Protein 3', 'Protein 4', 'Protein 5'],
        'p(LLPS)': [0.9, 0.8, 0.5, 0.3, 0.1],
        'Subcellular location [CC]': [
            'SUBCELLULAR LOCATION: Nucleus. Cytoplasm.',
            'Nucleus; Cytosol.',
            '',
            'Mitochondrion.',
            None,
        ],
        'Length': [500, 400, 300, 200, 100],
    })


@pytest.fixture()
def sample_interactions_df() -> pd.DataFrame:
    """Small synthetic interaction dataset (STRING preferredName format)."""
    return pd.DataFrame({
        'preferredName_A': ['GENE1', 'GENE1', 'GENE2', 'GENE3'],
        'preferredName_B': ['GENE2', 'GENE3', 'GENE4', 'GENE4'],
        'score': [900, 800, 750, 700],
    })
