"""
STRING Protein Interaction Functions

⚠️ DEPRECATION NOTICE:
This module has been consolidated into llps_functions.py.
Please update your imports to use:
    from llps_functions import fetch_string_interactions, ...
    
This file is kept for backward compatibility but may be removed in future versions.

Reusable functions for fetching and analyzing protein interactions from STRING database.
These functions can be used in both the Shiny app and Jupyter notebooks.

Author: LLPS Analysis Team
"""

import warnings
warnings.warn(
    "string_functions.py is deprecated. Please use llps_functions.py instead.",
    DeprecationWarning,
    stacklevel=2
)

import pandas as pd
import numpy as np
from pathlib import Path
import time
import requests
from typing import List, Dict, Tuple, Optional, Callable
from scipy.stats import chi2_contingency


def fetch_string_interactions(
    protein_ids: List[str],
    species: int = 9606,
    score_threshold: int = 700,
    batch_size: int = 100,
    progress_callback: Optional[Callable[[str], None]] = None,
    use_cache: bool = True
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Fetch protein-protein interactions from STRING database with optional caching.
    
    Parameters
    ----------
    protein_ids : List[str]
        List of UniProt protein IDs to fetch interactions for
    species : int, optional
        NCBI taxonomy ID (default: 9606 for human)
    score_threshold : int, optional
        Minimum combined score 0-1000 (default: 700 for high confidence)
        400=medium, 700=high, 900=highest confidence
    batch_size : int, optional
        Number of proteins per API request (default: 100)
    progress_callback : Callable, optional
        Function to call with progress updates
    use_cache : bool, optional
        Whether to try loading from cache first (default: True)
    
    Returns
    -------
    Tuple[pd.DataFrame, List[str]]
        DataFrame with interaction data and list of error messages
    
    Examples
    --------
    >>> protein_ids = ['P04637', 'P38398', 'P51587']  # p53, BRCA1, BRCA2
    >>> interactions_df, errors = fetch_string_interactions(protein_ids, score_threshold=700)
    >>> print(f"Found {len(interactions_df)} interactions")
    """
    string_api_url = "https://string-db.org/api/json/network"
    all_interactions = []
    errors = []
    
    # Try to load from cache first
    if use_cache:
        cache_file = Path(__file__).parent / "data" / f"string_cache_{score_threshold}.json"
        if cache_file.exists():
            try:
                import json
                with open(cache_file, 'r') as f:
                    cached_data = json.load(f)
                    # Filter cached data for requested proteins
                    all_interactions = [
                        i for i in cached_data 
                        if any(pid in str(i.get('preferredName_A', '')) + str(i.get('preferredName_B', ''))
                              for pid in protein_ids[:batch_size * 2])  # Check first few batches
                    ]
                    if all_interactions:
                        if progress_callback:
                            progress_callback(f"Loaded {len(all_interactions)} interactions from cache")
                        return pd.DataFrame(all_interactions), ["Using cached data"]
            except Exception as e:
                errors.append(f"Cache load failed: {str(e)}")
    
    total_batches = (len(protein_ids) + batch_size - 1) // batch_size
    
    # Fetch interactions from STRING API
    for i in range(0, len(protein_ids), batch_size):
        batch = protein_ids[i:i+batch_size]
        batch_num = i // batch_size + 1
        
        if progress_callback:
            progress_callback(f"Fetching batch {batch_num}/{total_batches}...")
        
        params = {
            "identifiers": "\r".join(batch),
            "species": species,
            "required_score": score_threshold,
            "network_type": "physical",  # Physical interactions only
            "caller_identity": "pllps_analysis"
        }
        
        try:
            response = requests.post(string_api_url, data=params, timeout=60)
            if response.status_code == 200:
                interactions = response.json()
                all_interactions.extend(interactions)
            elif response.status_code == 429:
                errors.append(f"Rate limited on batch {batch_num}, retrying...")
                time.sleep(30)
                response = requests.post(string_api_url, data=params, timeout=60)
                if response.status_code == 200:
                    interactions = response.json()
                    all_interactions.extend(interactions)
                else:
                    errors.append(f"Failed after retry: {response.status_code}")
            else:
                errors.append(f"Batch {batch_num}: HTTP {response.status_code}")
        except requests.Timeout:
            errors.append(f"Batch {batch_num}: Timeout")
        except requests.RequestException as e:
            errors.append(f"Batch {batch_num}: Network error - {str(e)}")
            # If network error, suggest using cached data
            if "Failed to resolve" in str(e) or "No address" in str(e):
                errors.append("Network unavailable. To use cached data:")
                errors.append("1. Run STRING analysis locally with network access")
                errors.append(f"2. Save results to data/string_cache_{score_threshold}.json")
                errors.append("3. Upload the cache file with your data")
        except ValueError as e:
            errors.append(f"Batch {batch_num}: JSON decode error - {str(e)}")
        
        time.sleep(1)  # Rate limiting
    
    if progress_callback:
        progress_callback(f"Completed: Retrieved {len(all_interactions)} interactions")
    
    return pd.DataFrame(all_interactions) if all_interactions else pd.DataFrame(), errors


def match_interactions_to_pllps(
    interactions_df: pd.DataFrame,
    pllps_df: pd.DataFrame,
    id_column: str = 'Entry'
) -> pd.DataFrame:
    """
    Match STRING interactions to pLLPS dataset and add pLLPS scores.
    
    Parameters
    ----------
    interactions_df : pd.DataFrame
        DataFrame from STRING with 'preferredName_A' and 'preferredName_B' columns
    pllps_df : pd.DataFrame
        DataFrame with protein data including pLLPS scores
    id_column : str, optional
        Column name in pllps_df containing protein IDs (default: 'Entry')
    
    Returns
    -------
    pd.DataFrame
        Interactions with added pLLPS scores for both partners
    
    Examples
    --------
    >>> matched_df = match_interactions_to_pllps(interactions_df, pllps_df)
    >>> print(f"Matched {len(matched_df)} interactions to pLLPS data")
    """
    if len(interactions_df) == 0:
        return pd.DataFrame()
    
    # Create lookup dictionary for pLLPS scores
    pllps_dict = pllps_df.set_index(id_column)['p(LLPS)'].to_dict()
    
    # Match interactions to pLLPS data
    matched = interactions_df.copy()
    matched['pllps_a'] = matched['preferredName_A'].map(pllps_dict)
    matched['pllps_b'] = matched['preferredName_B'].map(pllps_dict)
    
    # Keep only interactions where both partners are in the dataset
    matched = matched.dropna(subset=['pllps_a', 'pllps_b'])
    
    return matched


def analyze_interaction_enrichment(
    matched_df: pd.DataFrame,
    threshold: float = 0.7
) -> Optional[Dict]:
    """
    Analyze enrichment of interactions between high vs low pLLPS proteins.
    
    Performs chi-squared test to determine if high pLLPS proteins preferentially
    interact with each other compared to random expectation.
    
    Parameters
    ----------
    matched_df : pd.DataFrame
        DataFrame with 'pllps_A' and 'pllps_B' columns
    threshold : float, optional
        Threshold for classifying high vs low pLLPS (default: 0.7)
    
    Returns
    -------
    Optional[Dict]
        Dictionary with enrichment statistics including:
        - total: Total number of interactions
        - high_high, high_low, low_low: Counts of each interaction type
        - enrichment: Enrichment factor (observed/expected)
        - chi2, p_value: Chi-squared test statistics
        - expected_hh, expected_hl, expected_ll: Expected percentages
    
    Examples
    --------
    >>> results = analyze_interaction_enrichment(matched_df, threshold=0.7)
    >>> if results and results['p_value'] < 0.05:
    ...     print(f"Significant enrichment: {results['enrichment']:.2f}x")
    """
    if len(matched_df) == 0:
        return None
    
    # Classify interactions
    matched_df = matched_df.copy()
    matched_df['class_A'] = matched_df['pllps_A'] >= threshold
    matched_df['class_B'] = matched_df['pllps_B'] >= threshold
    
    # Count interaction types
    high_high = len(matched_df[(matched_df['class_A']) & (matched_df['class_B'])])
    high_low = len(matched_df[(matched_df['class_A']) & (~matched_df['class_B'])]) + \
               len(matched_df[(~matched_df['class_A']) & (matched_df['class_B'])])
    low_low = len(matched_df[(~matched_df['class_A']) & (~matched_df['class_B'])])
    
    total = high_high + high_low + low_low
    
    if total == 0:
        return None
    
    # Calculate proportions for proteins in the interaction network
    all_proteins = pd.concat([
        matched_df[['preferredName_A', 'pllps_A']].rename(columns={'preferredName_A': 'protein', 'pllps_A': 'pllps'}),
        matched_df[['preferredName_B', 'pllps_B']].rename(columns={'preferredName_B': 'protein', 'pllps_B': 'pllps'})
    ]).drop_duplicates()
    
    n_high = len(all_proteins[all_proteins['pllps'] >= threshold])
    n_total = len(all_proteins)
    p_high = n_high / n_total if n_total > 0 else 0
    p_low = 1 - p_high
    
    # Calculate expected counts under random interaction model
    expected_hh = (p_high ** 2) * 100
    expected_hl = (2 * p_high * p_low) * 100
    expected_ll = (p_low ** 2) * 100
    
    # Chi-squared test
    observed = np.array([high_high, high_low, low_low])
    expected = np.array([expected_hh * total / 100, expected_hl * total / 100, expected_ll * total / 100])
    
    # Avoid division by zero
    if np.any(expected == 0):
        chi2, p_value = None, None
    else:
        chi2, p_value = chi2_contingency([observed, expected])[:2]
    
    # Calculate enrichment factor
    enrichment = (high_high / total * 100) / expected_hh if expected_hh > 0 else 0
    
    return {
        'total': total,
        'high_high': high_high,
        'high_low': high_low,
        'low_low': low_low,
        'enrichment': enrichment,
        'chi2': chi2,
        'p_value': p_value,
        'expected_hh': expected_hh,
        'expected_hl': expected_hl,
        'expected_ll': expected_ll,
        'p_high': p_high,
        'p_low': p_low,
        'n_high_proteins': n_high,
        'n_total_proteins': n_total
    }


def save_interactions_to_cache(
    interactions_df: pd.DataFrame,
    score_threshold: int = 700,
    output_dir: str = "data"
) -> str:
    """
    Save STRING interactions to cache file for offline use.
    
    Parameters
    ----------
    interactions_df : pd.DataFrame
        DataFrame with interaction data from STRING
    score_threshold : int, optional
        Score threshold used (for filename)
    output_dir : str, optional
        Directory to save cache file (default: 'data')
    
    Returns
    -------
    str
        Path to saved cache file
    
    Examples
    --------
    >>> cache_path = save_interactions_to_cache(interactions_df, score_threshold=700)
    >>> print(f"Cache saved to: {cache_path}")
    """
    import json
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    cache_file = output_path / f"string_cache_{score_threshold}.json"
    
    # Convert DataFrame to list of dicts for JSON
    interactions_list = interactions_df.to_dict('records')
    
    with open(cache_file, 'w') as f:
        json.dump(interactions_list, f, indent=2)
    
    return str(cache_file)


# For backward compatibility with existing code
def get_string_interactions(*args, **kwargs):
    """Alias for fetch_string_interactions (backward compatibility)"""
    return fetch_string_interactions(*args, **kwargs)
