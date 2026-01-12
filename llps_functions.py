"""
LLPS Analysis Master Functions Module

This module consolidates all functions for analyzing Liquid-Liquid Phase Separation (LLPS) 
protein data and their interaction networks. All analysis notebooks and scripts should 
import from this master file.

USAGE:
------
In Python scripts:
    from llps_functions import load_llps_data, fetch_string_interactions
    
In Jupyter notebooks:
    import llps_functions as lf
    df = lf.load_llps_data('path/to/data.xlsx')

Organized sections:
1. Data Loading and Classification
2. Location Parsing and Analysis
3. STRING Database Interactions
4. Network Analysis
5. Enrichment Analysis
6. Visualization Functions
7. Export and Caching Functions
8. Result Saving and Loading Utilities
9. Functional Classification Utilities

Author: LLPS Analysis Team
Date: 2025
"""

import pandas as pd
import numpy as np
import requests
import time
import matplotlib.pyplot as plt
import seaborn as sns
import re
import json
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Callable, Any
from scipy.stats import chi2_contingency
from scipy import stats
import warnings


# =============================================================================
# 1. DATA LOADING AND CLASSIFICATION
# =============================================================================

def load_llps_data(filepath: str = None) -> pd.DataFrame:
    """
    Load LLPS protein data from Excel file.
    
    Parameters
    ----------
    filepath : str, optional
        Path to the Excel file. If None, tries to find the default data file.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with protein data including p(LLPS) scores.
    """
    if filepath is None:
        # Try to find default files
        possible_paths = [
            Path(__file__).parent / "Human Phase separation data.xlsx",
            Path(__file__).parent / "data" / "sample_data.xlsx",
            Path("Human Phase separation data.xlsx"),
            Path("data/sample_data.xlsx")
        ]
        for path in possible_paths:
            if path.exists():
                filepath = path
                break
    
    if filepath is None:
        raise FileNotFoundError("Could not find LLPS data file. Please specify filepath.")
    
    # Validate file exists
    if not Path(filepath).exists():
        raise FileNotFoundError(f"Data file not found: {filepath}")
    
    df = pd.read_excel(filepath, engine='openpyxl')
    print(f"✅ Loaded {len(df)} proteins from {filepath}")
    
    if 'p(LLPS)' in df.columns:
        print(f"   p(LLPS) range: {df['p(LLPS)'].min():.3f} - {df['p(LLPS)'].max():.3f}")
    
    return df


def load_and_classify_data(filepath, high_threshold=0.7, low_threshold=0.4):
    """
    Load the dataset and classify proteins into High, Medium, and Low pLLPS classes.
    
    Parameters
    ----------
    filepath : str
        Path to Excel file with protein data
    high_threshold : float
        Threshold for High pLLPS classification (default: 0.7)
    low_threshold : float
        Threshold for Low pLLPS classification (default: 0.4)
    
    Returns
    -------
    pd.DataFrame
        DataFrame with added pLLPS_class column
    """
    df = pd.read_excel(filepath)
    print(f"Total proteins: {len(df)}")
    print(f"p(LLPS) range: {df['p(LLPS)'].min():.3f} - {df['p(LLPS)'].max():.3f}")
    
    # Classify proteins
    conditions = [
        (df['p(LLPS)'] >= high_threshold),
        (df['p(LLPS)'] >= low_threshold) & (df['p(LLPS)'] < high_threshold),
        (df['p(LLPS)'] < low_threshold)
    ]
    choices = ['High', 'Medium', 'Low']
    df['pLLPS_class'] = np.select(conditions, choices, default='Low')
    
    # Count classes
    counts = df['pLLPS_class'].value_counts()
    print("Protein Classification Counts:")
    print(counts)
    print(f"\nHigh: >={high_threshold}")
    print(f"Medium: {low_threshold}-{high_threshold}")
    print(f"Low: <{low_threshold}")
    
    return df


def get_high_pllps_proteins(
    df: pd.DataFrame,
    threshold: float = 0.7,
    method: str = 'absolute'
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Identify high pLLPS proteins based on threshold.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with protein data.
    threshold : float
        Threshold for high pLLPS. Default is 0.7.
    method : str
        'absolute' for fixed threshold, 'percentile' for percentile-based.
        If 'percentile', threshold should be 0-100.
    
    Returns
    -------
    Tuple[pd.DataFrame, List[str]]
        DataFrame with pLLPS_class column added, and list of high pLLPS Entry IDs.
    """
    df = df.copy()
    
    if 'p(LLPS)' not in df.columns:
        raise ValueError("DataFrame must contain 'p(LLPS)' column")
    
    # Validate required columns
    required_cols = ['Entry', 'p(LLPS)']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    if method == 'percentile':
        actual_threshold = df['p(LLPS)'].quantile(threshold / 100)
        print(f"Using {threshold}th percentile: pLLPS >= {actual_threshold:.3f}")
        threshold = actual_threshold
    
    df['pLLPS_class'] = np.where(df['p(LLPS)'] >= threshold, 'high', 'low')
    
    high_count = (df['pLLPS_class'] == 'high').sum()
    low_count = (df['pLLPS_class'] == 'low').sum()
    
    print(f"\n📊 pLLPS Classification (threshold = {threshold}):")
    print(f"   High pLLPS: {high_count} proteins ({100*high_count/len(df):.1f}%)")
    print(f"   Low pLLPS:  {low_count} proteins ({100*low_count/len(df):.1f}%)")
    
    # Get high pLLPS proteins sorted by pLLPS (highest first)
    high_pllps_df = df[df['pLLPS_class'] == 'high'].sort_values('p(LLPS)', ascending=False)
    high_pllps_ids = high_pllps_df['Entry'].tolist()
    
    return df, high_pllps_ids


# =============================================================================
# 2. LOCATION PARSING AND ANALYSIS
# =============================================================================

def parse_location(location_str):
    """
    Parse a subcellular location string from UniProt and return a list of location terms.
    
    Parameters
    ----------
    location_str : str
        Subcellular location string from UniProt
    
    Returns
    -------
    list
        List of parsed location terms
    """
    if pd.isna(location_str) or location_str == '':
        return []
    
    location_str = str(location_str)
    
    # Remove curly bracket content like {ECO:0000269|PubMed:12345}
    location_str = re.sub(r'\{[^}]*\}', '', location_str)
    
    # Remove all square bracket content, optionally followed by colon
    # This handles [Isoform 1]: and [Note: ...]
    location_str = re.sub(r'\[[^\]]*\]:?', '', location_str)
    
    # Remove Isoform ...: tags (without brackets)
    location_str = re.sub(r'Isoform\s+[^:]+:\s*', '', location_str, flags=re.IGNORECASE)
    
    # Remove normal bracket content like (By similarity) or (Potential)
    location_str = re.sub(r'\([^)]*\)', '', location_str)
    
    # Remove "SUBCELLULAR LOCATION:" prefix if present
    location_str = re.sub(r'^SUBCELLULAR LOCATION:\s*', '', location_str, flags=re.IGNORECASE)
    
    # Remove everything after 'Note=' (case-insensitive)
    location_str = re.sub(r'\s*Note=.*', '', location_str, flags=re.IGNORECASE)
    
    # Split by common separators and clean up
    # UniProt uses semicolons, commas, and periods as separators
    parts = re.split(r'[;,.]', location_str)
    
    locations = []
    for part in parts:
        # Clean up whitespace and filter empty/very short strings
        cleaned = part.strip()
        if len(cleaned) < 2:  # Skip empty or single-char remnants
            continue
            
        # Normalize case: Ensure first letter is uppercase
        # This fixes issues with "cytosol" vs "Cytosol"
        cleaned = cleaned[0].upper() + cleaned[1:]
        
        if cleaned not in locations:
            locations.append(cleaned)
    
    return locations


def add_location_columns(df):
    """
    Add parsed location columns to the dataframe.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with 'Subcellular location [CC]' column
    
    Returns
    -------
    pd.DataFrame
        DataFrame with added 'Location Categories' column
    """
    if 'Subcellular location [CC]' not in df.columns:
        return df
    
    df = df.copy()
    # Add column with list of all parsed locations
    df['Location Categories'] = df['Subcellular location [CC]'].apply(parse_location)
    return df


def analyze_interactions_by_location(matched_df, full_dataset_df, locations, high_threshold=0.7, low_threshold=0.4):
    """
    Analyze interaction preferences for specific subcellular locations.
    
    Parameters
    ----------
    matched_df : pd.DataFrame
        DataFrame with matched interactions
    full_dataset_df : pd.DataFrame
        Full dataset DataFrame with all proteins
    locations : list
        List of location names to analyze
    high_threshold : float
        Threshold for High pLLPS classification
    low_threshold : float
        Threshold for Low pLLPS classification
    
    Returns
    -------
    dict
        Dictionary with location-specific results
    """
    results = {}
    
    # Ensure location columns exist
    if 'Location Categories' not in full_dataset_df.columns:
        full_dataset_df = add_location_columns(full_dataset_df)
        
    # We also need location info for the matched_df proteins
    # We can merge it from full_dataset_df
    if 'Location Categories' not in matched_df.columns:
        # Create a map from Entry -> Location Categories
        loc_map = dict(zip(full_dataset_df['Entry'], full_dataset_df['Location Categories']))
        
        # Helper to get locations for a Uniprot ID
        def get_locs(uid):
            return loc_map.get(uid, [])
            
        matched_df['locs_a'] = matched_df['uniprot_a'].apply(get_locs)
        matched_df['locs_b'] = matched_df['uniprot_b'].apply(get_locs)
    
    for loc in locations:
        print(f"\nAnalyzing location: {loc}")
        
        # Filter full dataset to this location (for background probs)
        loc_full_df = full_dataset_df[full_dataset_df['Location Categories'].apply(lambda x: loc in x)]
        
        if len(loc_full_df) < 10:
            print(f"  Skipping {loc}: Too few proteins ({len(loc_full_df)})")
            continue
            
        # Filter interactions where BOTH proteins are in this location
        # This represents the "Nucleus Interaction Network".
        
        loc_matched_df = matched_df[
            matched_df['locs_a'].apply(lambda x: loc in x) & 
            matched_df['locs_b'].apply(lambda x: loc in x)
        ].copy()
        
        print(f"  Proteins in location: {len(loc_full_df)}")
        print(f"  Interactions within location: {len(loc_matched_df)}")
        
        # Count entries per class in this location
        class_counts = loc_full_df['pLLPS_class'].value_counts()
        print(f"  Class breakdown: High={class_counts.get('High', 0)}, Medium={class_counts.get('Medium', 0)}, Low={class_counts.get('Low', 0)}")
        
        if len(loc_matched_df) < 10:
            print("  Not enough interactions for analysis.")
            continue
            
        # Run matrix analysis
        loc_results = analyze_interaction_matrix(loc_matched_df, loc_full_df, high_threshold, low_threshold)
        if loc_results:
            results[loc] = loc_results
            
    return results


# =============================================================================
# 3. STRING DATABASE INTERACTIONS
# =============================================================================

def get_string_mapping(protein_ids, species=9606, batch_size=2000):
    """
    Map Uniprot IDs to STRING preferred names (Gene Names) and STRING IDs.
    
    Parameters
    ----------
    protein_ids : list
        List of UniProt protein IDs
    species : int
        NCBI taxonomy ID (9606 = human)
    batch_size : int
        Number of proteins per API request
    
    Returns
    -------
    tuple
        Two dictionaries:
        1. string_to_uniprot: Maps STRING ID/Name -> Uniprot ID
        2. uniprot_to_string: Maps Uniprot ID -> STRING ID
    """
    url = "https://string-db.org/api/json/get_string_ids"
    string_to_uniprot = {}
    uniprot_to_string = {}
    
    total_batches = (len(protein_ids) + batch_size - 1) // batch_size
    print(f"Mapping {len(protein_ids)} proteins to STRING names in {total_batches} batches...")
    
    for i in range(0, len(protein_ids), batch_size):
        batch = protein_ids[i:i+batch_size]
        batch_num = i // batch_size + 1
        
        params = {
            "identifiers": "\r".join(batch),
            "species": species,
            "echo_query": 1
        }
        
        try:
            response = requests.post(url, data=params, timeout=60)
            if response.status_code == 200:
                data = response.json()
                for item in data:
                    # Map preferredName -> queryItem (Uniprot ID)
                    # Also map stringId -> queryItem
                    u_id = item['queryItem']
                    pref_name = item['preferredName']
                    string_id = item['stringId']
                    
                    string_to_uniprot[pref_name] = u_id
                    string_to_uniprot[string_id] = u_id
                    
                    # Map Uniprot -> STRING ID (preferred for querying)
                    uniprot_to_string[u_id] = string_id
                    
                print(f"  Batch {batch_num}/{total_batches}: Mapped {len(data)} items")
            else:
                print(f"  Batch {batch_num}/{total_batches}: Error {response.status_code}")
        except Exception as e:
            print(f"  Batch {batch_num}/{total_batches}: {e}")
            
        time.sleep(1)
        
    return string_to_uniprot, uniprot_to_string


def fetch_string_interactions(
    protein_ids: List[str],
    species: int = 9606,
    score_threshold: int = 700,
    batch_size: int = 100,
    progress_callback: Optional[Callable[[str], None]] = None,
    use_cache: bool = True,
    network_type: str = "physical"
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
    network_type : str, optional
        "physical" for experimentally determined, "functional" for all (default: "physical")
    
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
    
    if progress_callback:
        progress_callback(f"Fetching batch 1/{total_batches}...")
    
    # Fetch interactions from STRING API
    for i in range(0, len(protein_ids), batch_size):
        batch = protein_ids[i:i+batch_size]
        batch_num = i // batch_size + 1
        
        if progress_callback and batch_num > 1:
            progress_callback(f"Fetching batch {batch_num}/{total_batches}...")
        
        params = {
            "identifiers": "\r".join(batch),
            "species": species,
            "required_score": score_threshold,
            "network_type": network_type,
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


def fetch_interaction_partners(protein_ids, uniprot_map=None, species=9606, score_threshold=700, batch_size=100):
    """
    Fetch interaction partners for a list of proteins from STRING.
    
    Parameters
    ----------
    protein_ids : list
        List of UniProt IDs
    uniprot_map : dict, optional
        Dictionary mapping Uniprot IDs to STRING IDs (optional but recommended)
    species : int
        NCBI taxonomy ID (9606 = human)
    score_threshold : int
        Minimum confidence score (0-1000)
    batch_size : int
        Proteins per API request
    
    Returns
    -------
    pd.DataFrame
        DataFrame with interactions
    """
    # Use interaction_partners API to get ALL partners, not just interactions within the set
    string_api_url = "https://string-db.org/api/json/interaction_partners"
    all_interactions = []
    
    # Convert Uniprot IDs to STRING IDs if map provided
    # This improves matching success rate significantly
    query_ids = []
    if uniprot_map:
        for pid in protein_ids:
            if pid in uniprot_map:
                query_ids.append(uniprot_map[pid])
            else:
                # If no map, try using the ID directly (might fail if STRING doesn't recognize it)
                query_ids.append(pid)
    else:
        query_ids = protein_ids
        
    total_batches = (len(query_ids) + batch_size - 1) // batch_size
    print(f"Fetching interactions for {len(query_ids)} proteins in {total_batches} batches...")
    
    for i in range(0, len(query_ids), batch_size):
        batch = query_ids[i:i+batch_size]
        batch_num = i // batch_size + 1
        
        params = {
            "identifiers": "\r".join(batch),
            "species": species,
            "required_score": score_threshold,
            "caller_identity": "pllps_analysis",
            "limit": 50  # Limit partners per protein to avoid explosion
        }
        
        try:
            response = requests.post(string_api_url, data=params, timeout=60)
            if response.status_code == 200:
                interactions = response.json()
                all_interactions.extend(interactions)
                print(f"  Batch {batch_num}/{total_batches}: {len(interactions)} interactions")
            else:
                print(f"  Batch {batch_num}/{total_batches}: Error {response.status_code}")
        except Exception as e:
            print(f"  Batch {batch_num}/{total_batches}: {e}")
        
        time.sleep(1)  # Rate limiting
    
    if all_interactions:
        return pd.DataFrame(all_interactions)
    return pd.DataFrame()


def load_string_network_file(filepath: str) -> pd.DataFrame:
    """
    Load STRING network data from downloaded TSV file.
    
    Parameters
    ----------
    filepath : str
        Path to TSV file downloaded from STRING.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with interaction data.
    """
    df = pd.read_csv(filepath, sep='\t')
    print(f"✅ Loaded {len(df)} interactions from {filepath}")
    return df


# =============================================================================
# 4. NETWORK ANALYSIS
# =============================================================================

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


def match_interactors_to_pllps(interactions_df, pllps_df, string_map=None):
    """
    Match interaction partners to the pLLPS dataset.
    
    Returns DataFrame with both proteins' pLLPS values.
    
    Parameters
    ----------
    interactions_df : pd.DataFrame
        DataFrame with interaction data from STRING
    pllps_df : pd.DataFrame
        DataFrame with pLLPS scores
    string_map : dict, optional
        Dictionary mapping STRING IDs to UniProt IDs
    
    Returns
    -------
    pd.DataFrame
        DataFrame with matched pLLPS scores
    """
    if len(interactions_df) == 0:
        return pd.DataFrame()
    
    # Create lookup dictionary: Entry -> pLLPS
    pllps_lookup = dict(zip(pllps_df['Entry'], pllps_df['p(LLPS)']))
    
    # Also create lookup by Entry name (STRING often uses gene names)
    entry_name_lookup = dict(zip(pllps_df['Entry name'], pllps_df['p(LLPS)']))
    
    results = []
    
    for _, row in interactions_df.iterrows():
        protein_a_name = row.get('preferredName_A', '')
        protein_a_id = row.get('stringId_A', '')
        
        protein_b_name = row.get('preferredName_B', '')
        protein_b_id = row.get('stringId_B', '')
        
        score = row.get('score', 0)
        
        # Resolve Protein A to Uniprot ID
        uniprot_a = None
        if string_map:
            uniprot_a = string_map.get(protein_a_name) or string_map.get(protein_a_id)
        
        # Fallback if no map or not found
        if not uniprot_a:
            if protein_a_name in pllps_lookup:
                uniprot_a = protein_a_name
            elif protein_a_name in entry_name_lookup:
                uniprot_a = protein_a_name
        
        # Resolve Protein B to Uniprot ID
        uniprot_b = None
        if string_map:
            uniprot_b = string_map.get(protein_b_name) or string_map.get(protein_b_id)
            
        if not uniprot_b:
            if protein_b_name in pllps_lookup:
                uniprot_b = protein_b_name

        # Get pLLPS scores
        pllps_a = pllps_lookup.get(uniprot_a)
        if pllps_a is None and not string_map:
            pllps_a = pllps_lookup.get(protein_a_name) or entry_name_lookup.get(protein_a_name)

        pllps_b = pllps_lookup.get(uniprot_b)
        if pllps_b is None and not string_map:
            pllps_b = pllps_lookup.get(protein_b_name) or entry_name_lookup.get(protein_b_name)
        
        results.append({
            'protein_a': protein_a_name,
            'protein_b': protein_b_name,
            'uniprot_a': uniprot_a,
            'uniprot_b': uniprot_b,
            'score': score,
            'pllps_a': pllps_a,
            'pllps_b': pllps_b
        })
    
    return pd.DataFrame(results)


def analyze_network(
    interactions_df: pd.DataFrame,
    pllps_df: pd.DataFrame,
    high_threshold: float = 0.7
) -> Tuple[Dict, 'nx.Graph']:
    """
    Analyze the protein interaction network with pLLPS annotations.
    
    Parameters
    ----------
    interactions_df : pd.DataFrame
        DataFrame with interaction data (from STRING).
    pllps_df : pd.DataFrame
        DataFrame with protein pLLPS scores.
    high_threshold : float
        Threshold for high pLLPS classification.
    
    Returns
    -------
    Tuple[Dict, nx.Graph]
        Tuple containing:
        - Dict: Dictionary with network analysis results
        - nx.Graph: NetworkX graph object with pLLPS values as node attributes.
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx library required. Install with: pip install networkx")
    
    # Build network
    G = nx.Graph()
    
    # Determine column names for protein IDs
    if 'preferredName_A' in interactions_df.columns:
        col_a, col_b = 'preferredName_A', 'preferredName_B'
    elif 'stringId_A' in interactions_df.columns:
        col_a, col_b = 'stringId_A', 'stringId_B'
    else:
        # Try to find appropriate columns
        str_cols = interactions_df.select_dtypes(include=['object']).columns
        col_a, col_b = str_cols[0], str_cols[1]
        print(f"Using columns: {col_a}, {col_b}")
    
    # Determine score column
    if 'score' in interactions_df.columns:
        score_col = 'score'
    elif 'combined_score' in interactions_df.columns:
        score_col = 'combined_score'
    else:
        score_col = None
        print("⚠️  No score column found in interactions data. Edges will have no weight.")
    
    # Add edges
    for _, row in interactions_df.iterrows():
        edge_score = row[score_col] if score_col and pd.notna(row.get(score_col)) else None
        G.add_edge(
            row[col_a],
            row[col_b],
            score=edge_score
        )
    
    # Create pLLPS lookup
    pllps_dict = dict(zip(pllps_df['Entry'], pllps_df['p(LLPS)']))
    
    # Also try Entry name for matching
    if 'Entry name' in pllps_df.columns:
        entry_name_dict = dict(zip(pllps_df['Entry name'], pllps_df['p(LLPS)']))
        pllps_dict.update(entry_name_dict)
    
    # Add pLLPS as node attribute
    for node in G.nodes():
        # Try exact match first, then try to extract gene name from STRING ID
        pllps_value = pllps_dict.get(node)
        if pllps_value is None and '.' in str(node):
            # STRING IDs often have format SPECIES.PROTEIN
            potential_id = str(node).split('.')[-1]
            pllps_value = pllps_dict.get(potential_id)
        G.nodes[node]['pLLPS'] = pllps_value
    
    # Basic network metrics
    results = {
        'total_nodes': G.number_of_nodes(),
        'total_edges': G.number_of_edges(),
        'density': nx.density(G),
        'avg_clustering': nx.average_clustering(G) if G.number_of_nodes() > 0 else 0,
    }
    
    # Classify nodes
    high_pllps_nodes = []
    low_pllps_nodes = []
    unknown_nodes = []
    
    for node in G.nodes():
        pllps = G.nodes[node].get('pLLPS')
        if pllps is not None:
            if pllps >= high_threshold:
                high_pllps_nodes.append(node)
            else:
                low_pllps_nodes.append(node)
        else:
            unknown_nodes.append(node)
    
    results['high_pllps_nodes'] = len(high_pllps_nodes)
    results['low_pllps_nodes'] = len(low_pllps_nodes)
    results['unknown_pllps_nodes'] = len(unknown_nodes)
    
    # Analyze high pLLPS subnetwork
    if len(high_pllps_nodes) > 1:
        high_subgraph = G.subgraph(high_pllps_nodes)
        results['high_pllps_edges'] = high_subgraph.number_of_edges()
        results['high_pllps_density'] = nx.density(high_subgraph)
        results['high_pllps_avg_clustering'] = nx.average_clustering(high_subgraph)
    else:
        results['high_pllps_edges'] = 0
        results['high_pllps_density'] = 0
        results['high_pllps_avg_clustering'] = 0
    
    # Interaction type analysis
    high_set = set(high_pllps_nodes)
    low_set = set(low_pllps_nodes)
    
    high_high = 0
    high_low = 0
    low_low = 0
    
    for edge in G.edges():
        a, b = edge
        a_high = a in high_set
        b_high = b in high_set
        
        if a_high and b_high:
            high_high += 1
        elif a_high or b_high:
            high_low += 1
        elif a in low_set and b in low_set:
            low_low += 1
    
    results['high_high_interactions'] = high_high
    results['high_low_interactions'] = high_low
    results['low_low_interactions'] = low_low
    
    # Calculate expected vs observed
    n_classified = len(high_pllps_nodes) + len(low_pllps_nodes)
    if n_classified > 0:
        p_high = len(high_pllps_nodes) / n_classified
        expected_high_high_ratio = p_high * p_high
        
        total_classified_edges = high_high + high_low + low_low
        if total_classified_edges > 0:
            observed_high_high_ratio = high_high / total_classified_edges
            results['enrichment_ratio'] = (
                observed_high_high_ratio / expected_high_high_ratio 
                if expected_high_high_ratio > 0 else 0
            )
        else:
            results['enrichment_ratio'] = 0
    else:
        results['enrichment_ratio'] = 0
    
    # Degree analysis
    if G.number_of_nodes() > 0:
        degrees = dict(G.degree())
        
        high_degrees = [degrees[n] for n in high_pllps_nodes if n in degrees]
        low_degrees = [degrees[n] for n in low_pllps_nodes if n in degrees]
        
        results['avg_degree_high_pllps'] = np.mean(high_degrees) if high_degrees else 0
        results['avg_degree_low_pllps'] = np.mean(low_degrees) if low_degrees else 0
    else:
        results['avg_degree_high_pllps'] = 0
        results['avg_degree_low_pllps'] = 0
    
    return results, G


# =============================================================================
# 5. ENRICHMENT ANALYSIS
# =============================================================================

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
        DataFrame with 'pllps_a' and 'pllps_b' columns
    threshold : float, optional
        Threshold for classifying high vs low pLLPS (default: 0.7)
    
    Returns
    -------
    Optional[Dict]
        Dictionary with enrichment statistics
    
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
    
    # Handle both column name conventions
    if 'pllps_A' in matched_df.columns:
        matched_df['class_A'] = matched_df['pllps_A'] >= threshold
        matched_df['class_B'] = matched_df['pllps_B'] >= threshold
        pllps_a_col, pllps_b_col = 'pllps_A', 'pllps_B'
        name_a_col, name_b_col = 'preferredName_A', 'preferredName_B'
    else:
        matched_df['class_A'] = matched_df['pllps_a'] >= threshold
        matched_df['class_B'] = matched_df['pllps_b'] >= threshold
        pllps_a_col, pllps_b_col = 'pllps_a', 'pllps_b'
        name_a_col = 'protein_a' if 'protein_a' in matched_df.columns else 'preferredName_A'
        name_b_col = 'protein_b' if 'protein_b' in matched_df.columns else 'preferredName_B'
    
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
        matched_df[[name_a_col, pllps_a_col]].rename(columns={name_a_col: 'protein', pllps_a_col: 'pllps'}),
        matched_df[[name_b_col, pllps_b_col]].rename(columns={name_b_col: 'protein', pllps_b_col: 'pllps'})
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


def analyze_interaction_matrix(matched_df, full_dataset_df, high_threshold=0.7, low_threshold=0.4):
    """
    Analyze interaction preferences between High, Medium, and Low pLLPS classes.
    Creates a 3x3 Enrichment Matrix.
    
    Parameters
    ----------
    matched_df : pd.DataFrame
        DataFrame with matched interactions
    full_dataset_df : pd.DataFrame
        Full dataset with all proteins
    high_threshold : float
        Threshold for High pLLPS classification
    low_threshold : float
        Threshold for Low pLLPS classification
    
    Returns
    -------
    dict or None
        Dictionary with enrichment analysis results
    """
    # Filter to interactions where both proteins have pLLPS scores
    complete_df = matched_df.dropna(subset=['pllps_a', 'pllps_b']).copy()
    
    if len(complete_df) == 0:
        print("No interactions with complete pLLPS data found.")
        return None
    
    # Define classes for the interaction pairs
    def get_class(score):
        if score >= high_threshold: return 'High'
        elif score >= low_threshold: return 'Medium'
        else: return 'Low'
        
    complete_df['class_a'] = complete_df['pllps_a'].apply(get_class)
    complete_df['class_b'] = complete_df['pllps_b'].apply(get_class)
    
    # Calculate Genomic Background Probabilities
    total_proteins = len(full_dataset_df)
    class_counts = full_dataset_df['pLLPS_class'].value_counts()
    p_genome = {
        'High': class_counts.get('High', 0) / total_proteins,
        'Medium': class_counts.get('Medium', 0) / total_proteins,
        'Low': class_counts.get('Low', 0) / total_proteins
    }
    
    print(f"\n{'='*50}")
    print("GENOMIC BACKGROUND PROBABILITIES")
    print(f"{'='*50}")
    for cls, p in p_genome.items():
        print(f"P({cls}): {p:.3f}")
        
    # Initialize Matrices
    classes = ['High', 'Medium', 'Low']
    observed_counts = pd.DataFrame(0, index=classes, columns=classes)
    expected_counts = pd.DataFrame(0.0, index=classes, columns=classes)
    
    # Count Observed Interactions
    for _, row in complete_df.iterrows():
        c1, c2 = row['class_a'], row['class_b']
        observed_counts.loc[c1, c2] += 1
        observed_counts.loc[c2, c1] += 1
        
    # Calculate Expected Counts
    total_observations = observed_counts.sum().sum()
    
    for c1 in classes:
        for c2 in classes:
            expected_counts.loc[c1, c2] = total_observations * p_genome[c1] * p_genome[c2]
            
    # Calculate Enrichment (Observed / Expected)
    enrichment_matrix = observed_counts / expected_counts
    
    print(f"\n{'='*50}")
    print("INTERACTION ENRICHMENT MATRIX (Observed / Expected)")
    print(f"{'='*50}")
    print(enrichment_matrix.round(2))
    
    return {
        'observed': observed_counts,
        'expected': expected_counts,
        'enrichment': enrichment_matrix,
        'p_genome': p_genome,
        'total_interactions': len(complete_df)
    }


# =============================================================================
# 6. VISUALIZATION FUNCTIONS
# =============================================================================

def plot_interaction_heatmap(results, output_file='pllps_interaction_matrix.png'):
    """
    Visualize the interaction enrichment matrix as a heatmap.
    
    Parameters
    ----------
    results : dict
        Results dictionary from analyze_interaction_matrix
    output_file : str
        Path to save the output figure
    """
    if results is not None:
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot Heatmap of Enrichment
        sns.heatmap(results['enrichment'], annot=True, fmt=".2f", cmap="RdBu_r", center=1,
                    linewidths=.5, ax=ax, vmin=0, vmax=3)
        
        ax.set_title('Interaction Enrichment Matrix\n(>1 = Enriched, <1 = Depleted)', fontsize=14)
        ax.set_xlabel('Partner Protein Class', fontsize=12)
        ax.set_ylabel('Query Protein Class', fontsize=12)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=150)
        plt.show()
        
        print(f"\nFigure saved to: {output_file}")
        
        # Print interpretation
        print("\nInterpretation:")
        print("- Values > 1.0 (Red): These classes interact MORE than expected by chance.")
        print("- Values < 1.0 (Blue): These classes interact LESS than expected by chance.")
        print("- Diagonal (High-High, Med-Med, Low-Low): Shows homotypic preference.")


def plot_location_heatmaps(results):
    """
    Plot a grid of heatmaps for location-specific results.
    
    Parameters
    ----------
    results : dict
        Dictionary with location-specific enrichment results
    """
    if not results:
        print("No results to plot.")
        return
        
    n_locs = len(results)
    cols = 3
    rows = (n_locs + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows))
    if n_locs > 1:
        axes = axes.flatten()
    else:
        axes = [axes]
    
    for i, (loc, res) in enumerate(results.items()):
        ax = axes[i]
        sns.heatmap(res['enrichment'], annot=True, fmt=".2f", cmap="RdBu_r", center=1,
                    linewidths=.5, ax=ax, vmin=0, vmax=3, cbar=False)
        ax.set_title(f"{loc}\n(n={res['total_interactions']})")
        ax.set_xlabel('')
        ax.set_ylabel('')
        
    # Hide empty subplots
    for j in range(i+1, len(axes)):
        axes[j].axis('off')
        
    plt.tight_layout()
    plt.savefig('pllps_location_analysis.png', dpi=150)
    plt.show()
    print("Figure saved to: pllps_location_analysis.png")


def print_analysis_report(results: Dict):
    """
    Print formatted analysis report.
    
    Parameters
    ----------
    results : dict
        Dictionary with network analysis results
    """
    print("\n" + "="*60)
    print("📊 NETWORK ANALYSIS REPORT")
    print("="*60)
    
    print("\n🔗 Network Overview:")
    print(f"   Total nodes: {results['total_nodes']}")
    print(f"   Total edges: {results['total_edges']}")
    print(f"   Network density: {results['density']:.4f}")
    print(f"   Average clustering coefficient: {results['avg_clustering']:.4f}")
    
    print("\n🧬 pLLPS Classification:")
    print(f"   High pLLPS nodes: {results['high_pllps_nodes']}")
    print(f"   Low pLLPS nodes: {results['low_pllps_nodes']}")
    print(f"   Unknown pLLPS nodes: {results['unknown_pllps_nodes']}")
    
    print("\n🔬 High pLLPS Subnetwork:")
    print(f"   Edges among high pLLPS proteins: {results['high_pllps_edges']}")
    print(f"   Subnetwork density: {results['high_pllps_density']:.4f}")
    print(f"   Subnetwork avg clustering: {results['high_pllps_avg_clustering']:.4f}")
    
    print("\n⚡ Interaction Type Analysis:")
    total = results['high_high_interactions'] + results['high_low_interactions'] + results['low_low_interactions']
    if total > 0:
        print(f"   High-High: {results['high_high_interactions']} ({100*results['high_high_interactions']/total:.1f}%)")
        print(f"   High-Low:  {results['high_low_interactions']} ({100*results['high_low_interactions']/total:.1f}%)")
        print(f"   Low-Low:   {results['low_low_interactions']} ({100*results['low_low_interactions']/total:.1f}%)")
        print(f"\n   Enrichment of High-High interactions: {results['enrichment_ratio']:.2f}x")
    
    print("\n📈 Degree Analysis:")
    print(f"   Avg degree (high pLLPS): {results['avg_degree_high_pllps']:.2f}")
    print(f"   Avg degree (low pLLPS): {results['avg_degree_low_pllps']:.2f}")
    
    if results['avg_degree_low_pllps'] > 0:
        hub_ratio = results['avg_degree_high_pllps'] / results['avg_degree_low_pllps']
        print(f"   Ratio (high/low): {hub_ratio:.2f}x")
        if hub_ratio > 1.5:
            print("   ⚠️  High pLLPS proteins appear to be network hubs!")
    
    print("\n" + "="*60)


# =============================================================================
# 7. EXPORT AND CACHING FUNCTIONS
# =============================================================================

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


def export_protein_list(
    protein_ids: List[str],
    output_file: str = "high_pllps_proteins.txt"
) -> str:
    """
    Export protein IDs to a file for use with STRING web interface.
    
    Parameters
    ----------
    protein_ids : List[str]
        List of UniProt protein IDs.
    output_file : str
        Output filename.
    
    Returns
    -------
    str
        Path to the output file.
    """
    with open(output_file, 'w') as f:
        f.write('\n'.join(protein_ids))
    
    print(f"\n📄 Exported {len(protein_ids)} protein IDs to {output_file}")
    print(f"   Use this file with STRING web interface:")
    print(f"   1. Go to https://string-db.org/cgi/input")
    print(f"   2. Select 'Multiple proteins' tab")
    print(f"   3. Paste or upload the protein list")
    print(f"   4. Select 'Homo sapiens' as organism")
    print(f"   5. Download network as TSV for analysis")
    
    return output_file


# =============================================================================
# BACKWARD COMPATIBILITY ALIASES
# =============================================================================

# For backward compatibility with existing code
def get_string_interactions(*args, **kwargs):
    """Alias for fetch_string_interactions (backward compatibility)"""
    return fetch_string_interactions(*args, **kwargs)


# =============================================================================
# 8. RESULT SAVING AND LOADING UTILITIES
# =============================================================================

def save_analysis_result(data: Any, filename: str, results_dir: str = "results", 
                         format: str = "csv") -> Path:
    """
    Save analysis results to file.
    
    Parameters
    ----------
    data : Any
        Data to save (DataFrame, dict, list, etc.)
    filename : str
        Base filename (without extension)
    results_dir : str
        Directory to save results in (default: "results")
    format : str
        Format to save in: "csv", "json", "pickle" (default: "csv")
    
    Returns
    -------
    Path
        Path to saved file
    """
    results_path = Path(results_dir)
    results_path.mkdir(exist_ok=True)
    
    if format == "csv":
        filepath = results_path / f"{filename}.csv"
        if isinstance(data, pd.DataFrame):
            data.to_csv(filepath, index=False)
        else:
            raise ValueError("CSV format requires DataFrame input")
    elif format == "json":
        filepath = results_path / f"{filename}.json"
        with open(filepath, 'w') as f:
            if isinstance(data, (dict, list)):
                json.dump(data, f, indent=2, default=str)
            elif isinstance(data, pd.DataFrame):
                json.dump(data.to_dict(orient='records'), f, indent=2, default=str)
            else:
                raise ValueError("JSON format requires dict, list, or DataFrame")
    elif format == "pickle":
        filepath = results_path / f"{filename}.pkl"
        if isinstance(data, pd.DataFrame):
            data.to_pickle(filepath)
        else:
            import pickle
            with open(filepath, 'wb') as f:
                pickle.dump(data, f)
    else:
        raise ValueError(f"Unknown format: {format}")
    
    print(f"✅ Saved {format.upper()} to: {filepath}")
    return filepath


def load_analysis_result(filename: str, results_dir: str = "results", 
                         format: str = "csv") -> Any:
    """
    Load analysis results from file.
    
    Parameters
    ----------
    filename : str
        Base filename (with or without extension)
    results_dir : str
        Directory to load results from (default: "results")
    format : str
        Format to load: "csv", "json", "pickle" (default: "csv")
    
    Returns
    -------
    Any
        Loaded data (DataFrame, dict, list, etc.)
    """
    results_path = Path(results_dir)
    
    # Add extension if not present
    if not any(filename.endswith(ext) for ext in ['.csv', '.json', '.pkl']):
        if format == "csv":
            filename = f"{filename}.csv"
        elif format == "json":
            filename = f"{filename}.json"
        elif format == "pickle":
            filename = f"{filename}.pkl"
    
    filepath = results_path / filename
    
    if not filepath.exists():
        raise FileNotFoundError(f"Result file not found: {filepath}")
    
    if format == "csv" or filename.endswith('.csv'):
        data = pd.read_csv(filepath)
        print(f"✅ Loaded CSV from: {filepath} ({len(data)} rows)")
    elif format == "json" or filename.endswith('.json'):
        with open(filepath, 'r') as f:
            data = json.load(f)
        print(f"✅ Loaded JSON from: {filepath}")
    elif format == "pickle" or filename.endswith('.pkl'):
        import pickle
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
        print(f"✅ Loaded pickle from: {filepath}")
    else:
        raise ValueError(f"Cannot determine format for: {filename}")
    
    return data


def list_saved_results(results_dir: str = "results") -> List[str]:
    """
    List all saved result files.
    
    Parameters
    ----------
    results_dir : str
        Directory containing results (default: "results")
    
    Returns
    -------
    List[str]
        List of result filenames
    """
    results_path = Path(results_dir)
    if not results_path.exists():
        print(f"⚠️  Results directory not found: {results_path}")
        return []
    
    files = []
    for ext in ['*.csv', '*.json', '*.pkl']:
        files.extend([f.name for f in results_path.glob(ext)])
    
    if files:
        print(f"📁 Found {len(files)} result files in {results_path}:")
        for f in sorted(files):
            file_path = results_path / f
            size = file_path.stat().st_size
            size_str = f"{size/1024:.1f} KB" if size < 1024*1024 else f"{size/(1024*1024):.1f} MB"
            print(f"   - {f} ({size_str})")
    else:
        print(f"📁 No result files found in {results_path}")
    
    return sorted(files)


# =============================================================================
# 9. FUNCTIONAL CLASSIFICATION UTILITIES
# =============================================================================

def is_membrane_protein(function_str: str, protein_name_str: str = None, 
                       location_str: str = None) -> bool:
    """
    Determine if a protein is a membrane protein based on function annotation,
    protein name, and subcellular location.
    
    Parameters
    ----------
    function_str : str
        Function [CC] annotation
    protein_name_str : str, optional
        Protein names
    location_str : str, optional
        Subcellular location [CC] annotation
    
    Returns
    -------
    bool
        True if protein is classified as membrane protein
    """
    if not isinstance(function_str, str):
        function_str = ""
    if not isinstance(protein_name_str, str):
        protein_name_str = ""
    if not isinstance(location_str, str):
        location_str = ""
    
    # Combine all text for searching
    combined_text = f"{function_str} {protein_name_str} {location_str}".lower()
    
    # Membrane protein patterns
    membrane_patterns = [
        r'\bmulti-pass membrane protein\b',
        r'\bsingle-pass membrane protein\b',
        r'\bmembrane protein\b',
        r'\btransmembrane\b',
        r'\bintegral membrane\b',
        r'\bintrinsic membrane\b',
        r'\bmembrane-bound\b',
    ]
    
    for pattern in membrane_patterns:
        if re.search(pattern, combined_text):
            return True
    
    return False


def classify_protein_function(function_str: str, protein_name_str: str = None) -> List[str]:
    """
    Classify protein into functional categories based on annotations.
    
    Parameters
    ----------
    function_str : str
        Function [CC] annotation
    protein_name_str : str, optional
        Protein names
    
    Returns
    -------
    List[str]
        List of functional categories the protein belongs to
    """
    if not isinstance(function_str, str):
        function_str = ""
    if not isinstance(protein_name_str, str):
        protein_name_str = ""
    
    combined_text = f"{function_str} {protein_name_str}".lower()
    
    # Define functional categories
    functional_groups = {
        'Ion Channel': [
            r'\bion channel\b',
            r'\bsodium channel\b',
            r'\bpotassium channel\b',
            r'\bcalcium channel\b',
            r'\bchloride channel\b',
            r'\bchannel protein\b',
            r'\bvoltage-gated\b',
            r'\bligand-gated\b',
        ],
        'GPCR': [
            r'\bG protein-coupled receptor\b',
            r'\bGPCR\b',
            r'\breceptor protein-coupled\b',
            r'\b7-transmembrane\b',
        ],
        'Receptor Tyrosine Kinase': [
            r'\breceptor tyrosine kinase\b',
            r'\bRTK\b',
            r'\btyrosine-protein kinase receptor\b',
        ],
        'Transporter': [
            r'\btransporter\b',
            r'\bABC transporter\b',
            r'\bsolute carrier\b',
            r'\bSLC\b',
        ],
        'Enzyme': [
            r'\bkinase\b',
            r'\bphosphatase\b',
            r'\bprotease\b',
            r'\bligase\b',
            r'\btransferase\b',
        ],
        'Receptor': [
            r'\breceptor\b',
        ],
    }
    
    categories = []
    for category, patterns in functional_groups.items():
        for pattern in patterns:
            if re.search(pattern, combined_text):
                categories.append(category)
                break  # Only add category once
    
    return categories


def filter_membrane_proteins(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter DataFrame to only membrane proteins.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with protein data including 'Function [CC]', 'Protein names',
        and 'Subcellular location [CC]' columns
    
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with only membrane proteins
    """
    def check_membrane(row):
        return is_membrane_protein(
            row.get('Function [CC]', ''),
            row.get('Protein names', ''),
            row.get('Subcellular location [CC]', '')
        )
    
    membrane_mask = df.apply(check_membrane, axis=1)
    membrane_df = df[membrane_mask].copy()
    
    print(f"🔍 Filtered to {len(membrane_df)} membrane proteins ({len(membrane_df)/len(df)*100:.1f}%)")
    
    return membrane_df


def add_functional_categories(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add functional category columns to DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with protein data including 'Function [CC]' and 'Protein names' columns
    
    Returns
    -------
    pd.DataFrame
        DataFrame with added 'Functional_Categories' column
    """
    df = df.copy()
    
    df['Functional_Categories'] = df.apply(
        lambda row: classify_protein_function(
            row.get('Function [CC]', ''),
            row.get('Protein names', '')
        ),
        axis=1
    )
    
    # Add binary columns for each category
    all_categories = set()
    for cats in df['Functional_Categories']:
        all_categories.update(cats)
    
    for category in sorted(all_categories):
        col_name = f'Is_{category.replace(" ", "_")}'
        df[col_name] = df['Functional_Categories'].apply(lambda x: category in x)
    
    print(f"✅ Added functional categories")
    print(f"   Categories found: {sorted(all_categories)}")
    
    return df
