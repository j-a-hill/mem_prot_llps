"""
Section 3 - STRING Database Interactions.

Functions for querying the STRING protein interaction database.
"""

import requests
import time
import json
import warnings
import pandas as pd
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Callable

from llps.constants import (
    STRING_API_GET_IDS,
    STRING_API_NETWORK,
    STRING_API_PARTNERS,
    StringQueryConfig,
)


def get_string_mapping(protein_ids: List[str], species: int = 9606, batch_size: int = 2000) -> Tuple[Dict[str, str], Dict[str, str]]:
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
    url = STRING_API_GET_IDS
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
        except requests.Timeout:
            warnings.warn(f"  Batch {batch_num}/{total_batches}: Request timed out")
        except requests.RequestException as e:
            warnings.warn(f"  Batch {batch_num}/{total_batches}: Network error - {e}")
            
        time.sleep(1)
        
    return string_to_uniprot, uniprot_to_string


def fetch_string_interactions(
    protein_ids: List[str],
    config: Optional[StringQueryConfig] = None,
    progress_callback: Optional[Callable[[str], None]] = None,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Fetch protein-protein interactions from STRING database with optional caching.

    Parameters
    ----------
    protein_ids : List[str]
        List of UniProt protein IDs to fetch interactions for
    config : StringQueryConfig, optional
        Query configuration object. Uses default settings when not provided.
    progress_callback : Callable, optional
        Function to call with progress updates

    Returns
    -------
    Tuple[pd.DataFrame, List[str]]
        DataFrame with interaction data and list of error messages

    Examples
    --------
    >>> protein_ids = ['P04637', 'P38398', 'P51587']  # p53, BRCA1, BRCA2
    >>> cfg = StringQueryConfig(score_threshold=700)
    >>> interactions_df, errors = fetch_string_interactions(protein_ids, config=cfg)
    >>> print(f"Found {len(interactions_df)} interactions")
    """
    if config is None:
        config = StringQueryConfig()
    species = config.species
    score_threshold = config.score_threshold
    batch_size = config.batch_size
    use_cache = config.use_cache
    network_type = config.network_type
    string_api_url = STRING_API_NETWORK
    all_interactions = []
    errors = []
    
    # Try to load from cache first
    if use_cache:
        cache_file = Path(__file__).parent.parent / "data" / f"string_cache_{score_threshold}.json"
        if cache_file.exists():
            try:
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


def fetch_interaction_partners(
    protein_ids: List[str],
    uniprot_map: Optional[Dict[str, str]] = None,
    config: Optional[StringQueryConfig] = None,
) -> pd.DataFrame:
    """
    Fetch interaction partners for a list of proteins from STRING.

    Parameters
    ----------
    protein_ids : list
        List of UniProt IDs
    uniprot_map : dict, optional
        Dictionary mapping Uniprot IDs to STRING IDs (optional but recommended)
    config : StringQueryConfig, optional
        Query configuration object. Uses default settings when not provided.

    Returns
    -------
    pd.DataFrame
        DataFrame with interactions
    """
    if config is None:
        config = StringQueryConfig()
    species = config.species
    score_threshold = config.score_threshold
    batch_size = config.batch_size
    # Use interaction_partners API to get ALL partners, not just interactions within the set
    string_api_url = STRING_API_PARTNERS
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
        except requests.Timeout:
            warnings.warn(f"  Batch {batch_num}/{total_batches}: Request timed out")
        except requests.RequestException as e:
            warnings.warn(f"  Batch {batch_num}/{total_batches}: Network error - {e}")
        
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
