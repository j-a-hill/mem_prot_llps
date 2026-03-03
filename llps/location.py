"""
Section 2 - Location Parsing and Analysis.

Functions for parsing UniProt subcellular location strings and analysing
interaction preferences broken down by subcellular location.
"""

import re
import pandas as pd
from typing import List, Dict, Any, Union


def parse_location(location_str: Union[str, float, None]) -> List[str]:
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


def add_location_columns(df: pd.DataFrame) -> pd.DataFrame:
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


def analyze_interactions_by_location(matched_df: pd.DataFrame, full_dataset_df: pd.DataFrame, locations: List[str], high_threshold: float = 0.7, low_threshold: float = 0.4) -> Dict[str, Any]:
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
    # Import here to avoid circular dependency
    from llps.enrichment import analyze_interaction_matrix

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
