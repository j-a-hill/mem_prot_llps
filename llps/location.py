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
