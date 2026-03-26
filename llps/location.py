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


def categorize_location_to_compartment(
    location_str: Union[str, float, None],
    is_membrane: bool = False
) -> str:
    """
    Categorize a UniProt subcellular location string into major cellular compartments.

    Uses the is_membrane flag to disambiguate locations like "Endoplasmic reticulum"
    which could mean the membrane or the lumen. This avoids overcounting membrane proteins.

    Parameters
    ----------
    location_str : str
        Subcellular location string from UniProt
    is_membrane : bool, optional
        Whether the protein is classified as a membrane protein. Default is False.
        If True and location contains membrane-associated keywords, categorizes accordingly.

    Returns
    -------
    str
        Major compartment category (e.g., 'Plasma Membrane', 'Cytosol', 'ER Lumen', 'Mitochondrion')

    Notes
    -----
    - Compartments ending in "membrane" are only classified as membrane if is_membrane=True
    - Proteins in ER/Golgi/etc without the membrane flag are categorized as "_Lumen"
    - This prevents overcounting of membrane proteins
    """
    if pd.isna(location_str) or location_str == '':
        return 'Unknown'

    location_str = str(location_str).lower()

    # Check for cytosolic/nuclear proteins first (typically non-membrane)
    if any(term in location_str for term in ['cytoplasm', 'cytosol', 'nucleus', 'nucleoplasm']):
        return 'Cytosol'

    # Plasma/cell membrane - these are typically membrane proteins
    if any(term in location_str for term in ['plasma membrane', 'cell membrane']):
        return 'Plasma Membrane'

    # Mitochondrion - check if membrane or matrix/lumen
    if 'mitochondri' in location_str:
        if is_membrane:
            return 'Mitochondrial Membrane'
        else:
            return 'Mitochondrial Matrix'

    # Endoplasmic reticulum - disambiguate with is_membrane flag
    if 'endoplasmic reticulum' in location_str or 'er ' in location_str or location_str.endswith('er'):
        if is_membrane:
            return 'ER Membrane'
        else:
            return 'ER Lumen'

    # Golgi apparatus - disambiguate with is_membrane flag
    if 'golgi' in location_str:
        if is_membrane:
            return 'Golgi Membrane'
        else:
            return 'Golgi Lumen'

    # Peroxisome
    if 'peroxisom' in location_str:
        if is_membrane:
            return 'Peroxisomal Membrane'
        else:
            return 'Peroxisomal Matrix'

    # Lysosome/Vacuole
    if any(term in location_str for term in ['lysosom', 'vacuole']):
        if is_membrane:
            return 'Lysosomal Membrane'
        else:
            return 'Lysosomal Lumen'

    # General membrane (but not plasma) - only if is_membrane=True
    if 'membrane' in location_str:
        if is_membrane:
            return 'Other Membrane'
        else:
            return 'Membrane-Associated Lumen'

    # If nothing matched but is_membrane is True, probably a membrane protein
    if is_membrane:
        return 'Other Membrane'

    return 'Other'


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
