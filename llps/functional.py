"""
Section 9 - Functional Classification Utilities.

Functions for classifying proteins by functional category using YAML-defined
term patterns, and for identifying membrane proteins.
"""

import re
import yaml
import pandas as pd
from pathlib import Path
from typing import List, Dict


# Cache for loaded functional terms
_FUNCTIONAL_TERMS_CACHE = None


def _load_functional_terms(yaml_path: str = None) -> Dict:
    """
    Load functional classification terms from YAML file.
    
    Parameters
    ----------
    yaml_path : str, optional
        Path to YAML file. If None, uses default location.
    
    Returns
    -------
    Dict
        Dictionary with functional_groups and membrane_indicators
    """
    global _FUNCTIONAL_TERMS_CACHE
    
    if _FUNCTIONAL_TERMS_CACHE is not None:
        return _FUNCTIONAL_TERMS_CACHE
    
    if yaml_path is None:
        # Try to find the default YAML file
        possible_paths = [
            Path('data/functional_classification_terms.yaml'),
            Path(__file__).parent.parent / 'data' / 'functional_classification_terms.yaml',
        ]
        
        for path in possible_paths:
            if path.exists():
                yaml_path = path
                break
        
        if yaml_path is None:
            raise FileNotFoundError(
                "Could not find functional_classification_terms.yaml. "
                "Please provide the path explicitly."
            )
    
    with open(yaml_path, 'r') as f:
        terms = yaml.safe_load(f)
    
    _FUNCTIONAL_TERMS_CACHE = terms
    return terms


def parse_function_categories(
    function_str: str,
    protein_name_str: str = None,
    yaml_path: str = None
) -> List[str]:
    """
    Parse function and protein name strings to extract functional categories.
    
    Similar to parse_location(), this function:
    1. Cleans up the input text
    2. Matches against patterns defined in YAML file
    3. Returns list of matched functional categories
    
    Parameters
    ----------
    function_str : str
        Function [CC] annotation from UniProt
    protein_name_str : str, optional
        Protein names field from UniProt
    yaml_path : str, optional
        Path to YAML file with classification terms
    
    Returns
    -------
    List[str]
        List of functional categories that match the protein
    
    Examples
    --------
    >>> categories = parse_function_categories(
    ...     "Voltage-gated potassium channel",
    ...     "KCNA1_HUMAN Potassium voltage-gated channel"
    ... )
    >>> print(categories)
    ['Ion Channel']
    """
    if pd.isna(function_str) or function_str == '':
        function_str = ''
    if pd.isna(protein_name_str) or protein_name_str == '':
        protein_name_str = ''
    
    # Convert to strings and combine
    function_str = str(function_str)
    protein_name_str = str(protein_name_str)
    
    # Remove curly bracket content like {ECO:0000269|PubMed:12345}
    function_str = re.sub(r'\{[^}]*\}', '', function_str)
    protein_name_str = re.sub(r'\{[^}]*\}', '', protein_name_str)
    
    # Remove square bracket content
    function_str = re.sub(r'\[[^\]]*\]', '', function_str)
    protein_name_str = re.sub(r'\[[^\]]*\]', '', protein_name_str)
    
    # Remove normal bracket content like (By similarity) or (Probable)
    function_str = re.sub(r'\([^)]*\)', '', function_str)
    protein_name_str = re.sub(r'\([^)]*\)', '', protein_name_str)
    
    # Combine and lowercase for matching
    combined_text = f"{function_str} {protein_name_str}".lower()
    
    # Load functional terms from YAML
    terms = _load_functional_terms(yaml_path)
    functional_groups = terms.get('functional_groups', {})
    
    # Match against patterns
    categories = []
    for category, group_data in functional_groups.items():
        patterns = group_data.get('patterns', [])
        for pattern in patterns:
            if re.search(pattern, combined_text, re.IGNORECASE):
                categories.append(category)
                break  # Only add category once
    
    return categories


def is_membrane_protein(
    function_str: str,
    protein_name_str: str = None,
    location_str: str = None,
    yaml_path: str = None
) -> bool:
    """
    Check if a protein is a membrane protein based on annotations.
    
    Parameters
    ----------
    function_str : str
        Function [CC] annotation
    protein_name_str : str, optional
        Protein names
    location_str : str, optional
        Subcellular location [CC] annotation
    yaml_path : str, optional
        Path to YAML file with classification terms
    
    Returns
    -------
    bool
        True if protein is classified as a membrane protein
    """
    if not isinstance(function_str, str):
        function_str = ""
    if not isinstance(protein_name_str, str):
        protein_name_str = ""
    if not isinstance(location_str, str):
        location_str = ""
    
    # Clean up text
    function_str = re.sub(r'\{[^}]*\}', '', function_str)
    protein_name_str = re.sub(r'\{[^}]*\}', '', protein_name_str)
    location_str = re.sub(r'\{[^}]*\}', '', location_str)
    
    combined_text = f"{function_str} {protein_name_str} {location_str}".lower()
    
    # Load membrane indicators from YAML
    terms = _load_functional_terms(yaml_path)
    membrane_indicators = terms.get('membrane_indicators', {}).get('patterns', [])
    
    for pattern in membrane_indicators:
        if re.search(pattern, combined_text, re.IGNORECASE):
            return True
    
    return False


def classify_protein_function(function_str: str, protein_name_str: str = None) -> List[str]:
    """
    Classify protein into functional categories based on annotations.
    
    DEPRECATED: Use parse_function_categories() instead.
    Kept for backward compatibility.
    
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
    return parse_function_categories(function_str, protein_name_str)


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
    
    Similar to add_location_columns(), this function:
    1. Parses function and protein name fields
    2. Adds a 'Functional_Categories' column with list of categories
    3. Optionally adds binary columns for each category
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with protein data including 'Function [CC]' and 'Protein names' columns
    
    Returns
    -------
    pd.DataFrame
        DataFrame with added 'Functional_Categories' column
    
    Examples
    --------
    >>> df_with_categories = add_functional_categories(df)
    >>> print(df_with_categories['Functional_Categories'].head())
    """
    if 'Function [CC]' not in df.columns:
        print("⚠️  Warning: 'Function [CC]' column not found")
        return df
    
    df = df.copy()
    
    # Add column with list of all parsed functional categories
    df['Functional_Categories'] = df.apply(
        lambda row: parse_function_categories(
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
    
    print(f"✅ Added functional categories to {len(df)} proteins")
    print(f"   Categories found: {sorted(all_categories)}")
    print(f"   Total category assignments: {sum(len(cats) for cats in df['Functional_Categories'])}")
    
    return df
