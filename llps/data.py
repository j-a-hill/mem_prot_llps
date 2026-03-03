"""
Section 1 - Data Loading and Classification.

Functions for loading LLPS protein data and classifying proteins by pLLPS score.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Tuple


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
            Path(__file__).parent.parent / "Human Phase separation data.xlsx",
            Path(__file__).parent.parent / "data" / "sample_data.xlsx",
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


def load_and_classify_data(filepath: str, high_threshold: float = 0.7, low_threshold: float = 0.4) -> pd.DataFrame:
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
