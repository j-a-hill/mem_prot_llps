"""
Section 1 - Data Loading and Classification.

Functions for loading LLPS protein data and classifying proteins by pLLPS score.
"""

import time
import warnings
from io import StringIO
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import requests


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


_UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


def fetch_uniprot_tm_annotations(
    entry_ids: List[str],
    batch_size: int = 100,
    cache_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Fetch Transmembrane and Intramembrane annotations from the UniProt REST API.

    Parameters
    ----------
    entry_ids : List[str]
        UniProt accession IDs to look up.
    batch_size : int
        Proteins per API request. Keep ≤ 100 to avoid URL length limits.
    cache_path : str, optional
        If given, load from this CSV path on first call and save there after fetching.

    Returns
    -------
    pd.DataFrame
        Columns: Entry, Transmembrane, Intramembrane
    """
    if cache_path and Path(cache_path).exists():
        print(f"Loading TM annotation cache: {cache_path}")
        return pd.read_csv(cache_path)

    total_batches = (len(entry_ids) + batch_size - 1) // batch_size
    print(f"Fetching TM annotations for {len(entry_ids)} proteins in {total_batches} batches...")

    records: List[pd.DataFrame] = []
    for i in range(0, len(entry_ids), batch_size):
        batch = entry_ids[i : i + batch_size]
        batch_num = i // batch_size + 1
        query = "accession:(" + " OR ".join(batch) + ")"
        params = {
            "query": query,
            "fields": "accession,ft_transmem,ft_intramem",
            "format": "tsv",
            "size": batch_size,
        }
        try:
            resp = requests.get(_UNIPROT_SEARCH_URL, params=params, timeout=30)
            if resp.status_code == 429:
                warnings.warn(f"Rate limited (batch {batch_num}), waiting 10s...")
                time.sleep(10)
                resp = requests.get(_UNIPROT_SEARCH_URL, params=params, timeout=30)
            if resp.status_code == 200:
                chunk = pd.read_csv(StringIO(resp.text), sep="\t")
                records.append(chunk)
                print(f"  Batch {batch_num}/{total_batches}: {len(chunk)} proteins")
            else:
                warnings.warn(f"  Batch {batch_num}/{total_batches}: HTTP {resp.status_code}")
        except requests.RequestException as e:
            warnings.warn(f"  Batch {batch_num}/{total_batches}: {e}")
        time.sleep(0.5)

    if not records:
        return pd.DataFrame(columns=["Entry", "Transmembrane", "Intramembrane"])

    result = pd.concat(records, ignore_index=True)

    # Normalise column names (UniProt API field headers can vary)
    col_map: dict = {}
    for col in result.columns:
        low = col.lower()
        if low == "entry":
            col_map[col] = "Entry"
        elif "intramembrane" in low:
            col_map[col] = "Intramembrane"
        elif "transmembrane" in low:
            col_map[col] = "Transmembrane"
    result = result.rename(columns=col_map)
    for col in ("Transmembrane", "Intramembrane"):
        if col not in result.columns:
            result[col] = np.nan

    result = result[["Entry", "Transmembrane", "Intramembrane"]]

    if cache_path:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(cache_path, index=False)
        print(f"Saved TM annotation cache to: {cache_path}")

    return result


def add_tmd_count(
    df: pd.DataFrame,
    tm_annotations: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    Add a TMD_count column counting transmembrane + intramembrane domain spans.

    If *tm_annotations* is provided (as returned by fetch_uniprot_tm_annotations),
    it is left-joined onto *df* on Entry before counting. Otherwise the function
    uses any Transmembrane/Intramembrane columns already present in *df*.

    Parameters
    ----------
    df : pd.DataFrame
        Protein DataFrame with an Entry column.
    tm_annotations : pd.DataFrame, optional
        DataFrame with Entry, Transmembrane, Intramembrane columns.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with added TMD_count column.
    """
    from llps.functional import count_tm_domains

    df = df.copy()

    if tm_annotations is not None:
        merge_cols = ["Entry"] + [
            c for c in ("Transmembrane", "Intramembrane") if c in tm_annotations.columns
        ]
        df = df.merge(tm_annotations[merge_cols], on="Entry", how="left")

    tmd = (
        df["Transmembrane"].apply(count_tm_domains)
        if "Transmembrane" in df.columns
        else pd.Series(0, index=df.index)
    )
    imd = (
        df["Intramembrane"].apply(count_tm_domains)
        if "Intramembrane" in df.columns
        else pd.Series(0, index=df.index)
    )
    df["TMD_count"] = tmd + imd

    n_with = (df["TMD_count"] > 0).sum()
    print(
        f"TMD_count added: {n_with} proteins have ≥1 membrane domain "
        f"(max={df['TMD_count'].max()})"
    )
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
