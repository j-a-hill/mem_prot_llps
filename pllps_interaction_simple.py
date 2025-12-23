"""
High pLLPS Protein Interaction Analysis

⚠️ DEPRECATION NOTICE:
This module has been consolidated into llps_functions.py.
Please update your imports to use:
    from llps_functions import load_llps_data, get_high_pllps_proteins, ...
    
This file is kept for backward compatibility but may be removed in future versions.

Simple workflow to analyze if high pLLPS proteins preferentially interact with each other.

Steps:
1. Filter high pLLPS entries
2. Fetch interaction partners from STRING
3. Match interactors to pLLPS dataset
4. Test for high-high vs high-low over-representation

Usage:
    python pllps_interaction_simple.py
"""

import warnings
warnings.warn(
    "pllps_interaction_simple.py is deprecated. Please use llps_functions.py instead.",
    DeprecationWarning,
    stacklevel=2
)

import pandas as pd
import numpy as np
import requests
import time
from scipy import stats
from pathlib import Path


def load_and_filter_high_pllps(filepath: str, threshold: float = 0.7) -> tuple:
    """
    Step 1: Load data and filter high pLLPS proteins.
    
    Args:
        filepath: Path to Excel file with protein data
        threshold: pLLPS threshold for "high" classification (default: 0.7)
    
    Returns:
        tuple: (full dataframe, list of high pLLPS Entry IDs sorted by pLLPS descending)
    
    Raises:
        FileNotFoundError: If the specified file does not exist
        ValueError: If required columns are missing
    """
    # Validate file exists
    if not Path(filepath).exists():
        raise FileNotFoundError(f"Data file not found: {filepath}")
    
    df = pd.read_excel(filepath, engine='openpyxl')
    
    # Validate required columns
    required_cols = ['Entry', 'p(LLPS)']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    print(f"Loaded {len(df)} proteins")
    print(f"p(LLPS) range: {df['p(LLPS)'].min():.3f} - {df['p(LLPS)'].max():.3f}")
    
    df['pLLPS_class'] = np.where(df['p(LLPS)'] >= threshold, 'high', 'low')
    
    # Get high pLLPS proteins sorted by pLLPS (highest first)
    high_pllps_df = df[df['pLLPS_class'] == 'high'].sort_values('p(LLPS)', ascending=False)
    high_pllps_ids = high_pllps_df['Entry'].tolist()
    
    print(f"\nHigh pLLPS (>= {threshold}): {len(high_pllps_ids)} proteins")
    print(f"Low pLLPS: {len(df) - len(high_pllps_ids)} proteins")
    
    return df, high_pllps_ids


def fetch_interactions(protein_ids: list, species: int = 9606, 
                       score_threshold: int = 700, batch_size: int = 100) -> pd.DataFrame:
    """
    Step 2: Fetch interaction partners from STRING database.
    
    Args:
        protein_ids: List of UniProt IDs
        species: NCBI taxonomy ID (9606 = human)
        score_threshold: Minimum confidence (0-1000), 700 = high confidence
        batch_size: Proteins per request
    
    Returns:
        DataFrame with interactions
    
    Note:
        - Requires network access to string-db.org
        - Implements rate limiting (1 request/second)
        - Handles common API errors with informative messages
    """
    string_api_url = "https://string-db.org/api/json/network"
    all_interactions = []
    
    total_batches = (len(protein_ids) + batch_size - 1) // batch_size
    print(f"\nFetching interactions for {len(protein_ids)} proteins ({total_batches} batches)...")
    
    for i in range(0, len(protein_ids), batch_size):
        batch = protein_ids[i:i+batch_size]
        batch_num = i // batch_size + 1
        
        params = {
            "identifiers": "\r".join(batch),
            "species": species,
            "required_score": score_threshold,
            "caller_identity": "pllps_analysis"
        }
        
        try:
            response = requests.post(string_api_url, data=params, timeout=60)
            if response.status_code == 200:
                interactions = response.json()
                all_interactions.extend(interactions)
                print(f"  Batch {batch_num}/{total_batches}: {len(interactions)} interactions")
            elif response.status_code == 429:
                print(f"  Batch {batch_num}: Rate limited. Waiting 30s...")
                time.sleep(30)
                # Retry once
                response = requests.post(string_api_url, data=params, timeout=60)
                if response.status_code == 200:
                    interactions = response.json()
                    all_interactions.extend(interactions)
                    print(f"  Batch {batch_num}/{total_batches}: {len(interactions)} interactions (retry)")
            elif response.status_code == 400:
                print(f"  Batch {batch_num}: Bad request - check protein IDs format")
            elif response.status_code == 404:
                print(f"  Batch {batch_num}: No interactions found for these proteins")
            elif response.status_code >= 500:
                print(f"  Batch {batch_num}: STRING server error ({response.status_code}). Try again later.")
            else:
                print(f"  Batch {batch_num}/{total_batches}: Error {response.status_code}")
        except Exception as e:
            print(f"  Batch {batch_num}/{total_batches}: {e}")
        
        time.sleep(1)
    
    print(f"\nTotal interactions: {len(all_interactions)}")
    return pd.DataFrame(all_interactions) if all_interactions else pd.DataFrame()


def match_to_pllps(interactions_df: pd.DataFrame, pllps_df: pd.DataFrame) -> pd.DataFrame:
    """
    Step 3: Match interaction partners to pLLPS dataset.
    
    Returns:
        DataFrame with both proteins' pLLPS values
    """
    if len(interactions_df) == 0:
        return pd.DataFrame()
    
    # Create lookup dictionaries
    pllps_by_entry = dict(zip(pllps_df['Entry'], pllps_df['p(LLPS)']))
    pllps_by_name = dict(zip(pllps_df['Entry name'], pllps_df['p(LLPS)']))
    
    results = []
    for _, row in interactions_df.iterrows():
        protein_a = row.get('preferredName_A', row.get('stringId_A', ''))
        protein_b = row.get('preferredName_B', row.get('stringId_B', ''))
        score = row.get('score', 0)
        
        # Try multiple matching strategies
        pllps_a = pllps_by_entry.get(protein_a) or pllps_by_name.get(protein_a)
        pllps_b = pllps_by_entry.get(protein_b) or pllps_by_name.get(protein_b)
        
        results.append({
            'protein_a': protein_a,
            'protein_b': protein_b,
            'score': score,
            'pllps_a': pllps_a,
            'pllps_b': pllps_b
        })
    
    matched = pd.DataFrame(results)
    
    # Report matching statistics
    both_matched = (matched['pllps_a'].notna() & matched['pllps_b'].notna()).sum()
    print(f"\nMatching results:")
    print(f"  Both proteins matched to pLLPS: {both_matched}/{len(matched)} ({100*both_matched/len(matched):.1f}%)")
    
    return matched


def test_enrichment(matched_df: pd.DataFrame, threshold: float = 0.7) -> dict:
    """
    Step 4: Test for high-high vs high-low over-representation.
    
    Performs chi-squared test to determine if high pLLPS proteins
    preferentially interact with each other.
    
    Returns:
        dict with enrichment results
    """
    # Keep only interactions where both proteins have pLLPS scores
    complete = matched_df.dropna(subset=['pllps_a', 'pllps_b']).copy()
    
    if len(complete) == 0:
        print("No interactions with complete pLLPS data.")
        return None
    
    # Classify interactions
    complete['class_a'] = np.where(complete['pllps_a'] >= threshold, 'high', 'low')
    complete['class_b'] = np.where(complete['pllps_b'] >= threshold, 'high', 'low')
    
    def interaction_type(row):
        if row['class_a'] == 'high' and row['class_b'] == 'high':
            return 'high-high'
        elif row['class_a'] == 'low' and row['class_b'] == 'low':
            return 'low-low'
        return 'high-low'
    
    complete['interaction_type'] = complete.apply(interaction_type, axis=1)
    
    # Count
    counts = complete['interaction_type'].value_counts()
    total = len(complete)
    high_high = counts.get('high-high', 0)
    high_low = counts.get('high-low', 0)
    low_low = counts.get('low-low', 0)
    
    # Calculate expected (null hypothesis)
    all_proteins = pd.concat([
        complete[['protein_a', 'pllps_a']].rename(columns={'protein_a': 'p', 'pllps_a': 'v'}),
        complete[['protein_b', 'pllps_b']].rename(columns={'protein_b': 'p', 'pllps_b': 'v'})
    ]).drop_duplicates(subset='p')
    
    n_high = (all_proteins['v'] >= threshold).sum()
    p_high = n_high / len(all_proteins)
    p_low = 1 - p_high
    
    expected_hh = p_high ** 2
    expected_hl = 2 * p_high * p_low
    expected_ll = p_low ** 2
    
    observed_hh = high_high / total
    enrichment = observed_hh / expected_hh if expected_hh > 0 else 0
    
    # Print results
    print("\n" + "="*60)
    print("HIGH pLLPS INTERACTION ENRICHMENT ANALYSIS")
    print("="*60)
    print(f"\nTotal interactions analyzed: {total}")
    print(f"\nObserved:")
    print(f"  High-High: {high_high:4d} ({100*high_high/total:5.1f}%)")
    print(f"  High-Low:  {high_low:4d} ({100*high_low/total:5.1f}%)")
    print(f"  Low-Low:   {low_low:4d} ({100*low_low/total:5.1f}%)")
    print(f"\nExpected by chance:")
    print(f"  High-High: {expected_hh*100:5.1f}%")
    print(f"  High-Low:  {expected_hl*100:5.1f}%")
    print(f"  Low-Low:   {expected_ll*100:5.1f}%")
    print(f"\n{'='*60}")
    print(f"ENRICHMENT OF HIGH-HIGH INTERACTIONS: {enrichment:.2f}x")
    print("="*60)
    
    # Chi-squared test
    observed = [high_high, high_low, low_low]
    expected = [expected_hh * total, expected_hl * total, expected_ll * total]
    
    if all(e > 5 for e in expected):
        chi2, p_value = stats.chisquare(observed, expected)
        print(f"\nChi-squared test: χ² = {chi2:.2f}, p = {p_value:.2e}")
        
        if p_value < 0.05:
            if enrichment > 1:
                print("\n✅ SIGNIFICANT: High pLLPS proteins preferentially interact with each other!")
            else:
                print("\n✅ SIGNIFICANT: High pLLPS proteins tend to avoid each other.")
        else:
            print("\n⚠️ NOT SIGNIFICANT: No preference detected (p >= 0.05)")
    else:
        print("\n⚠️ Chi-squared test not valid (sample size too small)")
        p_value = None
        chi2 = None
    
    return {
        'high_high': high_high,
        'high_low': high_low,
        'low_low': low_low,
        'total': total,
        'enrichment': enrichment,
        'p_value': p_value,
        'chi2': chi2
    }


def main(data_file: str = None, threshold: float = 0.7, sample_size: int = 100):
    """
    Run the complete analysis pipeline.
    
    Args:
        data_file: Path to Excel file (default: 'Human Phase separation data.xlsx')
        threshold: pLLPS threshold (default: 0.7)
        sample_size: Number of top high pLLPS proteins to analyze (default: 100, None for all)
    """
    print("="*60)
    print("pLLPS PROTEIN INTERACTION ANALYSIS")
    print("="*60)
    
    # Use default if not specified
    if data_file is None:
        data_file = "Human Phase separation data.xlsx"
    
    # Step 1: Load and filter
    print("\n[Step 1] Loading data and filtering high pLLPS proteins...")
    df, high_pllps_ids = load_and_filter_high_pllps(data_file, threshold)
    
    # Step 2: Fetch interactions
    # Note: high_pllps_ids is already sorted by pLLPS descending
    print("\n[Step 2] Fetching interaction partners from STRING...")
    if sample_size:
        proteins_to_analyze = high_pllps_ids[:sample_size]
        print(f"(Analyzing top {sample_size} highest pLLPS proteins)")
    else:
        proteins_to_analyze = high_pllps_ids
    
    interactions = fetch_interactions(proteins_to_analyze)
    
    if len(interactions) == 0:
        print("\n❌ No interactions found. Check network connectivity.")
        return None
    
    # Step 3: Match to pLLPS data
    print("\n[Step 3] Matching interactors to pLLPS dataset...")
    matched = match_to_pllps(interactions, df)
    
    # Step 4: Test enrichment
    print("\n[Step 4] Testing for high-high enrichment...")
    results = test_enrichment(matched, threshold)
    
    if results:
        print("\n" + "="*60)
        print("ANALYSIS COMPLETE")
        print("="*60)
    
    return results


if __name__ == "__main__":
    main()
