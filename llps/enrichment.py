"""
Section 5 - Enrichment Analysis.

Functions for analysing enrichment of interactions between protein classes.
"""

import numpy as np
import pandas as pd
from typing import Any, Dict, List, Optional
from scipy.stats import chi2_contingency
from llps.location import add_location_columns


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
    
    # Handle multiple column name conventions
    if 'pllps_1' in matched_df.columns and 'pllps_2' in matched_df.columns:
        # New convention
        matched_df['class_A'] = matched_df['pllps_1'] >= threshold
        matched_df['class_B'] = matched_df['pllps_2'] >= threshold
        pllps_a_col, pllps_b_col = 'pllps_1', 'pllps_2'
        name_a_col = 'protein1' if 'protein1' in matched_df.columns else 'preferredName_A'
        name_b_col = 'protein2' if 'protein2' in matched_df.columns else 'preferredName_B'
    elif 'pllps_A' in matched_df.columns:
        matched_df['class_A'] = matched_df['pllps_A'] >= threshold
        matched_df['class_B'] = matched_df['pllps_B'] >= threshold
        pllps_a_col, pllps_b_col = 'pllps_A', 'pllps_B'
        name_a_col, name_b_col = 'preferredName_A', 'preferredName_B'
    else:
        # Old convention
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


def analyze_interaction_matrix(matched_df: pd.DataFrame, full_dataset_df: pd.DataFrame, high_threshold: float = 0.7, low_threshold: float = 0.4) -> Optional[Dict]:
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
    # Handle different column name conventions
    if 'pllps_1' in matched_df.columns and 'pllps_2' in matched_df.columns:
        pllps_cols = ['pllps_1', 'pllps_2']
        pllps_a_col, pllps_b_col = 'pllps_1', 'pllps_2'
    else:
        pllps_cols = ['pllps_a', 'pllps_b']
        pllps_a_col, pllps_b_col = 'pllps_a', 'pllps_b'
    
    complete_df = matched_df.dropna(subset=pllps_cols).copy()
    
    if len(complete_df) == 0:
        print("No interactions with complete pLLPS data found.")
        return None
    
    # Define classes for the interaction pairs
    def get_class(score):
        if score >= high_threshold: return 'High'
        elif score >= low_threshold: return 'Medium'
        else: return 'Low'
        
    complete_df['class_a'] = complete_df[pllps_a_col].apply(get_class)
    complete_df['class_b'] = complete_df[pllps_b_col].apply(get_class)
    
    # Calculate Genomic Background Probabilities
    total_proteins = len(full_dataset_df)
    # Check for both column name variations
    pllps_class_col = 'pLLPS_Class' if 'pLLPS_Class' in full_dataset_df.columns else 'pLLPS_class'
    class_counts = full_dataset_df[pllps_class_col].value_counts()
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



def analyze_interactions_by_location(
    matched_df: pd.DataFrame,
    full_dataset_df: pd.DataFrame,
    locations: list[str],
    high_threshold: float = 0.7,
    low_threshold: float = 0.4,
) -> dict[str, Any]:
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

    if 'Location Categories' not in full_dataset_df.columns:
        full_dataset_df = add_location_columns(full_dataset_df)

    if 'Location Categories' not in matched_df.columns:
        loc_map = dict(zip(full_dataset_df['Entry'], full_dataset_df['Location Categories']))

        def get_locs(uid: str) -> list[str]:
            return loc_map.get(uid, [])

        matched_df['locs_a'] = matched_df['uniprot_a'].apply(get_locs)
        matched_df['locs_b'] = matched_df['uniprot_b'].apply(get_locs)

    for loc in locations:
        print(f"\nAnalyzing location: {loc}")

        loc_full_df = full_dataset_df[
            full_dataset_df['Location Categories'].apply(lambda x: loc in x)
        ]

        if len(loc_full_df) < 10:
            print(f"  Skipping {loc}: Too few proteins ({len(loc_full_df)})")
            continue

        loc_matched_df = matched_df[
            matched_df['locs_a'].apply(lambda x: loc in x)
            & matched_df['locs_b'].apply(lambda x: loc in x)
        ].copy()

        print(f"  Proteins in location: {len(loc_full_df)}")
        print(f"  Interactions within location: {len(loc_matched_df)}")

        class_counts = loc_full_df['pLLPS_class'].value_counts()
        print(
            f"  Class breakdown: High={class_counts.get('High', 0)}, "
            f"Medium={class_counts.get('Medium', 0)}, Low={class_counts.get('Low', 0)}"
        )

        if len(loc_matched_df) < 10:
            print("  Not enough interactions for analysis.")
            continue

        loc_results = analyze_interaction_matrix(loc_matched_df, loc_full_df, high_threshold, low_threshold)
        if loc_results:
            results[loc] = loc_results

    return results
