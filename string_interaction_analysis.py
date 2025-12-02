"""
LLPS Protein Interaction Analysis Module

This module provides functions for analyzing protein-protein interactions
among high pLLPS (Liquid-Liquid Phase Separation) proteins using the STRING
database and network analysis techniques.

Usage:
    python string_interaction_analysis.py

Or import as a module:
    from string_interaction_analysis import (
        load_llps_data,
        get_high_pllps_proteins,
        fetch_string_interactions,
        analyze_network
    )
"""

import pandas as pd
import numpy as np
from pathlib import Path
import time
from typing import List, Dict, Tuple, Optional
import warnings

# Configuration
DEFAULT_PLLPS_THRESHOLD = 0.7
DEFAULT_STRING_SCORE_THRESHOLD = 400
SPECIES_ID = 9606  # Human


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
    
    df = pd.read_excel(filepath, engine='openpyxl')
    print(f"✅ Loaded {len(df)} proteins from {filepath}")
    
    if 'p(LLPS)' in df.columns:
        print(f"   p(LLPS) range: {df['p(LLPS)'].min():.3f} - {df['p(LLPS)'].max():.3f}")
    
    return df


def get_high_pllps_proteins(
    df: pd.DataFrame,
    threshold: float = DEFAULT_PLLPS_THRESHOLD,
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
    
    high_pllps_ids = df[df['pLLPS_class'] == 'high']['Entry'].tolist()
    
    return df, high_pllps_ids


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


def fetch_string_interactions(
    protein_ids: List[str],
    species: int = SPECIES_ID,
    score_threshold: int = DEFAULT_STRING_SCORE_THRESHOLD,
    batch_size: int = 100,
    network_type: str = "physical"
) -> pd.DataFrame:
    """
    Fetch protein-protein interactions from STRING database.
    
    Note: This function requires network access to string-db.org.
    If network access is unavailable, use the STRING web interface
    and export_protein_list() function instead.
    
    Parameters
    ----------
    protein_ids : List[str]
        List of UniProt protein IDs.
    species : int
        NCBI taxonomy ID (9606 for human).
    score_threshold : int
        Minimum combined score (0-1000). 400=medium, 700=high, 900=highest.
    batch_size : int
        Number of proteins per API request.
    network_type : str
        "physical" for experimentally determined, "functional" for all.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with interaction data.
    """
    try:
        import requests
    except ImportError:
        raise ImportError("requests library required. Install with: pip install requests")
    
    string_api_url = "https://string-db.org/api/json/network"
    
    all_interactions = []
    total_batches = (len(protein_ids) + batch_size - 1) // batch_size
    
    print(f"\n🔄 Fetching STRING interactions ({total_batches} batches)...")
    
    for i in range(0, len(protein_ids), batch_size):
        batch = protein_ids[i:i+batch_size]
        batch_num = i // batch_size + 1
        
        params = {
            "identifiers": "\r".join(batch),
            "species": species,
            "required_score": score_threshold,
            "network_type": network_type,
            "caller_identity": "llps_analysis"
        }
        
        try:
            response = requests.post(string_api_url, data=params, timeout=60)
            if response.status_code == 200:
                interactions = response.json()
                all_interactions.extend(interactions)
                print(f"   Batch {batch_num}/{total_batches}: {len(interactions)} interactions")
            else:
                print(f"   Batch {batch_num}/{total_batches}: Error {response.status_code}")
        except Exception as e:
            print(f"   Batch {batch_num}/{total_batches}: Exception - {e}")
        
        time.sleep(1)  # Rate limiting
    
    if not all_interactions:
        warnings.warn("No interactions retrieved. Check network connectivity to string-db.org")
        return pd.DataFrame()
    
    df_interactions = pd.DataFrame(all_interactions)
    print(f"\n✅ Retrieved {len(df_interactions)} total interactions")
    
    return df_interactions


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


def analyze_network(
    interactions_df: pd.DataFrame,
    pllps_df: pd.DataFrame,
    high_threshold: float = DEFAULT_PLLPS_THRESHOLD
) -> Dict:
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
    Dict
        Dictionary with network analysis results.
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
    
    # Add edges
    score_col = 'score' if 'score' in interactions_df.columns else 'combined_score'
    
    for _, row in interactions_df.iterrows():
        G.add_edge(
            row[col_a],
            row[col_b],
            score=row.get(score_col, 0)
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


def print_analysis_report(results: Dict):
    """Print formatted analysis report."""
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


def main():
    """Main function demonstrating the analysis workflow."""
    print("="*60)
    print("🧬 LLPS Protein Interaction Analysis")
    print("="*60)
    
    # Load LLPS data
    try:
        df = load_llps_data()
    except FileNotFoundError as e:
        print(f"\n❌ {e}")
        print("Please provide the path to your LLPS data file.")
        return
    
    # Get high pLLPS proteins
    df, high_pllps_ids = get_high_pllps_proteins(df, threshold=0.7)
    
    # Export for STRING web interface
    export_protein_list(high_pllps_ids)
    
    print("\n" + "-"*60)
    print("📌 NEXT STEPS:")
    print("-"*60)
    print("""
    Since direct API access may be limited in some environments,
    follow these steps to complete the analysis:
    
    1. Go to https://string-db.org/cgi/input
    2. Select 'Multiple proteins' tab
    3. Paste the protein IDs from 'high_pllps_proteins.txt'
       (or upload the file directly)
    4. Select 'Homo sapiens (Human)' as the organism
    5. Set confidence threshold (suggest: 0.4 for medium, 0.7 for high)
    6. Click 'Continue'
    7. Download the network:
       - Click 'Exports' tab
       - Download as TSV (save as 'string_network.tsv')
    
    8. Then run the analysis:
    
       from string_interaction_analysis import (
           load_llps_data, load_string_network_file, analyze_network
       )
       
       df = load_llps_data()
       interactions = load_string_network_file('string_network.tsv')
       results, G = analyze_network(interactions, df)
       print_analysis_report(results)
    """)
    
    # If there's an existing network file, analyze it
    network_files = list(Path('.').glob('*string*.tsv')) + list(Path('.').glob('*STRING*.tsv'))
    
    if network_files:
        print(f"\n📂 Found existing network file: {network_files[0]}")
        print("   Running analysis...")
        
        interactions_df = load_string_network_file(str(network_files[0]))
        results, G = analyze_network(interactions_df, df)
        print_analysis_report(results)


if __name__ == "__main__":
    main()
