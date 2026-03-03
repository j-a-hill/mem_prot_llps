"""
Section 4 - Network Analysis.

Functions for building and analysing protein interaction networks with
pLLPS annotations using NetworkX.
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional, Any

try:
    import networkx as nx
    _HAS_NETWORKX = True
except ImportError:
    _HAS_NETWORKX = False


def match_interactions_to_pllps(
    interactions_df: pd.DataFrame,
    pllps_df: pd.DataFrame,
    id_column: str = 'Entry'
) -> pd.DataFrame:
    """
    Match STRING interactions to pLLPS dataset and add pLLPS scores.
    
    STRING returns gene names (preferredName), but pLLPS dataset uses UniProt IDs.
    This function maps gene names to UniProt IDs using the Entry name column.
    
    Parameters
    ----------
    interactions_df : pd.DataFrame
        DataFrame from STRING with 'preferredName_A' and 'preferredName_B' columns
    pllps_df : pd.DataFrame
        DataFrame with protein data including pLLPS scores
    id_column : str, optional
        Column name in pllps_df containing protein IDs (default: 'Entry')
    
    Returns
    -------
    pd.DataFrame
        Interactions with added pLLPS scores for both partners and UniProt IDs
    
    Examples
    --------
    >>> matched_df = match_interactions_to_pllps(interactions_df, pllps_df)
    >>> print(f"Matched {len(matched_df)} interactions to pLLPS data")
    """
    if len(interactions_df) == 0:
        return pd.DataFrame()
    
    # Extract gene names from Entry name column (format: GENENAME_ORGANISM)
    pllps_df = pllps_df.copy()
    pllps_df['Gene'] = pllps_df['Entry name'].str.split('_').str[0]
    
    # Create lookup dictionaries mapping gene name to UniProt ID and pLLPS score
    gene_to_uniprot = dict(zip(pllps_df['Gene'], pllps_df[id_column]))
    gene_to_pllps = dict(zip(pllps_df['Gene'], pllps_df['p(LLPS)']))
    
    # Also create direct UniProt lookups (in case some are already mapped)
    uniprot_to_pllps = dict(zip(pllps_df[id_column], pllps_df['p(LLPS)']))
    
    # Match interactions to pLLPS data
    matched = interactions_df.copy()
    
    # Map gene names to UniProt IDs
    matched['protein1'] = matched['preferredName_A'].map(gene_to_uniprot)
    matched['protein2'] = matched['preferredName_B'].map(gene_to_uniprot)
    
    # If not found by gene name, try direct UniProt lookup
    matched['protein1'] = matched['protein1'].fillna(matched['preferredName_A'])
    matched['protein2'] = matched['protein2'].fillna(matched['preferredName_B'])
    
    # Map to pLLPS scores (try gene name first, then UniProt ID)
    matched['pllps_1'] = matched['preferredName_A'].map(gene_to_pllps).fillna(
        matched['protein1'].map(uniprot_to_pllps)
    )
    matched['pllps_2'] = matched['preferredName_B'].map(gene_to_pllps).fillna(
        matched['protein2'].map(uniprot_to_pllps)
    )
    
    # Add flag for whether both partners are in dataset
    matched['both_in_dataset'] = matched['pllps_1'].notna() & matched['pllps_2'].notna()
    
    # Keep only interactions where both partners are in the dataset
    matched = matched[matched['both_in_dataset']].copy()
    
    return matched


def match_interactors_to_pllps(interactions_df: pd.DataFrame, pllps_df: pd.DataFrame, string_map: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """
    Match interaction partners to the pLLPS dataset.
    
    Returns DataFrame with both proteins' pLLPS values.
    
    Parameters
    ----------
    interactions_df : pd.DataFrame
        DataFrame with interaction data from STRING
    pllps_df : pd.DataFrame
        DataFrame with pLLPS scores
    string_map : dict, optional
        Dictionary mapping STRING IDs to UniProt IDs
    
    Returns
    -------
    pd.DataFrame
        DataFrame with matched pLLPS scores
    """
    if len(interactions_df) == 0:
        return pd.DataFrame()
    
    # Create lookup dictionary: Entry -> pLLPS
    pllps_lookup = dict(zip(pllps_df['Entry'], pllps_df['p(LLPS)']))
    
    # Also create lookup by Entry name (STRING often uses gene names)
    entry_name_lookup = dict(zip(pllps_df['Entry name'], pllps_df['p(LLPS)']))
    
    results = []
    
    for _, row in interactions_df.iterrows():
        protein_a_name = row.get('preferredName_A', '')
        protein_a_id = row.get('stringId_A', '')
        
        protein_b_name = row.get('preferredName_B', '')
        protein_b_id = row.get('stringId_B', '')
        
        score = row.get('score', 0)
        
        # Resolve Protein A to Uniprot ID
        uniprot_a = None
        if string_map:
            uniprot_a = string_map.get(protein_a_name) or string_map.get(protein_a_id)
        
        # Fallback if no map or not found
        if not uniprot_a:
            if protein_a_name in pllps_lookup:
                uniprot_a = protein_a_name
            elif protein_a_name in entry_name_lookup:
                uniprot_a = protein_a_name
        
        # Resolve Protein B to Uniprot ID
        uniprot_b = None
        if string_map:
            uniprot_b = string_map.get(protein_b_name) or string_map.get(protein_b_id)
            
        if not uniprot_b:
            if protein_b_name in pllps_lookup:
                uniprot_b = protein_b_name

        # Get pLLPS scores
        pllps_a = pllps_lookup.get(uniprot_a)
        if pllps_a is None and not string_map:
            pllps_a = pllps_lookup.get(protein_a_name) or entry_name_lookup.get(protein_a_name)

        pllps_b = pllps_lookup.get(uniprot_b)
        if pllps_b is None and not string_map:
            pllps_b = pllps_lookup.get(protein_b_name) or entry_name_lookup.get(protein_b_name)
        
        results.append({
            'protein_a': protein_a_name,
            'protein_b': protein_b_name,
            'uniprot_a': uniprot_a,
            'uniprot_b': uniprot_b,
            'score': score,
            'pllps_a': pllps_a,
            'pllps_b': pllps_b
        })
    
    return pd.DataFrame(results)


def _build_pllps_graph(interactions_df: pd.DataFrame, pllps_df: pd.DataFrame) -> 'nx.Graph':
    """Build a NetworkX graph from interaction data with pLLPS node attributes."""
    G = nx.Graph()

    if 'protein1' in interactions_df.columns and 'protein2' in interactions_df.columns:
        col_a, col_b = 'protein1', 'protein2'
    elif 'preferredName_A' in interactions_df.columns:
        col_a, col_b = 'preferredName_A', 'preferredName_B'
    elif 'stringId_A' in interactions_df.columns:
        col_a, col_b = 'stringId_A', 'stringId_B'
    else:
        str_cols = interactions_df.select_dtypes(include=['object']).columns
        col_a, col_b = str_cols[0], str_cols[1]
        print(f"Using columns: {col_a}, {col_b}")

    if 'score' in interactions_df.columns:
        score_col: Optional[str] = 'score'
    elif 'combined_score' in interactions_df.columns:
        score_col = 'combined_score'
    else:
        score_col = None
        print("⚠️  No score column found in interactions data. Edges will have no weight.")

    for _, row in interactions_df.iterrows():
        edge_score = row[score_col] if score_col and pd.notna(row.get(score_col)) else None
        G.add_edge(row[col_a], row[col_b], score=edge_score)

    pllps_dict = dict(zip(pllps_df['Entry'], pllps_df['p(LLPS)']))
    if 'Entry name' in pllps_df.columns:
        pllps_dict.update(dict(zip(pllps_df['Entry name'], pllps_df['p(LLPS)'])))

    for node in G.nodes():
        pllps_value = pllps_dict.get(node)
        if pllps_value is None and '.' in str(node):
            pllps_value = pllps_dict.get(str(node).split('.')[-1])
        G.nodes[node]['pLLPS'] = pllps_value

    return G


def _classify_network_nodes(G: 'nx.Graph', high_threshold: float) -> Tuple[List, List, List]:
    """Classify graph nodes into high/low/unknown pLLPS categories."""
    high_pllps_nodes: List = []
    low_pllps_nodes: List = []
    unknown_nodes: List = []
    for node in G.nodes():
        pllps = G.nodes[node].get('pLLPS')
        if pllps is None:
            unknown_nodes.append(node)
        elif pllps >= high_threshold:
            high_pllps_nodes.append(node)
        else:
            low_pllps_nodes.append(node)
    return high_pllps_nodes, low_pllps_nodes, unknown_nodes


def _compute_network_metrics(
    G: 'nx.Graph',
    high_pllps_nodes: List,
    low_pllps_nodes: List,
    unknown_nodes: List,
) -> Dict[str, Any]:
    """Compute network topology and enrichment metrics."""
    results: Dict[str, Any] = {
        'total_nodes': G.number_of_nodes(),
        'total_edges': G.number_of_edges(),
        'density': nx.density(G),
        'avg_clustering': nx.average_clustering(G) if G.number_of_nodes() > 0 else 0,
        'high_pllps_nodes': len(high_pllps_nodes),
        'low_pllps_nodes': len(low_pllps_nodes),
        'unknown_pllps_nodes': len(unknown_nodes),
    }

    if len(high_pllps_nodes) > 1:
        high_subgraph = G.subgraph(high_pllps_nodes)
        high_pllps_edges = high_subgraph.number_of_edges()
        high_pllps_density = nx.density(high_subgraph)
        high_pllps_avg_clustering = nx.average_clustering(high_subgraph)
    else:
        high_pllps_edges = 0
        high_pllps_density = 0.0
        high_pllps_avg_clustering = 0.0
    results['high_pllps_edges'] = high_pllps_edges
    results['high_pllps_density'] = high_pllps_density
    results['high_pllps_avg_clustering'] = high_pllps_avg_clustering

    high_set = set(high_pllps_nodes)
    low_set = set(low_pllps_nodes)
    high_high = sum(1 for a, b in G.edges() if a in high_set and b in high_set)
    high_low = sum(1 for a, b in G.edges() if (a in high_set) != (b in high_set))
    low_low = sum(1 for a, b in G.edges() if a in low_set and b in low_set)
    results['high_high_interactions'] = high_high
    results['high_low_interactions'] = high_low
    results['low_low_interactions'] = low_low

    n_classified = len(high_pllps_nodes) + len(low_pllps_nodes)
    enrichment_ratio = 0.0
    if n_classified > 0:
        p_high = len(high_pllps_nodes) / n_classified
        expected_high_high_ratio = p_high * p_high
        total_classified_edges = high_high + high_low + low_low
        if total_classified_edges > 0 and expected_high_high_ratio > 0:
            enrichment_ratio = (high_high / total_classified_edges) / expected_high_high_ratio
    results['enrichment_ratio'] = enrichment_ratio

    degrees = dict(G.degree()) if G.number_of_nodes() > 0 else {}
    high_degrees = [degrees[n] for n in high_pllps_nodes if n in degrees]
    low_degrees = [degrees[n] for n in low_pllps_nodes if n in degrees]
    results['avg_degree_high_pllps'] = float(np.mean(high_degrees)) if high_degrees else 0.0
    results['avg_degree_low_pllps'] = float(np.mean(low_degrees)) if low_degrees else 0.0

    return results


def analyze_network(
    interactions_df: pd.DataFrame,
    pllps_df: pd.DataFrame,
    high_threshold: float = 0.7
) -> Tuple[Dict, 'nx.Graph']:
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
    Tuple[Dict, nx.Graph]
        Tuple containing:
        - Dict: Dictionary with network analysis results
        - nx.Graph: NetworkX graph object with pLLPS values as node attributes.
    """
    if not _HAS_NETWORKX:
        raise ImportError("networkx library required. Install with: pip install networkx")

    G = _build_pllps_graph(interactions_df, pllps_df)
    high_pllps_nodes, low_pllps_nodes, unknown_nodes = _classify_network_nodes(G, high_threshold)
    results = _compute_network_metrics(G, high_pllps_nodes, low_pllps_nodes, unknown_nodes)
    return results, G
