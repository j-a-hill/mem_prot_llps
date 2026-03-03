"""Generate all group networks programmatically"""
import pandas as pd
import networkx as nx
import numpy as np
import json
from pathlib import Path

# Setup - use absolute paths relative to project root
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
OUTPUT_DIR = PROJECT_ROOT / 'results/functional_group_networks'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load data
df_full = pd.read_csv(PROJECT_ROOT / 'results/full_dataset.csv')
df_groups = pd.read_csv(PROJECT_ROOT / 'results/functional_groups_with_pllps.csv')

# Groups to process (align with available groups)
GROUPS = [
    'Structural',
    'Ion Channel',
    'Transporter',
    'GPCR',
    'Receptor Tyrosine Kinase',
    'Ligase',
    'Receptor',
    'Nuclear Receptor',
    'Synthetase',
    'Guanine Nucleotide Exchange Factor'
]

for TARGET_GROUP in GROUPS:
    print(f"\n{'='*70}")
    print(f"PROCESSING: {TARGET_GROUP}")
    print(f"{'='*70}\n")
    
    # Load group data
    output_name = TARGET_GROUP.lower().replace(' ', '_')
    string_file = PROJECT_ROOT / 'results/string_networks_by_group' / f"{output_name}_interactions.csv"
    
    if not string_file.exists():
        print(f"⚠️  Skipping {TARGET_GROUP} - no interaction file found")
        continue
    
    interactions_df = pd.read_csv(string_file)
    print(f"✓ Loaded {len(interactions_df)} interactions")
    
    # Subset proteins for this group
    group_proteins = df_groups[df_groups['Functional Group'] == TARGET_GROUP].copy()
    protein_ids = group_proteins['Entry'].tolist()
    
    # Build graph
    G = nx.Graph()
    pllps_dict = df_full.set_index('Entry')['p(LLPS)'].to_dict()
    
    for _, row in interactions_df.iterrows():
        protein = row['protein']
        partner = row['partner']
        score = row['combined_score']
        score_norm = row['score_normalized']
        
        G.add_edge(protein, partner, 
                   weight=score_norm,
                   string_score=score,
                   string_score_norm=score_norm)
        
        if protein not in G.nodes:
            G.add_node(protein)
        if partner not in G.nodes:
            G.add_node(partner)
    
    # Add node attributes
    for node in G.nodes():
        pllps = pllps_dict.get(node, np.nan)
        G.nodes[node]['pllps'] = pllps
        G.nodes[node]['protein_name'] = node
        G.nodes[node]['uniprot_id'] = node
    
    # Add isolated group proteins not present in interactions
    added_isolated = 0
    for prot_id in protein_ids:
        if prot_id not in G.nodes:
            pllps_val = pllps_dict.get(prot_id, np.nan)
            prot_name = group_proteins.loc[group_proteins['Entry'] == prot_id, 'Entry name']
            prot_name = prot_name.values[0] if len(prot_name) > 0 else prot_id
            G.add_node(prot_id, pllps=pllps_val, protein_name=prot_name, uniprot_id=prot_id)
            added_isolated += 1
    if added_isolated:
        print(f"✓ Added {added_isolated} isolated nodes (no STRING edges at cutoff)")
    
    print(f"✓ Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    # Build gene name map
    gene_name_map = {}
    for _, row in interactions_df.iterrows():
        gene_name_map[row['protein']] = row['protein_name']
        gene_name_map[row['partner']] = row['partner_name']
    
    # Export nodes
    nodes_data = []
    for node in G.nodes():
        nodes_data.append({
            'uniprot_id': node,
            'gene_name': gene_name_map.get(node, node),
            'pllps_score': G.nodes[node]['pllps'],
            'degree': G.degree(node),
            'functional_group': TARGET_GROUP
        })
    
    nodes_df = pd.DataFrame(nodes_data)
    nodes_file = OUTPUT_DIR / f"{output_name}_nodes.csv"
    nodes_df.to_csv(nodes_file, index=False)
    print(f"✓ Saved {len(nodes_df)} nodes: {nodes_file}")
    
    # Export edges
    edges_data = []
    for u, v in G.edges():
        edges_data.append({
            'protein1_uniprot': u,
            'protein2_uniprot': v,
            'protein1_gene': gene_name_map.get(u, u),
            'protein2_gene': gene_name_map.get(v, v),
            'string_score': G[u][v]['string_score'],
            'string_score_normalized': G[u][v]['string_score_norm'],
            'pllps_1': G.nodes[u]['pllps'],
            'pllps_2': G.nodes[v]['pllps'],
            'pllps_mean': np.nanmean([G.nodes[u]['pllps'], G.nodes[v]['pllps']])
        })
    
    edges_df = pd.DataFrame(edges_data)
    edges_file = OUTPUT_DIR / f"{output_name}_edges.csv"
    edges_df.to_csv(edges_file, index=False)
    print(f"✓ Saved {len(edges_df)} edges: {edges_file}")
    
    # Export Cytoscape JSON (proper format for Cytoscape Desktop)
    cytoscape_data = {
        "data": {
            "name": f"{TARGET_GROUP} Network",
            "description": "STRING physical interactions colored by pLLPS score"
        },
        "elements": {
            "nodes": [],
            "edges": []
        }
    }
    
    for node in G.nodes():
        cytoscape_data["elements"]["nodes"].append({
            "data": {
                "id": node,
                "name": gene_name_map.get(node, node),
                "gene_name": gene_name_map.get(node, node),
                "pllps": float(G.nodes[node]['pllps']) if not np.isnan(G.nodes[node]['pllps']) else 0.0,
                "degree": G.degree(node),
                "uniprot_id": node
            }
        })
    
    for u, v in G.edges():
        cytoscape_data["elements"]["edges"].append({
            "data": {
                "id": f"{u}_{v}",
                "source": u,
                "target": v,
                "interaction": "interacts_with",
                "string_score": float(G[u][v]['string_score']),
                "string_score_norm": float(G[u][v]['string_score_norm']),
                "weight": float(G[u][v]['weight'])
            }
        })
    
    cytoscape_file = OUTPUT_DIR / f"{output_name}_network.json"
    with open(cytoscape_file, 'w') as f:
        json.dump(cytoscape_data, f, indent=2)
    print(f"✓ Saved Cytoscape format: {cytoscape_file}")
    
    # Export summary stats
    summary_stats = {
        'functional_group': TARGET_GROUP,
        'num_nodes': G.number_of_nodes(),
        'num_edges': G.number_of_edges(),
        'num_connected_components': nx.number_connected_components(G),
        'density': float(nx.density(G)),
        'mean_degree': float(np.mean([G.degree(n) for n in G.nodes()])),
        'mean_pllps': float(nodes_df['pllps_score'].mean()),
        'high_pllps_proteins': int((nodes_df['pllps_score'] > 0.7).sum()),
        'mean_string_score': float(edges_df['string_score'].mean()),
        'edge_confidence_cutoff': 600,
        'network_type': 'physical'
    }
    
    stats_file = OUTPUT_DIR / f"{output_name}_summary_stats.json"
    with open(stats_file, 'w') as f:
        json.dump(summary_stats, f, indent=2)
    print(f"✓ Saved summary: {stats_file}")
    
    print(f"\n✅ {TARGET_GROUP} complete!")

print(f"\n{'='*70}")
print("✅ ALL GROUPS EXPORTED SUCCESSFULLY")
print(f"{'='*70}")
