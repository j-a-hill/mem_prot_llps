# Exploring High pLLPS Protein Interactions: A Research Framework

## Overview

This document explores how proteins with high pLLPS (probability of Liquid-Liquid Phase Separation) interact with other proteins. The goal is to understand whether high pLLPS proteins:
1. Form biologically relevant clusters with other high pLLPS proteins
2. Act as central hubs in protein interaction networks
3. Exist in isolation from other phase-separating proteins

## Table of Contents
1. [Background](#background)
2. [Available APIs for Protein Interaction Data](#available-apis-for-protein-interaction-data)
3. [Analysis Strategies](#analysis-strategies)
4. [Implementation Approach](#implementation-approach)
5. [Sample Code](#sample-code)
6. [Expected Outcomes](#expected-outcomes)

---

## Background

### What is LLPS?
Liquid-Liquid Phase Separation (LLPS) is a process where proteins form condensates or droplets within cells. Proteins with high pLLPS scores are predicted to participate in phase separation, often forming:
- Stress granules
- P-bodies
- Nuclear speckles
- Nucleoli

### Key Questions
1. **Do high pLLPS proteins preferentially interact with each other?**
   - Hypothesis: Phase-separating proteins may co-localize and interact to form condensates
   
2. **Are high pLLPS proteins network hubs?**
   - Hypothesis: Phase-separating proteins might serve as scaffolds with many interaction partners
   
3. **Do high pLLPS proteins show functional clustering?**
   - Hypothesis: Interacting phase-separating proteins might share similar biological functions

---

## Available APIs for Protein Interaction Data

### 1. STRING MCP Server (Recommended for AI-Assisted Analysis)

**Repository:** https://github.com/Augmented-Nature/STRING-db-MCP-Server

**Description:** A Model Context Protocol (MCP) server that provides AI agents with direct access to STRING database functionality. This is the **recommended approach** for AI-assisted analysis as it integrates seamlessly with AI tools like Claude.

**Installation:**
```bash
npm install
npm run build
```

**Configuration (Claude Desktop):**
```json
{
  "mcpServers": {
    "string-server": {
      "command": "node",
      "args": ["/path/to/string-server/build/index.js"]
    }
  }
}
```

**Available Tools:**
| Tool | Description |
|------|-------------|
| `get_protein_interactions` | Get direct interaction partners with confidence scores |
| `get_interaction_network` | Build and analyze protein interaction networks |
| `get_functional_enrichment` | GO terms, KEGG pathways enrichment analysis |
| `get_protein_annotations` | Detailed protein annotations and functions |
| `find_homologs` | Find homologous proteins across species |
| `search_proteins` | Search proteins by name or identifier |

**Resource Templates:**
- `string://network/{protein_ids}` - Interaction network data
- `string://enrichment/{protein_ids}` - Functional enrichment results
- `string://interactions/{protein_id}` - Direct interaction partners
- `string://annotations/{protein_id}` - Protein annotations

**Example Usage with MCP:**
```
# Get interaction network for high pLLPS proteins
Use get_interaction_network with proteins: Q9Y6V0, P53350, Q9Y566, Q9Y520

# Perform functional enrichment on high pLLPS proteins
Use get_functional_enrichment with protein list and species 9606 (human)

# Find hub proteins
Use get_protein_interactions for each high pLLPS protein and count partners
```

**Advantages of MCP Approach:**
- ✅ Direct integration with AI assistants
- ✅ Natural language queries
- ✅ Automatic result interpretation
- ✅ No need to write API handling code
- ✅ Error handling built-in

---

### 2. STRING REST API (Direct Programmatic Access)

**Website:** https://string-db.org

**Description:** STRING is a comprehensive database of known and predicted protein-protein interactions, including:
- Physical interactions (experimentally determined)
- Functional associations (predicted)
- Text-mining derived associations

**API Base URL:** `https://string-db.org/api/`

**Key Endpoints:**

| Endpoint | Description |
|----------|-------------|
| `/api/json/get_string_ids` | Map protein identifiers to STRING IDs |
| `/api/json/network` | Get interaction network for proteins |
| `/api/json/interaction_partners` | Get interaction partners for proteins |
| `/api/json/enrichment` | Get functional enrichment |
| `/api/tsv/network` | Get network in TSV format |

**Confidence Scores:**
- Combined score (0-1000): Higher = more confident
- Recommended threshold: 400 (medium), 700 (high), 900 (highest)

**Rate Limits:**
- Free API access with reasonable rate limits
- For batch queries, they recommend max 200 proteins per request

### 3. BioGRID

**Website:** https://thebiogrid.org

**Description:** A curated database of protein, genetic, and chemical interactions.

**API Access:** REST API available with API key (free registration required)

**API Base URL:** `https://webservice.thebiogrid.org/`

**Key Features:**
- Focus on experimentally validated interactions
- Genetic interaction data (unique feature)
- Chemical-protein interactions

### 3. IntAct / PSICQUIC

**Website:** https://www.ebi.ac.uk/intact/

**Description:** EBI's molecular interaction database with standardized format.

**API Access:** PSICQUIC web services

**Key Features:**
- Standardized MI-TAB format
- Quality scores for interactions
- Cross-references to other databases

### 4. Reactome

**Website:** https://reactome.org

**Description:** Pathway database with interaction context.

**API Base URL:** `https://reactome.org/ContentService/`

**Key Features:**
- Pathway context for interactions
- Reaction-level detail
- Well-documented REST API

### 5. UniProt

**Website:** https://www.uniprot.org

**Description:** Primary protein database with interaction annotations.

**API Features:**
- Binary interactions from IntAct
- Subcomplexes information
- Can be queried via REST API

---

## Analysis Strategies

### Strategy 1: Network Topology Analysis

**Goal:** Understand the position of high pLLPS proteins in the interaction network

**Metrics to Calculate:**
1. **Degree Centrality:** Number of interaction partners
2. **Betweenness Centrality:** How often a protein lies on shortest paths between other proteins
3. **Clustering Coefficient:** How connected a protein's neighbors are to each other
4. **Hub Score:** Identification of highly connected proteins

**Implementation:**
```python
import networkx as nx

# Build network from interactions
G = nx.Graph()
for _, row in interactions_df.iterrows():
    G.add_edge(row['protein_a'], row['protein_b'], weight=row['score'])

# Calculate centrality metrics
degree_centrality = nx.degree_centrality(G)
betweenness_centrality = nx.betweenness_centrality(G)
clustering_coef = nx.clustering(G)
```

### Strategy 2: High pLLPS Subnetwork Analysis

**Goal:** Extract and analyze the subnetwork containing only high pLLPS proteins

**Steps:**
1. Define "high pLLPS" threshold (e.g., top 10%, > 0.7, > 0.9)
2. Extract subnetwork of high pLLPS proteins
3. Calculate network density and connectivity
4. Compare to random subnetworks of same size

**Key Questions:**
- Is the high pLLPS subnetwork more connected than expected by chance?
- Are there distinct clusters within high pLLPS proteins?

### Strategy 3: Enrichment Analysis

**Goal:** Determine if high pLLPS proteins are enriched for specific:
- GO terms (Biological Process, Molecular Function, Cellular Component)
- KEGG pathways
- Disease associations
- Protein domains

**Implementation:** Use STRING's enrichment endpoint or tools like:
- g:Profiler (https://biit.cs.ut.ee/gprofiler/)
- DAVID (https://david.ncifcrf.gov/)
- Enrichr (https://maayanlab.cloud/Enrichr/)

### Strategy 4: Community Detection

**Goal:** Identify functional modules/clusters in the interaction network

**Algorithms:**
1. **Louvain Algorithm:** Fast community detection
2. **Label Propagation:** Finds communities based on labels
3. **Spectral Clustering:** Uses eigenvalues of adjacency matrix

**Implementation:**
```python
import community.community_louvain as community_louvain

# Detect communities
partition = community_louvain.best_partition(G)

# Analyze pLLPS distribution in each community
for community_id in set(partition.values()):
    members = [node for node, comm in partition.items() if comm == community_id]
    avg_pllps = df[df['Entry'].isin(members)]['p(LLPS)'].mean()
    print(f"Community {community_id}: {len(members)} proteins, avg pLLPS: {avg_pllps:.3f}")
```

### Strategy 5: Comparison Analysis

**Goal:** Compare high pLLPS vs low pLLPS proteins

**Comparisons:**
1. Number of interaction partners (degree)
2. Types of interactions (physical vs predicted)
3. Interaction confidence scores
4. Functional annotations of partners

---

## Implementation Approach

### Phase 1: Data Preparation

```python
import pandas as pd

# Load LLPS data
df = pd.read_excel('Human Phase separation data.xlsx')

# Define high pLLPS threshold
high_pllps_threshold = 0.7  # or use percentile-based
high_pllps_proteins = df[df['p(LLPS)'] >= high_pllps_threshold]['Entry'].tolist()

print(f"Total proteins: {len(df)}")
print(f"High pLLPS proteins (>= {high_pllps_threshold}): {len(high_pllps_proteins)}")
```

### Phase 2: Fetch Interactions from STRING

```python
import requests
import time

def get_string_interactions(proteins, species=9606, score_threshold=400):
    """
    Fetch protein-protein interactions from STRING database.
    
    Parameters:
    -----------
    proteins : list
        List of protein identifiers (UniProt IDs)
    species : int
        NCBI taxonomy ID (9606 for human)
    score_threshold : int
        Minimum combined score (0-1000)
    
    Returns:
    --------
    pd.DataFrame : Interaction data
    """
    string_api_url = "https://string-db.org/api"
    output_format = "tsv"
    method = "network"
    
    # Process in batches of 200
    batch_size = 200
    all_interactions = []
    
    for i in range(0, len(proteins), batch_size):
        batch = proteins[i:i+batch_size]
        
        params = {
            "identifiers": "%0d".join(batch),  # Newline-separated
            "species": species,
            "caller_identity": "llps_analysis",
            "network_type": "physical",  # or "functional"
            "required_score": score_threshold
        }
        
        request_url = f"{string_api_url}/{output_format}/{method}"
        response = requests.post(request_url, data=params)
        
        if response.status_code == 200:
            # Parse TSV response
            from io import StringIO
            df_batch = pd.read_csv(StringIO(response.text), sep='\t')
            all_interactions.append(df_batch)
        
        time.sleep(1)  # Rate limiting
    
    return pd.concat(all_interactions, ignore_index=True)
```

### Phase 3: Network Analysis

```python
import networkx as nx
import matplotlib.pyplot as plt

def analyze_pllps_network(interactions_df, pllps_df, high_threshold=0.7):
    """
    Analyze the protein interaction network with pLLPS annotations.
    """
    # Build network
    G = nx.Graph()
    
    for _, row in interactions_df.iterrows():
        G.add_edge(
            row['preferredName_A'], 
            row['preferredName_B'],
            score=row['score']
        )
    
    # Add pLLPS as node attribute
    pllps_dict = dict(zip(pllps_df['Entry'], pllps_df['p(LLPS)']))
    nx.set_node_attributes(G, pllps_dict, 'pLLPS')
    
    # Calculate metrics
    results = {
        'total_nodes': G.number_of_nodes(),
        'total_edges': G.number_of_edges(),
        'density': nx.density(G),
        'avg_clustering': nx.average_clustering(G)
    }
    
    # Analyze high pLLPS subnetwork
    high_pllps_nodes = [n for n, attr in G.nodes(data=True) 
                        if attr.get('pLLPS', 0) >= high_threshold]
    
    high_subgraph = G.subgraph(high_pllps_nodes)
    results['high_pllps_nodes'] = len(high_pllps_nodes)
    results['high_pllps_edges'] = high_subgraph.number_of_edges()
    results['high_pllps_density'] = nx.density(high_subgraph) if len(high_pllps_nodes) > 1 else 0
    
    return G, results
```

### Phase 4: Visualization

```python
def visualize_pllps_network(G, high_threshold=0.7, output_file='pllps_network.html'):
    """
    Create interactive network visualization using Plotly.
    """
    import plotly.graph_objects as go
    
    # Get positions using spring layout
    pos = nx.spring_layout(G, k=2, iterations=50)
    
    # Create edge traces
    edge_x, edge_y = [], []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines'
    )
    
    # Create node traces
    node_x = [pos[node][0] for node in G.nodes()]
    node_y = [pos[node][1] for node in G.nodes()]
    node_pllps = [G.nodes[node].get('pLLPS', 0) for node in G.nodes()]
    node_text = [f"{node}<br>pLLPS: {G.nodes[node].get('pLLPS', 'N/A'):.3f}" 
                 for node in G.nodes()]
    
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        text=node_text,
        marker=dict(
            showscale=True,
            colorscale='RdYlBu_r',
            color=node_pllps,
            size=10,
            colorbar=dict(
                thickness=15,
                title='p(LLPS)',
                xanchor='left'
            )
        )
    )
    
    fig = go.Figure(data=[edge_trace, node_trace],
                   layout=go.Layout(
                       title='Protein Interaction Network Colored by pLLPS',
                       showlegend=False,
                       hovermode='closest',
                       xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
                   ))
    
    fig.write_html(output_file)
    return fig
```

---

## Sample Code

### Complete Analysis Pipeline

```python
"""
LLPS Protein Interaction Analysis Pipeline

This script analyzes protein-protein interactions for high pLLPS proteins
using the STRING database and network analysis techniques.
"""

import pandas as pd
import numpy as np
import requests
import networkx as nx
from collections import Counter
import time

# Configuration
PLLPS_THRESHOLD = 0.7  # Define high pLLPS
STRING_SCORE_THRESHOLD = 700  # High confidence interactions
SPECIES_ID = 9606  # Human

def load_llps_data(filepath):
    """Load LLPS protein data from Excel file."""
    df = pd.read_excel(filepath)
    print(f"Loaded {len(df)} proteins")
    print(f"pLLPS range: {df['p(LLPS)'].min():.3f} - {df['p(LLPS)'].max():.3f}")
    return df

def classify_proteins(df, threshold):
    """Classify proteins as high or low pLLPS."""
    df = df.copy()
    df['pLLPS_class'] = np.where(df['p(LLPS)'] >= threshold, 'high', 'low')
    
    high_count = (df['pLLPS_class'] == 'high').sum()
    low_count = (df['pLLPS_class'] == 'low').sum()
    
    print(f"High pLLPS (>= {threshold}): {high_count} proteins")
    print(f"Low pLLPS (< {threshold}): {low_count} proteins")
    
    return df

def fetch_interactions_batch(protein_ids, species=9606, score_threshold=400, batch_size=100):
    """
    Fetch interactions from STRING in batches.
    
    Note: In sandbox environments, this may not work due to network restrictions.
    For production use, ensure network access to string-db.org.
    """
    string_api_url = "https://string-db.org/api/json/network"
    
    all_interactions = []
    
    for i in range(0, len(protein_ids), batch_size):
        batch = protein_ids[i:i+batch_size]
        
        params = {
            "identifiers": "\r".join(batch),
            "species": species,
            "required_score": score_threshold,
            "caller_identity": "llps_analysis"
        }
        
        try:
            response = requests.post(string_api_url, data=params, timeout=30)
            if response.status_code == 200:
                interactions = response.json()
                all_interactions.extend(interactions)
                print(f"Batch {i//batch_size + 1}: Found {len(interactions)} interactions")
        except Exception as e:
            print(f"Error fetching batch {i//batch_size + 1}: {e}")
        
        time.sleep(1)  # Rate limiting
    
    return pd.DataFrame(all_interactions)

def build_network(interactions_df, pllps_df):
    """Build NetworkX graph from interactions."""
    G = nx.Graph()
    
    # Add edges
    for _, row in interactions_df.iterrows():
        G.add_edge(
            row['preferredName_A'],
            row['preferredName_B'],
            score=row.get('score', 0)
        )
    
    # Add pLLPS attribute to nodes
    pllps_dict = dict(zip(pllps_df['Entry'], pllps_df['p(LLPS)']))
    for node in G.nodes():
        G.nodes[node]['pLLPS'] = pllps_dict.get(node, None)
    
    return G

def analyze_hub_status(G, pllps_threshold=0.7, top_n=50):
    """Analyze whether high pLLPS proteins are network hubs."""
    # Calculate degree for all nodes
    degrees = dict(G.degree())
    
    # Get top hubs
    top_hubs = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:top_n]
    
    # Check pLLPS status of hubs
    high_pllps_hubs = 0
    hub_data = []
    
    for protein, degree in top_hubs:
        pllps = G.nodes[protein].get('pLLPS')
        is_high = pllps is not None and pllps >= pllps_threshold
        if is_high:
            high_pllps_hubs += 1
        hub_data.append({
            'protein': protein,
            'degree': degree,
            'pLLPS': pllps,
            'is_high_pLLPS': is_high
        })
    
    print(f"\nHub Analysis (top {top_n} by degree):")
    print(f"High pLLPS hubs: {high_pllps_hubs}/{top_n} ({100*high_pllps_hubs/top_n:.1f}%)")
    
    return pd.DataFrame(hub_data)

def analyze_interaction_enrichment(G, pllps_threshold=0.7):
    """
    Analyze if high pLLPS proteins preferentially interact with each other.
    """
    # Classify nodes
    high_pllps_nodes = set()
    low_pllps_nodes = set()
    
    for node in G.nodes():
        pllps = G.nodes[node].get('pLLPS')
        if pllps is not None:
            if pllps >= pllps_threshold:
                high_pllps_nodes.add(node)
            else:
                low_pllps_nodes.add(node)
    
    # Count interaction types
    high_high = 0
    high_low = 0
    low_low = 0
    
    for edge in G.edges():
        a, b = edge
        a_high = a in high_pllps_nodes
        b_high = b in high_pllps_nodes
        
        if a_high and b_high:
            high_high += 1
        elif a_high or b_high:
            high_low += 1
        else:
            low_low += 1
    
    total_edges = G.number_of_edges()
    
    print(f"\nInteraction Type Analysis:")
    print(f"High-High pLLPS: {high_high} ({100*high_high/total_edges:.1f}%)")
    print(f"High-Low pLLPS: {high_low} ({100*high_low/total_edges:.1f}%)")
    print(f"Low-Low pLLPS: {low_low} ({100*low_low/total_edges:.1f}%)")
    
    # Expected by chance (null model)
    n_high = len(high_pllps_nodes)
    n_low = len(low_pllps_nodes)
    n_total = n_high + n_low
    
    if n_total > 0:
        p_high = n_high / n_total
        expected_high_high = p_high * p_high
        expected_high_low = 2 * p_high * (1 - p_high)
        expected_low_low = (1 - p_high) ** 2
        
        print(f"\nExpected by chance:")
        print(f"High-High: {expected_high_high*100:.1f}%")
        print(f"High-Low: {expected_high_low*100:.1f}%")
        print(f"Low-Low: {expected_low_low*100:.1f}%")
        
        # Calculate enrichment
        observed_high_high = high_high / total_edges if total_edges > 0 else 0
        enrichment = observed_high_high / expected_high_high if expected_high_high > 0 else 0
        print(f"\nEnrichment of High-High interactions: {enrichment:.2f}x")
    
    return {
        'high_high': high_high,
        'high_low': high_low,
        'low_low': low_low,
        'total': total_edges
    }

def detect_communities(G, min_size=5):
    """Detect communities and analyze pLLPS distribution."""
    try:
        import community.community_louvain as community_louvain
        partition = community_louvain.best_partition(G)
    except ImportError:
        # Fallback to label propagation
        communities = list(nx.community.label_propagation_communities(G))
        partition = {}
        for i, community in enumerate(communities):
            for node in community:
                partition[node] = i
    
    # Analyze communities
    community_data = []
    
    for comm_id in set(partition.values()):
        members = [n for n, c in partition.items() if c == comm_id]
        if len(members) >= min_size:
            pllps_values = [G.nodes[n].get('pLLPS') for n in members if G.nodes[n].get('pLLPS') is not None]
            
            if pllps_values:
                community_data.append({
                    'community_id': comm_id,
                    'size': len(members),
                    'avg_pLLPS': np.mean(pllps_values),
                    'std_pLLPS': np.std(pllps_values),
                    'min_pLLPS': min(pllps_values),
                    'max_pLLPS': max(pllps_values),
                    'high_pLLPS_fraction': sum(1 for p in pllps_values if p >= 0.7) / len(pllps_values)
                })
    
    df_communities = pd.DataFrame(community_data)
    df_communities = df_communities.sort_values('avg_pLLPS', ascending=False)
    
    print(f"\nCommunity Analysis ({len(df_communities)} communities with size >= {min_size}):")
    print(df_communities.head(10).to_string())
    
    return df_communities, partition

# Main execution (example)
if __name__ == "__main__":
    # Load data
    df = load_llps_data('Human Phase separation data.xlsx')
    df = classify_proteins(df, PLLPS_THRESHOLD)
    
    # Get high pLLPS proteins
    high_pllps_proteins = df[df['pLLPS_class'] == 'high']['Entry'].tolist()
    print(f"\nFetching interactions for {len(high_pllps_proteins)} high pLLPS proteins...")
    
    # Note: This requires network access to STRING
    # interactions_df = fetch_interactions_batch(high_pllps_proteins)
    
    # For offline analysis, you can export protein IDs and use STRING web interface
    print("\nFor offline analysis:")
    print("1. Copy protein IDs to STRING web interface (https://string-db.org/)")
    print("2. Download the network in TSV format")
    print("3. Load the file and continue analysis")
```

---

## Expected Outcomes

### Scenario A: High pLLPS Proteins Form Clusters
**Findings:**
- High-High interactions are enriched above expectation
- Communities with high average pLLPS exist
- Subnetwork density is high

**Interpretation:** Phase-separating proteins may co-assemble into condensates

### Scenario B: High pLLPS Proteins are Network Hubs
**Findings:**
- High pLLPS proteins have above-average degree centrality
- High pLLPS proteins have high betweenness centrality
- Many low pLLPS proteins interact with high pLLPS proteins

**Interpretation:** Phase-separating proteins may act as scaffolds or organizers

### Scenario C: High pLLPS Proteins are Isolated
**Findings:**
- Low interaction density among high pLLPS proteins
- No enrichment of High-High interactions
- Random distribution across communities

**Interpretation:** Phase separation capability may not determine interaction partners

### Scenario D: Mixed Pattern
**Findings:**
- Some functional clusters of high pLLPS proteins
- Variable hub status depending on function
- Context-dependent interaction patterns

**Interpretation:** Multiple mechanisms at play; consider subcellular location and function

---

## Next Steps

1. **Data Collection:**
   - Export high pLLPS protein list
   - Fetch interactions using STRING API or web interface
   
2. **Network Construction:**
   - Build interaction network
   - Annotate nodes with pLLPS values
   
3. **Analysis:**
   - Calculate network metrics
   - Detect communities
   - Perform enrichment analysis
   
4. **Visualization:**
   - Create interactive network plots
   - Add to Streamlit dashboard
   
5. **Validation:**
   - Compare with known phase separation complexes
   - Validate against experimental data

---

## References

1. STRING Database: https://string-db.org
2. STRING MCP Server: https://github.com/Augmented-Nature/STRING-db-MCP-Server
3. BioGRID: https://thebiogrid.org
4. NetworkX Documentation: https://networkx.org
5. Phase Separation Database (PhaSepDB): http://db.phasep.pro/
6. DrLLPS Database: http://llps.biocuckoo.cn/

---

## Appendix: Using STRING MCP for LLPS Analysis

### Quick Start with STRING MCP

The STRING MCP server provides the most streamlined way to analyze protein interactions when working with AI assistants. Here's how to use it for LLPS analysis:

#### Step 1: Get High pLLPS Protein List
First, extract the high pLLPS proteins from your dataset:

```python
import pandas as pd

df = pd.read_excel('Human Phase separation data.xlsx')
high_pllps = df[df['p(LLPS)'] >= 0.9]['Entry'].head(50).tolist()
print(','.join(high_pllps))
```

#### Step 2: Query Interaction Network via MCP

With the STRING MCP server configured, you can ask the AI assistant:

```
Please use get_interaction_network to build an interaction network for these 
high pLLPS proteins: Q9Y6V0, P53350, Q9Y566, Q9Y520, Q9Y4H2
Set the species to 9606 (human) and minimum score to 700 (high confidence).
```

#### Step 3: Perform Functional Enrichment

```
Use get_functional_enrichment to analyze GO term and KEGG pathway enrichment
for the following high pLLPS proteins: [list of proteins]
```

#### Step 4: Identify Hub Proteins

```
For each of these high pLLPS proteins, use get_protein_interactions to find 
their interaction partners. Then count which proteins have the most partners 
to identify potential hub proteins.
```

### Sample MCP Queries for LLPS Analysis

| Question | MCP Query |
|----------|-----------|
| Do high pLLPS proteins interact with each other? | `get_interaction_network` with high pLLPS protein list |
| What functions are enriched in high pLLPS proteins? | `get_functional_enrichment` with protein list |
| Which high pLLPS proteins are hubs? | Multiple `get_protein_interactions` calls |
| Are there known phase separation annotations? | `get_protein_annotations` for each protein |
| Do mouse orthologs also have high pLLPS partners? | `find_homologs` for each protein |

### Advantages of MCP Approach for This Analysis

1. **Natural Language Queries**: Ask questions about protein interactions directly
2. **Automatic Batching**: MCP handles large protein lists automatically
3. **Integrated Analysis**: Combine interaction data with enrichment in one session
4. **Error Handling**: MCP provides clear error messages for invalid proteins
5. **Reproducibility**: MCP queries can be documented and replicated
