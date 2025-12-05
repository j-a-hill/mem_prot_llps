import pandas as pd
import numpy as np
import requests
import time
import matplotlib.pyplot as plt
import seaborn as sns
import re

def parse_location(location_str):
    """
    Parse a subcellular location string from UniProt and return a list of location terms.
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

def add_location_columns(df):
    """Add parsed location columns to the dataframe."""
    if 'Subcellular location [CC]' not in df.columns:
        return df
    
    df = df.copy()
    # Add column with list of all parsed locations
    df['Location Categories'] = df['Subcellular location [CC]'].apply(parse_location)
    return df

def analyze_interactions_by_location(matched_df, full_dataset_df, locations, high_threshold=0.7, low_threshold=0.4):
    """
    Analyze interaction preferences for specific subcellular locations.
    """
    results = {}
    
    # Ensure location columns exist
    if 'Location Categories' not in full_dataset_df.columns:
        full_dataset_df = add_location_columns(full_dataset_df)
        
    # We also need location info for the matched_df proteins
    # We can merge it from full_dataset_df
    if 'Location Categories' not in matched_df.columns:
        # Create a map from Entry -> Location Categories
        loc_map = dict(zip(full_dataset_df['Entry'], full_dataset_df['Location Categories']))
        
        # Helper to get locations for a Uniprot ID
        def get_locs(uid):
            return loc_map.get(uid, [])
            
        matched_df['locs_a'] = matched_df['uniprot_a'].apply(get_locs)
        matched_df['locs_b'] = matched_df['uniprot_b'].apply(get_locs)
    
    for loc in locations:
        print(f"\nAnalyzing location: {loc}")
        
        # Filter full dataset to this location (for background probs)
        loc_full_df = full_dataset_df[full_dataset_df['Location Categories'].apply(lambda x: loc in x)]
        
        if len(loc_full_df) < 10:
            print(f"  Skipping {loc}: Too few proteins ({len(loc_full_df)})")
            continue
            
        # Filter interactions where BOTH proteins are in this location
        # This represents the "Nucleus Interaction Network".
        
        loc_matched_df = matched_df[
            matched_df['locs_a'].apply(lambda x: loc in x) & 
            matched_df['locs_b'].apply(lambda x: loc in x)
        ].copy()
        
        print(f"  Proteins in location: {len(loc_full_df)}")
        print(f"  Interactions within location: {len(loc_matched_df)}")
        
        # Count entries per class in this location
        class_counts = loc_full_df['pLLPS_class'].value_counts()
        print(f"  Class breakdown: High={class_counts.get('High', 0)}, Medium={class_counts.get('Medium', 0)}, Low={class_counts.get('Low', 0)}")
        
        if len(loc_matched_df) < 10:
            print("  Not enough interactions for analysis.")
            continue
            
        # Run matrix analysis
        loc_results = analyze_interaction_matrix(loc_matched_df, loc_full_df, high_threshold, low_threshold)
        if loc_results:
            results[loc] = loc_results
            
    return results

def plot_location_heatmaps(results):
    """
    Plot a grid of heatmaps for location-specific results.
    """
    if not results:
        print("No results to plot.")
        return
        
    n_locs = len(results)
    cols = 3
    rows = (n_locs + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows))
    if n_locs > 1:
        axes = axes.flatten()
    else:
        axes = [axes]
    
    for i, (loc, res) in enumerate(results.items()):
        ax = axes[i]
        sns.heatmap(res['enrichment'], annot=True, fmt=".2f", cmap="RdBu_r", center=1,
                    linewidths=.5, ax=ax, vmin=0, vmax=3, cbar=False)
        ax.set_title(f"{loc}\n(n={res['total_interactions']})")
        ax.set_xlabel('')
        ax.set_ylabel('')
        
    # Hide empty subplots
    for j in range(i+1, len(axes)):
        axes[j].axis('off')
        
    plt.tight_layout()
    plt.savefig('pllps_location_analysis.png', dpi=150)
    plt.show()
    print("Figure saved to: pllps_location_analysis.png")

def load_and_classify_data(filepath, high_threshold=0.7, low_threshold=0.4):
    """
    Load the dataset and classify proteins into High, Medium, and Low pLLPS classes.
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

def get_string_mapping(protein_ids, species=9606, batch_size=2000):
    """
    Map Uniprot IDs to STRING preferred names (Gene Names) and STRING IDs.
    Returns two dictionaries:
    1. string_to_uniprot: Maps STRING ID/Name -> Uniprot ID
    2. uniprot_to_string: Maps Uniprot ID -> STRING ID
    """
    url = "https://string-db.org/api/json/get_string_ids"
    string_to_uniprot = {}
    uniprot_to_string = {}
    
    total_batches = (len(protein_ids) + batch_size - 1) // batch_size
    print(f"Mapping {len(protein_ids)} proteins to STRING names in {total_batches} batches...")
    
    for i in range(0, len(protein_ids), batch_size):
        batch = protein_ids[i:i+batch_size]
        batch_num = i // batch_size + 1
        
        params = {
            "identifiers": "\r".join(batch),
            "species": species,
            "echo_query": 1
        }
        
        try:
            response = requests.post(url, data=params, timeout=60)
            if response.status_code == 200:
                data = response.json()
                for item in data:
                    # Map preferredName -> queryItem (Uniprot ID)
                    # Also map stringId -> queryItem
                    u_id = item['queryItem']
                    pref_name = item['preferredName']
                    string_id = item['stringId']
                    
                    string_to_uniprot[pref_name] = u_id
                    string_to_uniprot[string_id] = u_id
                    
                    # Map Uniprot -> STRING ID (preferred for querying)
                    uniprot_to_string[u_id] = string_id
                    
                print(f"  Batch {batch_num}/{total_batches}: Mapped {len(data)} items")
            else:
                print(f"  Batch {batch_num}/{total_batches}: Error {response.status_code}")
        except Exception as e:
            print(f"  Batch {batch_num}/{total_batches}: {e}")
            
        time.sleep(1)
        
    return string_to_uniprot, uniprot_to_string

def fetch_interaction_partners(protein_ids, uniprot_map=None, species=9606, score_threshold=700, batch_size=100):
    """
    Fetch interaction partners for a list of proteins from STRING.
    
    Parameters:
    - protein_ids: List of UniProt IDs
    - uniprot_map: Dictionary mapping Uniprot IDs to STRING IDs (optional but recommended)
    - species: NCBI taxonomy ID (9606 = human)
    - score_threshold: Minimum confidence score (0-1000)
    - batch_size: Proteins per API request
    
    Returns:
    - DataFrame with interactions
    """
    # Use interaction_partners API to get ALL partners, not just interactions within the set
    string_api_url = "https://string-db.org/api/json/interaction_partners"
    all_interactions = []
    
    # Convert Uniprot IDs to STRING IDs if map provided
    # This improves matching success rate significantly
    query_ids = []
    if uniprot_map:
        for pid in protein_ids:
            if pid in uniprot_map:
                query_ids.append(uniprot_map[pid])
            else:
                # If no map, try using the ID directly (might fail if STRING doesn't recognize it)
                query_ids.append(pid)
    else:
        query_ids = protein_ids
        
    total_batches = (len(query_ids) + batch_size - 1) // batch_size
    print(f"Fetching interactions for {len(query_ids)} proteins in {total_batches} batches...")
    
    for i in range(0, len(query_ids), batch_size):
        batch = query_ids[i:i+batch_size]
        batch_num = i // batch_size + 1
        
        params = {
            "identifiers": "\r".join(batch),
            "species": species,
            "required_score": score_threshold,
            "caller_identity": "pllps_analysis",
            "limit": 50  # Limit partners per protein to avoid explosion
        }
        
        try:
            response = requests.post(string_api_url, data=params, timeout=60)
            if response.status_code == 200:
                interactions = response.json()
                all_interactions.extend(interactions)
                print(f"  Batch {batch_num}/{total_batches}: {len(interactions)} interactions")
            else:
                print(f"  Batch {batch_num}/{total_batches}: Error {response.status_code}")
        except Exception as e:
            print(f"  Batch {batch_num}/{total_batches}: {e}")
        
        time.sleep(1)  # Rate limiting
    
    if all_interactions:
        return pd.DataFrame(all_interactions)
    return pd.DataFrame()

def match_interactors_to_pllps(interactions_df, pllps_df, string_map=None):
    """
    Match interaction partners to the pLLPS dataset.
    
    Returns DataFrame with both proteins' pLLPS values.
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
             if protein_a_name in pllps_lookup: uniprot_a = protein_a_name
             elif protein_a_name in entry_name_lookup: uniprot_a = protein_a_name # This is actually pLLPS value lookup, wait.
        
        # Resolve Protein B to Uniprot ID
        uniprot_b = None
        if string_map:
            uniprot_b = string_map.get(protein_b_name) or string_map.get(protein_b_id)
            
        if not uniprot_b:
             if protein_b_name in pllps_lookup: uniprot_b = protein_b_name

        # Get pLLPS scores
        # If uniprot_a is found in map, it is a Uniprot ID.
        # If not found in map, we tried to use name as ID.
        
        pllps_a = pllps_lookup.get(uniprot_a)
        if pllps_a is None and not string_map: # Fallback to old logic if map not provided
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

def analyze_interaction_matrix(matched_df, full_dataset_df, high_threshold=0.7, low_threshold=0.4):
    """
    Analyze interaction preferences between High, Medium, and Low pLLPS classes.
    Creates a 3x3 Enrichment Matrix.
    """
    # Filter to interactions where both proteins have pLLPS scores
    complete_df = matched_df.dropna(subset=['pllps_a', 'pllps_b']).copy()
    
    if len(complete_df) == 0:
        print("No interactions with complete pLLPS data found.")
        return None
    
    # Define classes for the interaction pairs
    def get_class(score):
        if score >= high_threshold: return 'High'
        elif score >= low_threshold: return 'Medium'
        else: return 'Low'
        
    complete_df['class_a'] = complete_df['pllps_a'].apply(get_class)
    complete_df['class_b'] = complete_df['pllps_b'].apply(get_class)
    
    # Calculate Genomic Background Probabilities
    # P(High), P(Medium), P(Low) in the whole genome
    total_proteins = len(full_dataset_df)
    class_counts = full_dataset_df['pLLPS_class'].value_counts()
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
    # Note: Interactions are undirected. A-B is the same as B-A.
    # We will count each interaction twice (A->B and B->A) to make the matrix symmetric
    # and easier to interpret as "Given a protein of Class X, what is the prob of partner Class Y?"
    
    for _, row in complete_df.iterrows():
        c1, c2 = row['class_a'], row['class_b']
        observed_counts.loc[c1, c2] += 1
        observed_counts.loc[c2, c1] += 1
        
    # Calculate Expected Counts
    # Expected(A, B) = Total_Interactions * P(A) * P(B)
    # Since we double-counted above, Total_Interactions is 2 * len(complete_df)
    total_observations = observed_counts.sum().sum()
    
    for c1 in classes:
        for c2 in classes:
            # The expected probability of picking a pair (c1, c2) by chance
            # is P(c1) * P(c2)
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

def plot_interaction_heatmap(results, output_file='pllps_interaction_matrix.png'):
    """
    Visualize the interaction enrichment matrix as a heatmap.
    """
    if results is not None:
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot Heatmap of Enrichment
        sns.heatmap(results['enrichment'], annot=True, fmt=".2f", cmap="RdBu_r", center=1,
                    linewidths=.5, ax=ax, vmin=0, vmax=3)
        
        ax.set_title('Interaction Enrichment Matrix\n(>1 = Enriched, <1 = Depleted)', fontsize=14)
        ax.set_xlabel('Partner Protein Class', fontsize=12)
        ax.set_ylabel('Query Protein Class', fontsize=12)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=150)
        plt.show()
        
        print(f"\nFigure saved to: {output_file}")
        
        # Print interpretation
        print("\nInterpretation:")
        print("- Values > 1.0 (Red): These classes interact MORE than expected by chance.")
        print("- Values < 1.0 (Blue): These classes interact LESS than expected by chance.")
        print("- Diagonal (High-High, Med-Med, Low-Low): Shows homotypic preference.")
