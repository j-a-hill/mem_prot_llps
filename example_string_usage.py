"""
Example: Using LLPS Functions Module

This script demonstrates how to use the llps_functions module
to fetch and analyze protein interactions.
"""

import pandas as pd
from llps_functions import (
    fetch_string_interactions,
    match_interactions_to_pllps,
    analyze_interaction_enrichment,
    save_interactions_to_cache
)

# Example 1: Fetch interactions for a few proteins
print("=" * 60)
print("Example 1: Fetch STRING Interactions")
print("=" * 60)

# List of protein IDs (examples: p53, BRCA1, BRCA2)
protein_ids = ['P04637', 'P38398', 'P51587']

print(f"\nFetching interactions for {len(protein_ids)} proteins...")
interactions_df, errors = fetch_string_interactions(
    protein_ids,
    species=9606,  # Human
    score_threshold=700,  # High confidence
    progress_callback=lambda msg: print(f"  {msg}")
)

print(f"\nRetrieved {len(interactions_df)} interactions")
if errors:
    print("\nErrors:")
    for err in errors:
        print(f"  - {err}")

if len(interactions_df) > 0:
    print("\nFirst few interactions:")
    print(interactions_df[['preferredName_A', 'preferredName_B', 'score']].head())


# Example 2: Match interactions to pLLPS dataset
print("\n" + "=" * 60)
print("Example 2: Match to pLLPS Dataset")
print("=" * 60)

# Load your pLLPS data
pllps_data_file = "Human Phase separation data.xlsx"
try:
    pllps_df = pd.read_excel(pllps_data_file, engine='openpyxl')
    print(f"\nLoaded pLLPS data: {len(pllps_df)} proteins")
    
    # Get high pLLPS proteins
    high_pllps = pllps_df[pllps_df['p(LLPS)'] >= 0.7]
    print(f"High pLLPS proteins (>=0.7): {len(high_pllps)}")
    
    # Use top 50 for faster testing
    test_ids = high_pllps['Entry'].head(50).tolist()
    print(f"\nFetching interactions for top {len(test_ids)} high pLLPS proteins...")
    
    interactions_df, errors = fetch_string_interactions(
        test_ids,
        score_threshold=700,
        progress_callback=lambda msg: print(f"  {msg}")
    )
    
    if len(interactions_df) > 0:
        print(f"\nFound {len(interactions_df)} interactions")
        
        # Match to pLLPS dataset
        matched_df = match_interactions_to_pllps(interactions_df, pllps_df)
        print(f"Matched {len(matched_df)} interactions to pLLPS data")
        
        # Example 3: Analyze enrichment
        print("\n" + "=" * 60)
        print("Example 3: Enrichment Analysis")
        print("=" * 60)
        
        results = analyze_interaction_enrichment(matched_df, threshold=0.7)
        
        if results:
            print(f"\nTotal interactions analyzed: {results['total']}")
            print(f"  High-High: {results['high_high']} ({results['high_high']/results['total']*100:.1f}%)")
            print(f"  High-Low:  {results['high_low']} ({results['high_low']/results['total']*100:.1f}%)")
            print(f"  Low-Low:   {results['low_low']} ({results['low_low']/results['total']*100:.1f}%)")
            print(f"\nExpected under random model:")
            print(f"  High-High: {results['expected_hh']:.1f}%")
            print(f"  High-Low:  {results['expected_hl']:.1f}%")
            print(f"  Low-Low:   {results['expected_ll']:.1f}%")
            print(f"\nEnrichment factor: {results['enrichment']:.2f}x")
            if results['p_value'] is not None:
                print(f"Chi-squared test: χ² = {results['chi2']:.2f}, p = {results['p_value']:.2e}")
                if results['p_value'] < 0.05:
                    if results['enrichment'] > 1:
                        print("✅ Significant: High pLLPS proteins preferentially interact!")
                    else:
                        print("✅ Significant: High pLLPS proteins avoid each other.")
                else:
                    print("⚠️  Not significant (p >= 0.05)")
        
        # Example 4: Save to cache
        print("\n" + "=" * 60)
        print("Example 4: Save to Cache")
        print("=" * 60)
        
        cache_file = save_interactions_to_cache(
            interactions_df,
            score_threshold=700,
            output_dir="data"
        )
        print(f"\nCache saved to: {cache_file}")
        print("This cache can be used for offline analysis")
        
except FileNotFoundError:
    print(f"\nNote: {pllps_data_file} not found.")
    print("Place your LLPS data file in the same directory to run full examples.")


print("\n" + "=" * 60)
print("Examples Complete!")
print("=" * 60)
print("\nTo use these functions in your own code:")
print("  from llps_functions import fetch_string_interactions")
print("\nTo use in Jupyter notebook:")
print("  %matplotlib inline")
print("  import llps_functions as lf")
print("\nTo run the Shiny app:")
print("  shiny run shiny_app.py --port 8000")
print("=" * 60)
