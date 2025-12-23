#!/usr/bin/env python3
"""
Generate STRING Interaction Cache

This script fetches protein interactions from STRING database for high pLLPS
proteins and saves them to a cache file. This is useful for deploying the
Shiny app in environments without network access to string-db.org.

Usage:
    python generate_string_cache.py [--threshold 0.7] [--score 700] [--max-proteins 500]

The cache file will be saved to: data/string_cache_{score}.json
"""

import argparse
import json
import pandas as pd
from pathlib import Path
from llps_functions import (
    load_llps_data,
    get_high_pllps_proteins,
    fetch_string_interactions,
    save_interactions_to_cache
)


def main():
    parser = argparse.ArgumentParser(description='Generate STRING interaction cache')
    parser.add_argument('--threshold', type=float, default=0.7,
                       help='High pLLPS threshold (default: 0.7)')
    parser.add_argument('--score', type=int, default=700,
                       help='STRING confidence score threshold (default: 700)')
    parser.add_argument('--max-proteins', type=int, default=500,
                       help='Maximum proteins to analyze (default: 500)')
    parser.add_argument('--output-dir', type=str, default='data',
                       help='Output directory (default: data)')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("STRING Interaction Cache Generator")
    print("=" * 60)
    
    # Load data
    print("\n1. Loading LLPS data...")
    df = load_llps_data()
    
    # Get high pLLPS proteins
    print(f"\n2. Identifying high pLLPS proteins (threshold={args.threshold})...")
    df_classified, high_pllps_ids = get_high_pllps_proteins(df, threshold=args.threshold)
    print(f"   Found {len(high_pllps_ids)} high pLLPS proteins")
    
    # Limit to max proteins
    if len(high_pllps_ids) > args.max_proteins:
        print(f"   Limiting to top {args.max_proteins} proteins")
        # Sort by pLLPS score and take top proteins
        high_pllps_df = df_classified[df_classified['Entry'].isin(high_pllps_ids)]
        high_pllps_df = high_pllps_df.sort_values('p(LLPS)', ascending=False)
        high_pllps_ids = high_pllps_df['Entry'].head(args.max_proteins).tolist()
    
    # Fetch interactions
    print(f"\n3. Fetching STRING interactions (score>={args.score})...")
    print("   This may take several minutes...")
    
    def progress_print(msg):
        print(f"   {msg}")
    
    interactions_df, errors = fetch_string_interactions(
        high_pllps_ids,
        score_threshold=args.score,
        batch_size=100,
        progress_callback=progress_print,
        use_cache=False  # Don't use cache when generating cache
    )
    
    if errors:
        print("\n   Errors encountered:")
        for err in errors:
            print(f"     - {err}")
    
    print(f"\n4. Retrieved {len(interactions_df)} interactions")
    
    # Save to cache
    cache_file = save_interactions_to_cache(
        interactions_df,
        score_threshold=args.score,
        output_dir=args.output_dir
    )
    
    print(f"\n5. Cache saved to: {cache_file}")
    file_size = Path(cache_file).stat().st_size / 1024
    print(f"   File size: {file_size:.1f} KB")
    
    print("\n" + "=" * 60)
    print("✅ Cache generation complete!")
    print("=" * 60)
    print("\nTo use this cache in the Shiny app:")
    print(f"1. Ensure {cache_file} is in the data/ directory")
    print("2. The app will automatically load cached data if network is unavailable")
    print()


if __name__ == "__main__":
    main()
