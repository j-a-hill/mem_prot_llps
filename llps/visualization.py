"""
Section 6 - Visualization Functions.

Functions for plotting interaction heatmaps and printing analysis reports.
"""

import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Any, Optional


def plot_interaction_heatmap(results: Optional[Dict], output_file: str = 'pllps_interaction_matrix.png') -> None:
    """
    Visualize the interaction enrichment matrix as a heatmap.
    
    Parameters
    ----------
    results : dict
        Results dictionary from analyze_interaction_matrix
    output_file : str
        Path to save the output figure
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


def plot_location_heatmaps(results: Dict) -> None:
    """
    Plot a grid of heatmaps for location-specific results.
    
    Parameters
    ----------
    results : dict
        Dictionary with location-specific enrichment results
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


def print_analysis_report(results: Dict[str, Any]) -> None:
    """
    Print formatted analysis report.
    
    Parameters
    ----------
    results : dict
        Dictionary with network analysis results
    """
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
