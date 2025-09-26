#!/usr/bin/env python3
"""
UMAP Parameter Optimization Script

This script loads the merged spatial data saved by merge_hippo.py and runs UMAP embedding
with different combinations of min_dist, spread, and gamma parameters to find optimal
parameters for data visualization and analysis.
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product
from pathlib import Path
import time

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Configure scanpy settings
sc.settings.verbosity = 3  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')

def load_merged_data():
    """Load the merged spatial data from the h5ad file."""
    data_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_standard.h5ad'
    if not os.path.exists(data_path):
        print(f"Error: Data file not found at {data_path}")
        print("Please run merge_hippo_standard.py first to generate the merged data.")
        sys.exit(1)
    
    print(f"Loading merged data from {data_path}...")
    adata = sc.read_h5ad(data_path)
    print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} genes")
    return adata

def setup_output_directory():
    """Create output directory structure for UMAP optimization results."""
    output_dir = Path('MERGE_SPATIAL_DATA/umap_optimization')
    output_dir.mkdir(exist_ok=True)
    
    # Create subdirectories for different parameter combinations
    plots_dir = output_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)
    
    # Create subdirectory for summary results
    summary_dir = output_dir / 'summary'
    summary_dir.mkdir(exist_ok=True)
    
    return output_dir, plots_dir, summary_dir

def run_umap_with_params(adata, min_dist, spread, gamma, random_state=42):
    """Run UMAP with specified parameters."""
    print(f"Running UMAP with min_dist={min_dist}, spread={spread}, gamma={gamma}")
    
    # Create a copy to avoid modifying the original data
    adata_copy = adata.copy()
    
    # Run UMAP with specified parameters
    sc.tl.umap(adata_copy,
               min_dist=min_dist,
               spread=spread,
               gamma=gamma,
               n_components=2,
               random_state=random_state)
    
    return adata_copy

def plot_umap_results(adata, min_dist, spread, gamma, plots_dir):
    """Generate and save UMAP plots for different metadata categories."""
    param_str = f"min_dist_{min_dist}_spread_{spread}_gamma_{gamma}"
    param_dir = plots_dir / param_str
    param_dir.mkdir(exist_ok=True)
    
    # Set up the figure parameters
    sc.settings.figdir = str(param_dir)
    
    plot_configs = [
        {'color': 'time_point', 'title': 'Time Points'},
        {'color': 'region', 'title': 'Regions'},
        {'color': 'type', 'title': 'Sample Types'},
        {'color': 'sample_id', 'title': 'Sample IDs'}
    ]
    
    # Check if leiden clustering exists, if not run it
    if 'leiden' not in adata.obs.columns:
        print("Running Leiden clustering...")
        sc.tl.leiden(adata, resolution=0.5)
    
    plot_configs.append({'color': 'leiden', 'title': 'Leiden Clusters'})
    
    plot_files = []
    for config in plot_configs:
        filename = f"umap_{config['color']}_{param_str}.png"
        plot_files.append(filename)
        
        # Create the plot
        sc.pl.umap(adata, 
                   color=config['color'], 
                   title=f"UMAP - {config['title']} (min_dist={min_dist}, spread={spread}, gamma={gamma})",
                   save=f"_{config['color']}_{param_str}.png",
                   show=False)
    
    return plot_files, param_dir

def calculate_umap_metrics(adata):
    """Calculate quality metrics for UMAP embedding."""
    # Basic metrics
    umap_coords = adata.obsm['X_umap']
    
    # Calculate spread of points
    x_range = np.ptp(umap_coords[:, 0])  # peak-to-peak range
    y_range = np.ptp(umap_coords[:, 1])
    total_spread = x_range + y_range
    
    # Calculate cluster separation (if leiden clustering exists)
    cluster_separation = None
    if 'leiden' in adata.obs.columns:
        # Calculate average distance between cluster centroids
        clusters = adata.obs['leiden'].unique()
        centroids = []
        for cluster in clusters:
            cluster_mask = adata.obs['leiden'] == cluster
            centroid = np.mean(umap_coords[cluster_mask], axis=0)
            centroids.append(centroid)
        
        if len(centroids) > 1:
            centroids = np.array(centroids)
            # Calculate pairwise distances between centroids
            from scipy.spatial.distance import pdist
            distances = pdist(centroids)
            cluster_separation = np.mean(distances)
    
    # Calculate local structure preservation (approximate)
    # This is a simplified metric - could be improved with more sophisticated measures
    n_neighbors = min(15, adata.n_obs - 1)
    from sklearn.neighbors import NearestNeighbors
    
    # Original space (PCA coordinates)
    if 'X_pca' in adata.obsm:
        original_space = adata.obsm['X_pca'][:, :n_neighbors]
    else:
        original_space = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    
    # UMAP space
    umap_space = umap_coords
    
    # Find neighbors in both spaces
    nbrs_orig = NearestNeighbors(n_neighbors=n_neighbors).fit(original_space)
    nbrs_umap = NearestNeighbors(n_neighbors=n_neighbors).fit(umap_space)
    
    # Sample a subset of points for efficiency
    n_sample = min(1000, adata.n_obs)
    sample_idx = np.random.choice(adata.n_obs, n_sample, replace=False)
    
    neighborhood_preservation = 0
    for idx in sample_idx:
        orig_neighbors = set(nbrs_orig.kneighbors([original_space[idx]], return_distance=False)[0])
        umap_neighbors = set(nbrs_umap.kneighbors([umap_space[idx]], return_distance=False)[0])
        overlap = len(orig_neighbors.intersection(umap_neighbors))
        neighborhood_preservation += overlap / n_neighbors
    
    neighborhood_preservation /= n_sample
    
    return {
        'x_range': x_range,
        'y_range': y_range, 
        'total_spread': total_spread,
        'cluster_separation': cluster_separation,
        'neighborhood_preservation': neighborhood_preservation
    }

def main():
    """Main function to run UMAP parameter optimization."""
    print("Starting UMAP parameter optimization...")
    
    # Load the merged data
    adata = load_merged_data()
    
    # Setup output directories
    output_dir, plots_dir, summary_dir = setup_output_directory()
    
    # Define parameter combinations to test
    min_dist_values = [0.01, 0.1, 0.3, 0.5, 0.8]
    spread_values = [0.5, 1.0, 1.5, 2.0]
    gamma_values = [1.0, 2.0, 5.0, 10.0]
    
    print(f"Testing {len(min_dist_values)} × {len(spread_values)} × {len(gamma_values)} = {len(min_dist_values) * len(spread_values) * len(gamma_values)} parameter combinations")
    
    # Store results for summary
    results = []
    
    # Test all parameter combinations
    total_combinations = len(min_dist_values) * len(spread_values) * len(gamma_values)
    current_combination = 0
    
    for min_dist, spread, gamma in product(min_dist_values, spread_values, gamma_values):
        current_combination += 1
        print(f"\n=== Combination {current_combination}/{total_combinations} ===")
        
        start_time = time.time()
        
        try:
            # Run UMAP with current parameters
            adata_umap = run_umap_with_params(adata, min_dist, spread, gamma)
            
            # Generate plots
            plot_files, param_dir = plot_umap_results(adata_umap, min_dist, spread, gamma, plots_dir)
            
            # Calculate metrics
            metrics = calculate_umap_metrics(adata_umap)
            
            # Store results
            result = {
                'min_dist': min_dist,
                'spread': spread,
                'gamma': gamma,
                'param_string': f"min_dist_{min_dist}_spread_{spread}_gamma_{gamma}",
                'runtime': time.time() - start_time,
                **metrics
            }
            results.append(result)
            
            print(f"Completed in {result['runtime']:.2f} seconds")
            print(f"Metrics - Total spread: {metrics['total_spread']:.2f}, "
                  f"Neighborhood preservation: {metrics['neighborhood_preservation']:.3f}")
            
        except Exception as e:
            print(f"Error with parameters min_dist={min_dist}, spread={spread}, gamma={gamma}: {str(e)}")
            continue
    
    # Create summary DataFrame and save results
    results_df = pd.DataFrame(results)
    
    # Save summary CSV
    summary_csv = summary_dir / 'umap_parameter_optimization_results.csv'
    results_df.to_csv(summary_csv, index=False)
    print(f"\nResults saved to {summary_csv}")
    
    # Create summary plots
    create_summary_plots(results_df, summary_dir)
    
    # Print top recommendations
    print_recommendations(results_df)
    
    print(f"\nUMAP parameter optimization completed!")
    print(f"Results saved in: {output_dir}")
    print(f"Individual plots saved in: {plots_dir}")
    print(f"Summary files saved in: {summary_dir}")

def create_summary_plots(results_df, summary_dir):
    """Create summary plots comparing different parameter combinations."""
    print("\nCreating summary plots...")
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # 1. Heatmap of total spread
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # For each gamma value, create a heatmap
    gamma_values = sorted(results_df['gamma'].unique())
    
    for i, gamma in enumerate(gamma_values[:4]):  # Show first 4 gamma values
        ax = axes[i//2, i%2]
        
        # Pivot data for heatmap
        subset = results_df[results_df['gamma'] == gamma]
        heatmap_data = subset.pivot(index='min_dist', columns='spread', values='total_spread')
        
        sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap='viridis', ax=ax)
        ax.set_title(f'Total Spread (gamma={gamma})')
        ax.set_xlabel('Spread')
        ax.set_ylabel('Min Distance')
    
    plt.tight_layout()
    plt.savefig(summary_dir / 'total_spread_heatmaps.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Neighborhood preservation heatmaps
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    for i, gamma in enumerate(gamma_values[:4]):
        ax = axes[i//2, i%2]
        
        subset = results_df[results_df['gamma'] == gamma]
        heatmap_data = subset.pivot(index='min_dist', columns='spread', values='neighborhood_preservation')
        
        sns.heatmap(heatmap_data, annot=True, fmt='.3f', cmap='plasma', ax=ax)
        ax.set_title(f'Neighborhood Preservation (gamma={gamma})')
        ax.set_xlabel('Spread')
        ax.set_ylabel('Min Distance')
    
    plt.tight_layout()
    plt.savefig(summary_dir / 'neighborhood_preservation_heatmaps.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Scatter plot: Total spread vs Neighborhood preservation
    plt.figure(figsize=(10, 8))
    
    for gamma in gamma_values:
        subset = results_df[results_df['gamma'] == gamma]
        plt.scatter(subset['total_spread'], subset['neighborhood_preservation'], 
                   label=f'gamma={gamma}', alpha=0.7, s=60)
    
    plt.xlabel('Total Spread')
    plt.ylabel('Neighborhood Preservation')
    plt.title('UMAP Parameter Trade-offs: Spread vs Neighborhood Preservation')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(summary_dir / 'spread_vs_preservation_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Runtime analysis
    plt.figure(figsize=(10, 6))
    results_df_sorted = results_df.sort_values('runtime')
    plt.plot(range(len(results_df_sorted)), results_df_sorted['runtime'], 'o-', alpha=0.7)
    plt.xlabel('Parameter Combination (sorted by runtime)')
    plt.ylabel('Runtime (seconds)')
    plt.title('UMAP Runtime for Different Parameter Combinations')
    plt.grid(True, alpha=0.3)
    plt.savefig(summary_dir / 'runtime_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Summary plots created successfully!")

def print_recommendations(results_df):
    """Print top parameter recommendations based on different criteria."""
    print("\n=== TOP RECOMMENDATIONS ===")
    
    # Best neighborhood preservation
    best_preservation = results_df.loc[results_df['neighborhood_preservation'].idxmax()]
    print(f"\nBest Neighborhood Preservation ({best_preservation['neighborhood_preservation']:.3f}):")
    print(f"  min_dist={best_preservation['min_dist']}, spread={best_preservation['spread']}, gamma={best_preservation['gamma']}")
    
    # Best total spread
    best_spread = results_df.loc[results_df['total_spread'].idxmax()]
    print(f"\nBest Total Spread ({best_spread['total_spread']:.2f}):")
    print(f"  min_dist={best_spread['min_dist']}, spread={best_spread['spread']}, gamma={best_spread['gamma']}")
    
    # Balanced approach (high preservation + reasonable spread)
    # Normalize metrics and create a balanced score
    results_df_norm = results_df.copy()
    results_df_norm['preservation_norm'] = (results_df['neighborhood_preservation'] - results_df['neighborhood_preservation'].min()) / (results_df['neighborhood_preservation'].max() - results_df['neighborhood_preservation'].min())
    results_df_norm['spread_norm'] = (results_df['total_spread'] - results_df['total_spread'].min()) / (results_df['total_spread'].max() - results_df['total_spread'].min())
    results_df_norm['balanced_score'] = 0.6 * results_df_norm['preservation_norm'] + 0.4 * results_df_norm['spread_norm']
    
    best_balanced = results_df_norm.loc[results_df_norm['balanced_score'].idxmax()]
    print(f"\nBest Balanced Score ({best_balanced['balanced_score']:.3f}):")
    print(f"  min_dist={best_balanced['min_dist']}, spread={best_balanced['spread']}, gamma={best_balanced['gamma']}")
    print(f"  Neighborhood preservation: {best_balanced['neighborhood_preservation']:.3f}")
    print(f"  Total spread: {best_balanced['total_spread']:.2f}")

if __name__ == "__main__":
    main()