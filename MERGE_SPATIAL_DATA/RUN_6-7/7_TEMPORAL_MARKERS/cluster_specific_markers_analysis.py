#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cluster-specific Temporal Markers Analysis

This script performs differential expression analysis between ctrl and mut samples
within specific leiden_new clusters at each time point. Only clusters with at least
20 cells in both ctrl and mut conditions are analyzed.

Author: Generated for SRF Linda project
Date: 2024
"""

# %%
import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# %%
# Set up project paths
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# %%
# Configure scanpy
sc.settings.verbosity = 3  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')

# %%
# Load merged spatial data
print("Loading merged spatial data...")
merged_data_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
if not os.path.exists(merged_data_path):
    print(f"Error: Merged data file not found at {merged_data_path}")
    print("Please run the merge script first.")
    sys.exit(1)

adata = sc.read_h5ad(merged_data_path)
print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} genes")

# %%
# Set up output directory
output_dir = 'MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/cluster_results'
os.makedirs(output_dir, exist_ok=True)
print(f"Results will be saved to: {output_dir}")

# %%
# Set scanpy figure directory
sc.settings.figdir = output_dir

# %%
# Check available time points and clusters
time_points = adata.obs['time_point'].unique()
print(f"Available time points: {time_points}")

if 'leiden_new' not in adata.obs.columns:
    print("Error: 'leiden_new' column not found in adata.obs")
    print("Available columns:", adata.obs.columns.tolist())
    sys.exit(1)

clusters = adata.obs['leiden_new'].unique()
print(f"Available leiden_new clusters: {sorted(clusters)}")

# %%
# Function to check cluster validity (at least 20 cells in both ctrl and mut)
def check_cluster_validity(adata_subset, cluster, min_cells=50):
    """
    Check if a cluster has at least min_cells in both ctrl and mut conditions
    """
    cluster_data = adata_subset[adata_subset.obs['leiden_new'] == cluster]
    
    ctrl_count = len(cluster_data[cluster_data.obs['type'] == 'ctrl'])
    mut_count = len(cluster_data[cluster_data.obs['type'] == 'mut'])
    
    return ctrl_count >= min_cells and mut_count >= min_cells, ctrl_count, mut_count

# %%
# Function to perform cluster-specific differential expression analysis
def perform_cluster_de_analysis(adata_subset, time_point, cluster, output_dir):
    """
    Perform differential expression analysis between ctrl and mut for a specific cluster at a time point
    """
    print(f"\n=== Analyzing cluster {cluster} at time point {time_point} ===")
    
    # Filter for specific cluster
    cluster_data = adata_subset[adata_subset.obs['leiden_new'] == cluster].copy()
    
    # Check if we have both ctrl and mut
    types = cluster_data.obs['type'].unique()
    if len(types) < 2:
        print(f"Skipping cluster {cluster}: insufficient sample types")
        return None
    
    # Perform differential expression using Wilcoxon rank-sum test
    sc.tl.rank_genes_groups(
        cluster_data, 
        groupby='type', 
        groups=['mut'], 
        reference='ctrl',
        method='wilcoxon',
        key_added=f'de_{time_point}_cluster_{cluster}'
    )
    
    # Extract results
    de_results = sc.get.rank_genes_groups_df(
        cluster_data, 
        group='mut', 
        key=f'de_{time_point}_cluster_{cluster}'
    )
    
    # Add metadata
    de_results['time_point'] = time_point
    de_results['cluster'] = cluster
    
    # Calculate additional statistics
    de_results['abs_logfoldchanges'] = np.abs(de_results['logfoldchanges'])
    de_results['-log10_pvals_adj'] = -np.log10(de_results['pvals_adj'] + 1e-300)
    
    # Define significance thresholds
    pval_threshold = 0.05
    logfc_threshold = 0.5
    
    de_results['significant'] = (
        (de_results['pvals_adj'] < pval_threshold) & 
        (de_results['abs_logfoldchanges'] > logfc_threshold)
    )
    
    # Classify genes as upregulated or downregulated
    de_results['regulation'] = 'not_significant'
    de_results.loc[
        (de_results['significant']) & (de_results['logfoldchanges'] > 0), 
        'regulation'
    ] = 'upregulated'
    de_results.loc[
        (de_results['significant']) & (de_results['logfoldchanges'] < 0), 
        'regulation'
    ] = 'downregulated'
    
    ctrl_count = len(cluster_data[cluster_data.obs['type'] == 'ctrl'])
    mut_count = len(cluster_data[cluster_data.obs['type'] == 'mut'])
    
    print(f"Cluster {cluster} - Ctrl cells: {ctrl_count}, Mut cells: {mut_count}")
    print(f"Total genes analyzed: {len(de_results)}")
    print(f"Significantly upregulated: {sum(de_results['regulation'] == 'upregulated')}")
    print(f"Significantly downregulated: {sum(de_results['regulation'] == 'downregulated')}")
    
    # Save results
    results_file = os.path.join(output_dir, f'de_results_{time_point}_cluster_{cluster}_mut_vs_ctrl.csv')
    de_results.to_csv(results_file, index=False)
    print(f"Results saved to: {results_file}")
    
    return de_results

# %%
# Function to create cluster-specific volcano plot
def create_cluster_volcano_plot(de_results, time_point, cluster, output_dir):
    """
    Create volcano plot for cluster-specific differential expression results
    """
    plt.figure(figsize=(10, 8))
    
    # Create volcano plot
    colors = []
    for _, row in de_results.iterrows():
        if row['regulation'] == 'upregulated':
            colors.append('red')
        elif row['regulation'] == 'downregulated':
            colors.append('blue')
        else:
            colors.append('gray')
    
    plt.scatter(de_results['logfoldchanges'], de_results['-log10_pvals_adj'], 
                c=colors, alpha=0.6, s=20)
    
    # Add significance thresholds
    plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.7)
    plt.axvline(x=0.5, color='gray', linestyle='--', alpha=0.7)
    plt.axvline(x=-0.5, color='gray', linestyle='--', alpha=0.7)
    
    # Labels and title
    plt.xlabel('Log2 Fold Change (mut vs ctrl)')
    plt.ylabel('-Log10 Adjusted P-value')
    plt.title(f'Volcano Plot - {time_point.upper()} Cluster {cluster} (mut vs ctrl)')
    
    # Add legend
    upregulated_count = sum(de_results['regulation'] == 'upregulated')
    downregulated_count = sum(de_results['regulation'] == 'downregulated')
    not_significant_count = sum(de_results['regulation'] == 'not_significant')
    
    legend_elements = [
        plt.scatter([], [], c='gray', alpha=0.6, s=50, label=f'Not Significant ({not_significant_count})'),
        plt.scatter([], [], c='blue', alpha=0.6, s=50, label=f'Downregulated ({downregulated_count})'),
        plt.scatter([], [], c='red', alpha=0.6, s=50, label=f'Upregulated ({upregulated_count})')
    ]
    
    plt.legend(handles=legend_elements, loc='upper right')
    
    # Add threshold lines info
    plt.text(0.02, 0.98, f'p-adj = 0.05\n|log2FC| = 0.5', 
             transform=plt.gca().transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    # Save plot
    volcano_file = os.path.join(output_dir, f'volcano_plot_{time_point}_cluster_{cluster}.png')
    plt.savefig(volcano_file, dpi=300, bbox_inches='tight')
    print(f"Volcano plot saved to: {volcano_file}")
    plt.close()

# %%
# Function to create cluster-specific heatmap
def create_cluster_heatmap(adata_subset, de_results, time_point, cluster, output_dir, n_genes=20):
    """
    Create heatmap of top differentially expressed genes for a specific cluster
    """
    # Filter for specific cluster
    cluster_data = adata_subset[adata_subset.obs['leiden_new'] == cluster].copy()
    
    # Get top upregulated and downregulated genes
    significant_genes = de_results[de_results['significant']].copy()
    
    if len(significant_genes) == 0:
        print(f"No significant genes found for cluster {cluster} at {time_point}")
        return
    
    # Sort by absolute log fold change and take top genes
    significant_genes = significant_genes.sort_values('abs_logfoldchanges', ascending=False)
    top_genes = significant_genes.head(min(n_genes, len(significant_genes)))['names'].tolist()
    
    if len(top_genes) == 0:
        print(f"No genes to plot for cluster {cluster} at {time_point}")
        return
    
    # Create heatmap
    sc.pl.heatmap(
        cluster_data, 
        var_names=top_genes, 
        groupby='type',
        dendrogram=False,
        save=f'_heatmap_{time_point}_cluster_{cluster}.png',
        show=False,
        figsize=(8, max(6, len(top_genes) * 0.3))
    )
    
    print(f"Heatmap saved for cluster {cluster} at {time_point} with {len(top_genes)} genes")

# %%
# Main analysis loop
print("\n=== Starting Cluster-specific Analysis ===")

# Find valid comparisons (time points with both ctrl and mut)
valid_comparisons = []
all_statistics = []  # Store detailed statistics for each comparison

for tp in time_points:
    tp_data = adata[adata.obs['time_point'] == tp]
    types_available = tp_data.obs['type'].unique()
    if 'ctrl' in types_available and 'mut' in types_available:
        valid_comparisons.append(tp)
        ctrl_cells = len(tp_data[tp_data.obs['type'] == 'ctrl'])
        mut_cells = len(tp_data[tp_data.obs['type'] == 'mut'])
        print(f"Time point {tp}: {ctrl_cells} ctrl cells, {mut_cells} mut cells")
        
        # Store time point statistics
        all_statistics.append({
            'time_point': tp,
            'cluster': 'ALL',
            'comparison_type': 'timepoint_overview',
            'ctrl_cells': ctrl_cells,
            'mut_cells': mut_cells,
            'total_cells': ctrl_cells + mut_cells,
            'analysis_status': 'valid',
            'reason': 'Both ctrl and mut samples available'
        })
    else:
        print(f"Skipping time point {tp}: missing ctrl or mut samples")
        # Store skipped time point statistics
        all_statistics.append({
            'time_point': tp,
            'cluster': 'ALL',
            'comparison_type': 'timepoint_overview',
            'ctrl_cells': len(tp_data[tp_data.obs['type'] == 'ctrl']) if 'ctrl' in types_available else 0,
            'mut_cells': len(tp_data[tp_data.obs['type'] == 'mut']) if 'mut' in types_available else 0,
            'total_cells': len(tp_data),
            'analysis_status': 'skipped',
            'reason': f'Missing sample types: {list(set(["ctrl", "mut"]) - set(types_available))}'
        })

print(f"\nWill analyze time points: {valid_comparisons}")

# %%
# Analyze each time point and cluster combination
all_de_results = []
valid_cluster_info = []

for time_point in valid_comparisons:
    print(f"\n{'='*60}")
    print(f"Processing time point: {time_point}")
    print(f"{'='*60}")
    
    # Subset data for this time point
    adata_subset = adata[adata.obs['time_point'] == time_point].copy()
    
    # Check each cluster for validity
    time_point_clusters = adata_subset.obs['leiden_new'].unique()
    
    for cluster in sorted(time_point_clusters):
        is_valid, ctrl_count, mut_count = check_cluster_validity(adata_subset, cluster)
        
        # Store cluster statistics regardless of validity
        cluster_stats = {
            'time_point': time_point,
            'cluster': cluster,
            'comparison_type': 'cluster_specific',
            'ctrl_cells': ctrl_count,
            'mut_cells': mut_count,
            'total_cells': ctrl_count + mut_count,
            'analysis_status': 'analyzed' if is_valid else 'skipped',
            'reason': 'Sufficient cells (≥20 in both conditions)' if is_valid else f'Insufficient cells (ctrl: {ctrl_count}, mut: {mut_count})'
        }
        all_statistics.append(cluster_stats)
        
        if is_valid:
            print(f"\nCluster {cluster} is valid - Ctrl: {ctrl_count}, Mut: {mut_count}")
            
            # Perform differential expression analysis
            de_results = perform_cluster_de_analysis(adata_subset, time_point, cluster, output_dir)
            
            if de_results is not None:
                all_de_results.append(de_results)
                
                # Create visualizations
                create_cluster_volcano_plot(de_results, time_point, cluster, output_dir)
                create_cluster_heatmap(adata_subset, de_results, time_point, cluster, output_dir)
                
                # Store cluster info
                valid_cluster_info.append({
                    'time_point': time_point,
                    'cluster': cluster,
                    'ctrl_cells': ctrl_count,
                    'mut_cells': mut_count,
                    'total_genes': len(de_results),
                    'significant_genes': sum(de_results['significant']),
                    'upregulated': sum(de_results['regulation'] == 'upregulated'),
                    'downregulated': sum(de_results['regulation'] == 'downregulated')
                })
                
                # Update cluster statistics with DE results
                cluster_stats.update({
                    'total_genes_analyzed': len(de_results),
                    'significant_genes': sum(de_results['significant']),
                    'upregulated_genes': sum(de_results['regulation'] == 'upregulated'),
                    'downregulated_genes': sum(de_results['regulation'] == 'downregulated'),
                    'pct_significant': (sum(de_results['significant']) / len(de_results)) * 100 if len(de_results) > 0 else 0
                })
        else:
            print(f"Skipping cluster {cluster} - Ctrl: {ctrl_count}, Mut: {mut_count} (< 20 cells)")
            # Add empty DE statistics for skipped clusters
            cluster_stats.update({
                'total_genes_analyzed': 0,
                'significant_genes': 0,
                'upregulated_genes': 0,
                'downregulated_genes': 0,
                'pct_significant': 0
            })

# %%
# Create comprehensive statistics file
print("\n=== Creating Comprehensive Statistics ===")

# Create detailed statistics DataFrame
statistics_df = pd.DataFrame(all_statistics)

# Add additional metadata
statistics_df['analysis_date'] = pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
statistics_df['min_cells_threshold'] = 20
statistics_df['meets_threshold'] = statistics_df.apply(
    lambda row: (row['ctrl_cells'] >= 20 and row['mut_cells'] >= 20) if row['comparison_type'] == 'cluster_specific' else True, 
    axis=1
)

# Calculate additional statistics
statistics_df['ctrl_mut_ratio'] = statistics_df.apply(
    lambda row: row['ctrl_cells'] / row['mut_cells'] if row['mut_cells'] > 0 else float('inf'), 
    axis=1
)

# Save comprehensive statistics
statistics_file = os.path.join(output_dir, 'analysis_statistics_comprehensive.csv')
statistics_df.to_csv(statistics_file, index=False)
print(f"Comprehensive statistics saved to: {statistics_file}")

# Create summary statistics by time point
summary_by_timepoint = statistics_df[statistics_df['comparison_type'] == 'cluster_specific'].groupby('time_point').agg({
    'cluster': 'count',
    'ctrl_cells': 'sum',
    'mut_cells': 'sum',
    'total_cells': 'sum',
    'meets_threshold': 'sum',
    'total_genes_analyzed': 'sum',
    'significant_genes': 'sum',
    'upregulated_genes': 'sum',
    'downregulated_genes': 'sum'
}).rename(columns={
    'cluster': 'total_clusters',
    'meets_threshold': 'analyzed_clusters'
}).reset_index()

summary_by_timepoint['pct_clusters_analyzed'] = (summary_by_timepoint['analyzed_clusters'] / summary_by_timepoint['total_clusters']) * 100
summary_by_timepoint['avg_ctrl_cells_per_cluster'] = summary_by_timepoint['ctrl_cells'] / summary_by_timepoint['total_clusters']
summary_by_timepoint['avg_mut_cells_per_cluster'] = summary_by_timepoint['mut_cells'] / summary_by_timepoint['total_clusters']

# Save timepoint summary
timepoint_summary_file = os.path.join(output_dir, 'timepoint_summary_statistics.csv')
summary_by_timepoint.to_csv(timepoint_summary_file, index=False)
print(f"Timepoint summary saved to: {timepoint_summary_file}")

# %%
# Create summary results
if all_de_results:
    print("\n=== Creating DE Results Summary ===")
    
    # Combine all results
    combined_results = pd.concat(all_de_results, ignore_index=True)
    
    # Save combined results
    combined_file = os.path.join(output_dir, 'combined_cluster_de_results_all_timepoints.csv')
    combined_results.to_csv(combined_file, index=False)
    print(f"Combined results saved to: {combined_file}")
    
    # Create cluster summary
    cluster_summary_df = pd.DataFrame(valid_cluster_info)
    cluster_summary_file = os.path.join(output_dir, 'cluster_summary_statistics.csv')
    cluster_summary_df.to_csv(cluster_summary_file, index=False)
    print(f"Cluster summary saved to: {cluster_summary_file}")
    
    # Display summaries
    print("\n=== CLUSTER ANALYSIS SUMMARY ===")
    print(cluster_summary_df.to_string(index=False))
    
    print("\n=== TIMEPOINT SUMMARY ===")
    print(summary_by_timepoint.to_string(index=False))
    
    print("\n=== ANALYSIS STATISTICS OVERVIEW ===")
    total_comparisons = len(statistics_df[statistics_df['comparison_type'] == 'cluster_specific'])
    analyzed_comparisons = len(statistics_df[(statistics_df['comparison_type'] == 'cluster_specific') & (statistics_df['analysis_status'] == 'analyzed')])
    skipped_comparisons = total_comparisons - analyzed_comparisons
    
    print(f"Total cluster comparisons evaluated: {total_comparisons}")
    print(f"Cluster comparisons analyzed: {analyzed_comparisons}")
    print(f"Cluster comparisons skipped: {skipped_comparisons}")
    print(f"Analysis success rate: {(analyzed_comparisons/total_comparisons)*100:.1f}%")
    
    # Show distribution of cell counts
    cluster_stats = statistics_df[statistics_df['comparison_type'] == 'cluster_specific']
    print(f"\nCell count distribution (ctrl):")
    print(f"  Min: {cluster_stats['ctrl_cells'].min()}")
    print(f"  Max: {cluster_stats['ctrl_cells'].max()}")
    print(f"  Mean: {cluster_stats['ctrl_cells'].mean():.1f}")
    print(f"  Median: {cluster_stats['ctrl_cells'].median():.1f}")
    
    print(f"\nCell count distribution (mut):")
    print(f"  Min: {cluster_stats['mut_cells'].min()}")
    print(f"  Max: {cluster_stats['mut_cells'].max()}")
    print(f"  Mean: {cluster_stats['mut_cells'].mean():.1f}")
    print(f"  Median: {cluster_stats['mut_cells'].median():.1f}")
    
    # Create summary visualization
    if len(cluster_summary_df) > 0:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        
        # Plot 1: Number of significant genes per cluster
        cluster_summary_df['cluster_label'] = cluster_summary_df['time_point'] + '_C' + cluster_summary_df['cluster'].astype(str)
        ax1.bar(range(len(cluster_summary_df)), cluster_summary_df['significant_genes'], color='skyblue')
        ax1.set_xlabel('Cluster')
        ax1.set_ylabel('Significant Genes')
        ax1.set_title('Significant Genes per Cluster')
        ax1.set_xticks(range(len(cluster_summary_df)))
        ax1.set_xticklabels(cluster_summary_df['cluster_label'], rotation=45, fontsize=6)
        
        # Plot 2: Up vs Down regulated genes
        x = np.arange(len(cluster_summary_df))
        width = 0.35
        ax2.bar(x - width/2, cluster_summary_df['upregulated'], width, label='Upregulated', color='red', alpha=0.7)
        ax2.bar(x + width/2, cluster_summary_df['downregulated'], width, label='Downregulated', color='blue', alpha=0.7)
        ax2.set_xlabel('Cluster')
        ax2.set_ylabel('Number of Genes')
        ax2.set_title('Up vs Down Regulated Genes per Cluster')
        ax2.set_xticks(x)
        ax2.set_xticklabels(cluster_summary_df['cluster_label'], rotation=45, fontsize=6)
        ax2.legend()
        
        # Plot 3: Cell counts per cluster
        ax3.bar(x - width/2, cluster_summary_df['ctrl_cells'], width, label='Ctrl', color='green', alpha=0.7)
        ax3.bar(x + width/2, cluster_summary_df['mut_cells'], width, label='Mut', color='orange', alpha=0.7)
        ax3.set_xlabel('Cluster')
        ax3.set_ylabel('Number of Cells')
        ax3.set_title('Cell Counts per Cluster')
        ax3.set_xticks(x)
        ax3.set_xticklabels(cluster_summary_df['cluster_label'], rotation=45, fontsize=6)
        ax3.legend()
        
        # Plot 4: Percentage of significant genes
        cluster_summary_df['pct_significant'] = (cluster_summary_df['significant_genes'] / cluster_summary_df['total_genes']) * 100
        ax4.bar(range(len(cluster_summary_df)), cluster_summary_df['pct_significant'], color='purple', alpha=0.7)
        ax4.set_xlabel('Cluster')
        ax4.set_ylabel('% Significant Genes')
        ax4.set_title('Percentage of Significant Genes per Cluster')
        ax4.set_xticks(range(len(cluster_summary_df)))
        ax4.set_xticklabels(cluster_summary_df['cluster_label'], rotation=45, fontsize=6)
        
        plt.tight_layout()
        
        # Save summary plot
        summary_plot_file = os.path.join(output_dir, 'cluster_analysis_summary.png')
        plt.savefig(summary_plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Summary plot saved to: {summary_plot_file}")

else:
    print("No valid cluster differential expression results generated.")

print("\n=== Cluster-specific Analysis Complete ===")
print(f"All results saved to: {output_dir}")