#!/usr/bin/env python3
"""
Gene Set Scoring and Heatmap Visualization

This script performs gene set scoring at single cell level and creates heatmap visualizations
showing gene set scores across different clusters, conditions, regions, and time points.
Additionally, it creates temporal heatmaps showing how gene scores change across time points
(e16, p0, p7, p30) for different conditions, regions, and clusters.
"""

import scanpy as sc
import squidpy as sq
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Define paths
merged_h5ad_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
gene_sets_path = 'DATA/GENE_LISTS/output_combined/import_ready.csv'
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/gene_set_heatmaps'

# Create output directory
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

def load_gene_sets(gene_sets_path):
    """
    Load gene sets from CSV file.
    
    Parameters:
    -----------
    gene_sets_path : str
        Path to the gene sets CSV file
    
    Returns:
    --------
    dict
        Dictionary with gene set names as keys and lists of genes as values
    """
    print(f"Loading gene sets from: {gene_sets_path}")
    
    # Read the CSV file
    gene_sets_df = pd.read_csv(gene_sets_path, header=None, 
                              names=['gene_set', 'gene', 'col3', 'col4'])
    
    # Group genes by gene set
    gene_sets_dict = {}
    for gene_set, group in gene_sets_df.groupby('gene_set'):
        gene_sets_dict[gene_set] = group['gene'].tolist()
    
    print(f"Loaded {len(gene_sets_dict)} gene sets:")
    for gene_set, genes in gene_sets_dict.items():
        print(f"  - {gene_set}: {len(genes)} genes")
    
    return gene_sets_dict

def calculate_gene_set_scores(adata, gene_sets_dict):
    """
    Calculate gene set scores using scanpy's score_genes function.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object containing the gene expression data
    gene_sets_dict : dict
        Dictionary with gene set names as keys and lists of genes as values
    
    Returns:
    --------
    AnnData
        Updated AnnData object with gene set scores in .obs
    """
    print("Calculating gene set scores...")
    
    # Make a copy to avoid modifying the original
    adata_copy = adata.copy()
    
    # Calculate scores for each gene set
    for gene_set_name, genes in gene_sets_dict.items():
        print(f"  Scoring gene set: {gene_set_name}")
        
        # Filter genes to only those present in the dataset
        genes_in_data = [gene for gene in genes if gene in adata_copy.var.index]
        genes_missing = [gene for gene in genes if gene not in adata_copy.var.index]
        
        if len(genes_missing) > 0:
            print(f"    Warning: {len(genes_missing)} genes not found in dataset: {genes_missing[:5]}{'...' if len(genes_missing) > 5 else ''}")
        
        if len(genes_in_data) == 0:
            print(f"    Error: No genes from {gene_set_name} found in dataset. Skipping.")
            continue
        
        print(f"    Using {len(genes_in_data)}/{len(genes)} genes from {gene_set_name}")
        
        # Calculate gene set score using scanpy
        # The score will be stored in adata.obs with the gene set name
        sc.tl.score_genes(adata_copy, gene_list=genes_in_data, 
                         score_name=f"{gene_set_name}_score", use_raw=False)
    
    return adata_copy

def create_cluster_heatmap(adata, gene_sets_dict, time_point, output_dir):
    """
    Create heatmap showing gene set scores for leiden_new clusters within a given time point,
    across different conditions (ctrl and mut) and regions.
    All subplots in the generated file use the same color scale.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names as keys and lists of genes as values
    time_point : str
        Time point to analyze (e.g., 'e16', 'p0', 'p7', 'p30')
    output_dir : str
        Directory to save the plots
    """
    print(f"Creating cluster heatmap for time point: {time_point}")
    
    # Filter data for the specified time point
    time_adata = adata[adata.obs['time_point'] == time_point].copy()
    
    if len(time_adata) == 0:
        print(f"  No data found for time point: {time_point}")
        return
    
    # Get unique conditions and regions
    conditions = time_adata.obs['type'].unique()
    regions = time_adata.obs['region'].unique()
    
    # First, collect all data to calculate min/max for consistent color scaling within this file
    all_scores_for_timepoint = []
    all_heatmap_data = {}
    
    for region in regions:
        # Prepare data for heatmap
        heatmap_data = []
        
        # Get samples for this region
        region_samples = [s for s in time_adata.obs['sample_id'].unique() if region in str(s)]
        
        for condition in conditions:
            # Get samples for this condition and region
            condition_samples = [s for s in region_samples if condition in str(s)]
            
            if not condition_samples:
                continue
                
            # Get data for this condition and region
            condition_adata = time_adata[time_adata.obs['sample_id'].isin(condition_samples)]
            
            # Group by leiden_new clusters and calculate mean gene set scores
            for cluster in sorted(condition_adata.obs['leiden_new'].unique()):
                cluster_data = condition_adata[condition_adata.obs['leiden_new'] == cluster]
                
                if len(cluster_data) == 0:
                    continue
                
                # Calculate mean scores for each gene set
                row_data = {
                    'Condition': condition,
                    'Region': region,
                    'Cluster': cluster
                }
                
                for gene_set_name in gene_sets_dict.keys():
                    score_column = f"{gene_set_name}_score"
                    if score_column in cluster_data.obs.columns:
                        score = cluster_data.obs[score_column].mean()
                        row_data[gene_set_name] = score
                        all_scores_for_timepoint.append(score)
                
                heatmap_data.append(row_data)
        
        all_heatmap_data[region] = heatmap_data
    
    # Calculate min and max for this time point
    if all_scores_for_timepoint:
        file_min = min(all_scores_for_timepoint)
        file_max = max(all_scores_for_timepoint)
        print(f"  Using color scale range for this file: [{file_min:.3f}, {file_max:.3f}]")
    else:
        file_min = -1
        file_max = 1
        print(f"  No scores found, using default range [-1, 1]")
    
    # Create a figure with subplots for each region
    fig, axes = plt.subplots(len(regions), 1, figsize=(12, 8 * len(regions)))
    if len(regions) == 1:
        axes = [axes]
    
    for region_idx, region in enumerate(regions):
        ax = axes[region_idx]
        
        heatmap_data = all_heatmap_data[region]
        
        if not heatmap_data:
            print(f"  No data available for {region} at time point {time_point}")
            continue
        
        # Convert to DataFrame
        heatmap_df = pd.DataFrame(heatmap_data)
        
        # Set up the data for heatmap
        plot_data = heatmap_df.pivot_table(
            index=['Condition', 'Region', 'Cluster'],
            values=[col for col in heatmap_df.columns if col not in ['Condition', 'Region', 'Cluster']],
            aggfunc='mean'
        )
        
        # Create heatmap with consistent color scale for this file
        sns.heatmap(plot_data, ax=ax, cmap='coolwarm', center=0,
                   vmin=file_min, vmax=file_max,
                   cbar_kws={'label': 'Gene Set Score'})
        
        ax.set_title(f'Gene Set Scores by Cluster - Region {region}\nTime Point: {time_point}')
        ax.set_xlabel('Gene Sets')
        ax.set_ylabel('Condition, Region, Cluster')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'cluster_heatmap_timepoint_{time_point}.png'),
               bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"  Saved: cluster_heatmap_timepoint_{time_point}.png")

def create_celltype_heatmap(adata, gene_sets_dict, cluster, output_dir):
    """
    Create heatmap showing gene set scores for a given leiden_new cluster across different samples,
    conditions (ctrl and mut), regions, and time points.
    The heatmap uses a color scale specific to this cluster's data range.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names as keys and lists of genes as values
    cluster : str
        Leiden cluster to analyze
    output_dir : str
        Directory to save the plots
    """
    print(f"Creating cluster heatmap for: {cluster}")
    
    # Filter data for the specified cluster
    cluster_adata = adata[adata.obs['leiden_new'] == cluster].copy()
    
    if len(cluster_adata) == 0:
        print(f"  No data found for cluster: {cluster}")
        return
    
    # Define the desired time point order
    time_point_order = ["e16", "p0", "p7", "p30"]
    
    # Get unique conditions and regions
    conditions = cluster_adata.obs['type'].unique()
    regions = cluster_adata.obs['region'].unique()
    
    # Prepare data for heatmap and collect all scores for color scaling
    heatmap_data = []
    all_scores_for_cluster = []
    
    # Use the defined time point order instead of alphabetical sorting
    for time_point in time_point_order:
        # Only process time points that exist in the data
        if time_point not in cluster_adata.obs['time_point'].unique():
            continue
            
        for condition in conditions:
            for region in regions:
                # Get samples for this combination
                samples = [s for s in cluster_adata.obs['sample_id'].unique()
                          if time_point in str(s) and condition in str(s) and region in str(s)]
                
                if not samples:
                    continue
                
                # Get data for this combination
                combo_adata = cluster_adata[cluster_adata.obs['sample_id'].isin(samples)]
                
                if len(combo_adata) == 0:
                    continue
                
                # Calculate mean scores for each gene set
                row_data = {
                    'Time_Point': time_point,
                    'Condition': condition,
                    'Region': region
                }
                
                for gene_set_name in gene_sets_dict.keys():
                    score_column = f"{gene_set_name}_score"
                    if score_column in combo_adata.obs.columns:
                        score = combo_adata.obs[score_column].mean()
                        row_data[gene_set_name] = score
                        all_scores_for_cluster.append(score)
                
                heatmap_data.append(row_data)
    
    # Calculate min and max for this cluster
    if all_scores_for_cluster:
        file_min = min(all_scores_for_cluster)
        file_max = max(all_scores_for_cluster)
        print(f"  Using color scale range for this cluster: [{file_min:.3f}, {file_max:.3f}]")
    else:
        file_min = -1
        file_max = 1
        print(f"  No scores found, using default range [-1, 1]")
    
    # Convert to DataFrame
    heatmap_df = pd.DataFrame(heatmap_data)
    
    if heatmap_df.empty:
        print(f"  No data available for cluster: {cluster}")
        return
    
    # Convert Time_Point to categorical with the correct order to maintain sorting
    heatmap_df['Time_Point'] = pd.Categorical(
        heatmap_df['Time_Point'],
        categories=time_point_order,
        ordered=True
    )
    
    # Set up the data for heatmap
    plot_data = heatmap_df.pivot_table(
        index=['Time_Point', 'Condition', 'Region'],
        values=[col for col in heatmap_df.columns if col not in ['Time_Point', 'Condition', 'Region']],
        aggfunc='mean'
    )
    
    # Create heatmap with color scale specific to this cluster
    plt.figure(figsize=(12, 8))
    sns.heatmap(plot_data, cmap='coolwarm', center=0,
               vmin=file_min, vmax=file_max,
               cbar_kws={'label': 'Gene Set Score'})
    
    plt.title(f'Gene Set Scores for Cluster {cluster}\nAcross Time Points, Conditions, and Regions')
    plt.xlabel('Gene Sets')
    plt.ylabel('Time Point, Condition, Region')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'cluster_heatmap_{cluster}.png'),
               bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"  Saved: cluster_heatmap_{cluster}.png")

def create_temporal_heatmap(adata, gene_sets_dict, output_dir):
    """
    Create heatmap showing gene set scores with time points, conditions, and regions on x-axis
    and clusters on y-axis. One plot per gene set, combining all regions.
    All heatmaps use the same color scale for consistent comparison.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names as keys and lists of genes as values
    output_dir : str
        Directory to save the plots
    """
    print("Creating temporal heatmaps showing gene set scores across time points, conditions, and regions...")

    # Define the desired time point order
    time_point_order = ["e16", "p0", "p7", "p30"]

    # Get unique regions and clusters
    regions = sorted(adata.obs['region'].unique())
    clusters = sorted(adata.obs['leiden_new'].unique())
    conditions = sorted(adata.obs['type'].unique())

    # First pass: calculate global min and max values across all gene sets for consistent color scaling
    all_scores = []
    for gene_set_name in gene_sets_dict.keys():
        score_column = f"{gene_set_name}_score"

        if score_column not in adata.obs.columns:
            continue

        # Collect all non-NaN scores for this gene set
        gene_set_scores = adata.obs[score_column].dropna()
        all_scores.extend(gene_set_scores.tolist())

    # Calculate global min and max for consistent color scaling
    if all_scores:
        global_min = min(all_scores)
        global_max = max(all_scores)
        print(f"  Global score range for consistent color scaling: [{global_min:.3f}, {global_max:.3f}]")
    else:
        global_min = -1
        global_max = 1
        print("  No scores found, using default range [-1, 1]")

    # Create a separate heatmap for each gene set (combining all regions)
    for gene_set_name in gene_sets_dict.keys():
        score_column = f"{gene_set_name}_score"

        if score_column not in adata.obs.columns:
            print(f"  Score column {score_column} not found, skipping")
            continue

        print(f"  Processing gene set: {gene_set_name}")

        # Prepare data for heatmap
        heatmap_data = []

        # For each cluster, calculate mean scores across time points, conditions, and regions
        for cluster in clusters:
            cluster_adata = adata[adata.obs['leiden_new'] == cluster]

            if len(cluster_adata) == 0:
                continue

            # For each time point, condition, and region combination
            for time_point in time_point_order:
                for condition in conditions:
                    for region in regions:
                        # Get data for this specific combination
                        combo_adata = cluster_adata[
                            (cluster_adata.obs['time_point'] == time_point) &
                            (cluster_adata.obs['type'] == condition) &
                            (cluster_adata.obs['region'] == region)
                        ]

                        if len(combo_adata) == 0:
                            # Add NaN for missing combinations
                            heatmap_data.append({
                                'Cluster': cluster,
                                'Time_Condition_Region': f"{time_point}_{condition}_{region}",
                                'Score': np.nan
                            })
                        else:
                            # Calculate mean score
                            mean_score = combo_adata.obs[score_column].mean()
                            heatmap_data.append({
                                'Cluster': cluster,
                                'Time_Condition_Region': f"{time_point}_{condition}_{region}",
                                'Score': mean_score
                            })

        # Convert to DataFrame
        heatmap_df = pd.DataFrame(heatmap_data)

        if heatmap_df.empty:
            print(f"  No data available for gene set {gene_set_name}")
            continue

        # Pivot to create matrix for heatmap
        plot_data = heatmap_df.pivot(
            index='Cluster',
            columns='Time_Condition_Region',
            values='Score'
        )

        # Reorder columns to match time_point_order, conditions, and regions
        ordered_columns = [f"{tp}_{cond}_{region}" for tp in time_point_order for cond in conditions for region in regions]
        # Keep only columns that exist in the data
        ordered_columns = [col for col in ordered_columns if col in plot_data.columns]
        plot_data = plot_data[ordered_columns]

        # Create heatmap with consistent color scale
        plt.figure(figsize=(24, 12))  # Increased width to accommodate more columns
        sns.heatmap(plot_data, cmap='coolwarm', center=0,
                   vmin=global_min, vmax=global_max,  # Use global min/max for consistent scaling
                   cbar_kws={'label': 'Gene Set Score'},
                   annot=False, fmt='.2f')

        plt.title(f'Gene Set Scores: {gene_set_name}\nAll Regions Combined')
        plt.xlabel('Time Point & Condition & Region')
        plt.ylabel('Cluster')

        # Rotate x-axis labels for better readability
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)

        plt.tight_layout()

        # Create filename
        safe_gene_set = gene_set_name.replace('/', '_').replace('\\', '_')
        filename = f'temporal_heatmap_all_regions_{safe_gene_set}.png'
        plt.savefig(os.path.join(output_dir, filename),
                   bbox_inches='tight', dpi=300)
        plt.close()

        print(f"  Saved: {filename}")

    print("Temporal heatmaps creation completed!")

def create_summary_report(adata, gene_sets_dict, output_dir):
    """
    Create a summary report of gene set scores.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names as keys and lists of genes as values
    output_dir : str
        Directory to save the report
    """
    print("Creating summary report...")
    
    # Create summary statistics
    summary_data = []
    
    for gene_set_name in gene_sets_dict.keys():
        score_column = f"{gene_set_name}_score"
        
        if score_column not in adata.obs.columns:
            continue
        
        scores = adata.obs[score_column]
        
        # Overall statistics
        overall_stats = {
            'Gene_Set': gene_set_name,
            'Sample': 'All',
            'N_Cells': len(scores),
            'Mean_Score': scores.mean(),
            'Std_Score': scores.std(),
            'Min_Score': scores.min(),
            'Max_Score': scores.max(),
            'N_Genes_in_Set': len(gene_sets_dict[gene_set_name]),
            'N_Genes_Found': len([g for g in gene_sets_dict[gene_set_name] if g in adata.var.index])
        }
        summary_data.append(overall_stats)
        
        # Per-sample statistics
        for sample_id in adata.obs['sample_id'].unique():
            sample_scores = adata.obs[adata.obs['sample_id'] == sample_id][score_column]
            
            sample_stats = {
                'Gene_Set': gene_set_name,
                'Sample': sample_id,
                'N_Cells': len(sample_scores),
                'Mean_Score': sample_scores.mean(),
                'Std_Score': sample_scores.std(),
                'Min_Score': sample_scores.min(),
                'Max_Score': sample_scores.max(),
                'N_Genes_in_Set': len(gene_sets_dict[gene_set_name]),
                'N_Genes_Found': len([g for g in gene_sets_dict[gene_set_name] if g in adata.var.index])
            }
            summary_data.append(sample_stats)
    
    # Convert to DataFrame and save
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(output_dir, 'gene_set_scores_summary.csv'), index=False)
    
    # Create a formatted text report
    with open(os.path.join(output_dir, 'gene_set_scores_report.txt'), 'w') as f:
        f.write("Gene Set Scoring Summary Report\n")
        f.write("=" * 50 + "\n\n")
        
        f.write(f"Total gene sets analyzed: {len(gene_sets_dict)}\n")
        f.write(f"Total cells: {len(adata)}\n")
        f.write(f"Total samples: {len(adata.obs['sample_id'].unique())}\n\n")
        
        # Gene set overview
        f.write("Gene Sets Overview:\n")
        f.write("-" * 20 + "\n")
        for gene_set_name, genes in gene_sets_dict.items():
            score_column = f"{gene_set_name}_score"
            if score_column in adata.obs.columns:
                genes_found = len([g for g in genes if g in adata.var.index])
                f.write(f"{gene_set_name}: {genes_found}/{len(genes)} genes found\n")
            else:
                f.write(f"{gene_set_name}: SKIPPED (no genes found)\n")
        
        f.write(f"\nDetailed statistics saved to: gene_set_scores_summary.csv\n")
        f.write(f"Heatmaps saved to: {output_dir}\n")

def main():
    """Main function to run the gene set scoring and heatmap visualization pipeline."""
    
    print("Gene Set Scoring and Heatmap Visualization Pipeline")
    print("=" * 60)
    
    # Load the merged AnnData object
    if os.path.exists(merged_h5ad_path):
        print(f"Loading merged data from: {merged_h5ad_path}")
        adata = sc.read_h5ad(merged_h5ad_path)
        print(f"Loaded data: {adata.n_obs} cells x {adata.n_vars} genes")
    else:
        print(f"Error: Merged data file not found at {merged_h5ad_path}")
        sys.exit(1)
    
    # Load gene sets
    if os.path.exists(gene_sets_path):
        gene_sets_dict = load_gene_sets(gene_sets_path)
    else:
        print(f"Error: Gene sets file not found at {gene_sets_path}")
        sys.exit(1)
    
    # Calculate gene set scores
    adata_scored = calculate_gene_set_scores(adata, gene_sets_dict)
    
    # Create cluster heatmaps for each time point
    time_point_order = ["e16", "p0", "p7", "p30"]
    time_points = [tp for tp in time_point_order if tp in adata_scored.obs['time_point'].unique()]
    for time_point in time_points:
        create_cluster_heatmap(adata_scored, gene_sets_dict, time_point, output_dir)
    
    # Create cell type heatmaps for each leiden_new cluster (as proxy for cell types)
    # Note: If actual cell type annotations are available, replace this section
    unique_clusters = sorted(adata_scored.obs['leiden_new'].unique())
    for cluster in unique_clusters:  # Create for all clusters
        create_celltype_heatmap(adata_scored, gene_sets_dict, cluster, output_dir)
    
    # Create temporal heatmap showing gene set scores across time points
    # create_temporal_heatmap(adata_scored, gene_sets_dict, output_dir)
    
    # Create summary report
    create_summary_report(adata_scored, gene_sets_dict, output_dir)
    
    print("\n" + "=" * 60)
    print("Pipeline completed successfully!")
    print(f"Results saved to: {output_dir}")
    print("=" * 60)

if __name__ == "__main__":
    main()