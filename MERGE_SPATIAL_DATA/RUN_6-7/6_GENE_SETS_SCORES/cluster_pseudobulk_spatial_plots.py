#!/usr/bin/env python3
"""
Cluster-Level Pseudo-Bulk Gene Set Scoring and Spatial Visualization

This script calculates gene set scores at the cluster level (pseudo-bulk approach)
and displays these cluster-specific scores on spatial plots. All cells from the
same cluster will show the same score value.

The workflow:
1. Load spatial transcriptomics data
2. Load gene sets
3. Calculate pseudo-bulk expression for each cluster
4. Score gene sets on pseudo-bulk cluster profiles
5. Assign cluster scores to all cells in that cluster
6. Visualize on spatial plots
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
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/cluster_pseudobulk_spatial_plots'

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

def create_pseudobulk_by_cluster(adata, cluster_key='leiden_new'):
    """
    Create pseudo-bulk expression profiles for each cluster.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object containing the gene expression data
    cluster_key : str
        Key in adata.obs containing cluster assignments

    Returns:
    --------
    AnnData
        AnnData object where each observation is a cluster pseudo-bulk profile
    """
    print(f"Creating pseudo-bulk profiles for clusters...")

    clusters = sorted(adata.obs[cluster_key].unique())
    print(f"  Found {len(clusters)} clusters: {clusters}")

    # Calculate mean expression for each cluster
    pseudobulk_data = []
    cluster_names = []

    for cluster in clusters:
        cluster_cells = adata[adata.obs[cluster_key] == cluster]
        # Calculate mean expression across all cells in the cluster
        mean_expr = cluster_cells.X.mean(axis=0)
        if hasattr(mean_expr, 'A1'):  # If sparse matrix
            mean_expr = mean_expr.A1
        pseudobulk_data.append(mean_expr)
        cluster_names.append(cluster)

    # Create new AnnData object with pseudo-bulk data
    pseudobulk_data = np.array(pseudobulk_data)
    pseudobulk_adata = sc.AnnData(
        X=pseudobulk_data,
        obs=pd.DataFrame(index=cluster_names, data={cluster_key: cluster_names}),
        var=adata.var.copy()
    )

    print(f"  Created pseudo-bulk AnnData: {pseudobulk_adata.n_obs} clusters x {pseudobulk_adata.n_vars} genes")

    return pseudobulk_adata

def calculate_cluster_gene_set_scores(pseudobulk_adata, gene_sets_dict):
    """
    Calculate gene set scores for each cluster using pseudo-bulk profiles.

    Parameters:
    -----------
    pseudobulk_adata : AnnData
        Pseudo-bulk AnnData object with one observation per cluster
    gene_sets_dict : dict
        Dictionary with gene set names as keys and lists of genes as values

    Returns:
    --------
    dict
        Dictionary mapping cluster names to gene set scores
    """
    print("Calculating gene set scores for clusters...")

    # Calculate scores for each gene set
    for gene_set_name, genes in gene_sets_dict.items():
        print(f"  Scoring gene set: {gene_set_name}")

        # Filter genes to only those present in the dataset
        genes_in_data = [gene for gene in genes if gene in pseudobulk_adata.var.index]
        genes_missing = [gene for gene in genes if gene not in pseudobulk_adata.var.index]

        if len(genes_missing) > 0:
            print(f"    Warning: {len(genes_missing)} genes not found in dataset: {genes_missing[:5]}{'...' if len(genes_missing) > 5 else ''}")

        if len(genes_in_data) == 0:
            print(f"    Error: No genes from {gene_set_name} found in dataset. Skipping.")
            continue

        print(f"    Using {len(genes_in_data)}/{len(genes)} genes from {gene_set_name}")

        # Calculate gene set score using scanpy on pseudo-bulk data
        sc.tl.score_genes(pseudobulk_adata, gene_list=genes_in_data,
                         score_name=f"{gene_set_name}_score", use_raw=False)

    return pseudobulk_adata

def assign_cluster_scores_to_cells(adata, pseudobulk_adata, gene_sets_dict, cluster_key='leiden_new'):
    """
    Assign cluster-level gene set scores to individual cells.

    Parameters:
    -----------
    adata : AnnData
        Original single-cell AnnData object
    pseudobulk_adata : AnnData
        Pseudo-bulk AnnData with cluster scores
    gene_sets_dict : dict
        Dictionary with gene set names as keys
    cluster_key : str
        Key in adata.obs containing cluster assignments

    Returns:
    --------
    AnnData
        Updated AnnData object with cluster-level scores assigned to cells
    """
    print("Assigning cluster scores to individual cells...")

    adata_copy = adata.copy()

    # Create mapping from cluster to scores
    cluster_score_map = {}
    for gene_set_name in gene_sets_dict.keys():
        score_column = f"{gene_set_name}_score"
        if score_column not in pseudobulk_adata.obs.columns:
            continue

        cluster_score_map[score_column] = dict(zip(
            pseudobulk_adata.obs.index,
            pseudobulk_adata.obs[score_column]
        ))

    # Assign scores to cells based on their cluster
    for score_column, cluster_scores in cluster_score_map.items():
        # Map scores and ensure they are numeric (not categorical)
        mapped_scores = adata_copy.obs[cluster_key].map(cluster_scores)
        adata_copy.obs[f"cluster_{score_column}"] = pd.to_numeric(mapped_scores, errors='coerce')
        print(f"  Assigned {score_column} to cells")

    return adata_copy

def create_spatial_plots(adata, gene_sets_dict, output_dir, cluster_key='leiden_new'):
    """
    Create spatial plots showing cluster-level gene set scores with consistent color scales per time point.

    Parameters:
    -----------
    adata : AnnData
        AnnData object with cluster scores assigned to cells
    gene_sets_dict : dict
        Dictionary with gene set names
    output_dir : str
        Directory to save plots
    cluster_key : str
        Key in adata.obs containing cluster assignments
    """
    print("Creating spatial plots...")

    # Ensure spatial coordinates exist
    if 'spatial' not in adata.obsm:
        adata.obsm['spatial'] = np.c_[adata.obs['center_x'], adata.obs['center_y']]

    samples = sorted(adata.obs['sample_id'].unique())

    # Extract time points from sample IDs
    def get_time_point(sample_id):
        for tp in ['e16', 'p0', 'p7', 'p30']:
            if tp in str(sample_id).lower():
                return tp
        return 'unknown'

    # Group samples by time point
    time_point_groups = {}
    for sample_id in samples:
        tp = get_time_point(sample_id)
        if tp not in time_point_groups:
            time_point_groups[tp] = []
        time_point_groups[tp].append(sample_id)

    for gene_set_name in gene_sets_dict.keys():
        score_column = f"cluster_{gene_set_name}_score"

        # Skip if score column doesn't exist
        if score_column not in adata.obs.columns:
            print(f"Skipping {gene_set_name} - no cluster scores calculated")
            continue

        print(f"Creating spatial plots for gene set: {gene_set_name}")

        # Clean gene set name for filename
        clean_name = gene_set_name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')

        # Calculate color scales per time point
        time_point_scales = {}
        for time_point, tp_samples in time_point_groups.items():
            tp_adata = adata[adata.obs['sample_id'].isin(tp_samples)]
            tp_scores = tp_adata.obs[score_column].dropna()

            if len(tp_scores) > 0:
                vmin = tp_scores.quantile(0.01)
                vmax = tp_scores.quantile(0.99)
                time_point_scales[time_point] = (vmin, vmax)

        # Plot for each sample individually with consistent scales per time point
        for sample_id in samples:
            sample_adata = adata[adata.obs['sample_id'] == sample_id].copy()

            if len(sample_adata) == 0:
                continue

            time_point = get_time_point(sample_id)
            if time_point not in time_point_scales:
                continue

            vmin, vmax = time_point_scales[time_point]

            # Create figure
            fig, ax = plt.subplots(figsize=(12, 10))

            # Get spatial coordinates and scores
            spatial_coords = sample_adata.obsm['spatial']
            scores = sample_adata.obs[score_column]

            # Calculate color scale based on all non-NaN scores
            valid_scores = scores.dropna()
            if len(valid_scores) == 0:
                print(f"  Warning: No valid scores for {gene_set_name} in sample {sample_id}")
                plt.close(fig)
                continue

            # Create scatter plot with coolwarm colormap (diverging around 0)
            scatter = ax.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                               c=scores, cmap='coolwarm', s=5, alpha=0.8,
                               vmin=vmin, vmax=vmax)

            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
            cbar.set_label(f'{gene_set_name} Score\n(Cluster Level)', rotation=270, labelpad=20)

            # Set labels and title
            ax.set_xlabel('X coordinate')
            ax.set_ylabel('Y coordinate')
            ax.set_title(f'{gene_set_name} - Sample {sample_id}\n'
                       f'Cluster-level pseudo-bulk scores\n'
                       f'Scale: {vmin:.3f} to {vmax:.3f} ({time_point} range)')
            ax.set_aspect('equal')

            # Save the plot
            plot_filename = f"{clean_name}_cluster_pseudobulk_sample_{sample_id}.png"
            fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight', dpi=300)
            plt.close(fig)

            print(f"  Saved: {plot_filename}")

def create_cluster_score_summary(adata, pseudobulk_adata, gene_sets_dict, output_dir, cluster_key='leiden_new'):
    """
    Create summary table showing cluster-level scores.

    Parameters:
    -----------
    adata : AnnData
        Original single-cell data
    pseudobulk_adata : AnnData
        Pseudo-bulk data with cluster scores
    gene_sets_dict : dict
        Dictionary with gene set names
    output_dir : str
        Directory to save summary
    cluster_key : str
        Key for cluster assignments
    """
    print("Creating cluster score summary...")

    summary_data = []

    for cluster in pseudobulk_adata.obs.index:
        n_cells = len(adata[adata.obs[cluster_key] == cluster])

        row = {
            'Cluster': cluster,
            'N_Cells': n_cells
        }

        for gene_set_name in gene_sets_dict.keys():
            score_column = f"{gene_set_name}_score"
            if score_column in pseudobulk_adata.obs.columns:
                row[f"{gene_set_name}_score"] = pseudobulk_adata.obs.loc[cluster, score_column]

        summary_data.append(row)

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(output_dir, 'cluster_pseudobulk_scores_summary.csv'), index=False)

    print(f"  Saved: cluster_pseudobulk_scores_summary.csv")

    # Create a heatmap of cluster scores
    score_columns = [col for col in summary_df.columns if col.endswith('_score')]
    if len(score_columns) > 0:
        heatmap_data = summary_df.set_index('Cluster')[score_columns]

        # Rename columns to remove '_score' suffix for cleaner display
        heatmap_data.columns = [col.replace('_score', '') for col in heatmap_data.columns]

        # Create heatmap
        fig, ax = plt.subplots(figsize=(12, max(6, len(heatmap_data) * 0.3)))
        sns.heatmap(heatmap_data, cmap='coolwarm', center=0,
                   annot=True, fmt='.3f', linewidths=0.5,
                   cbar_kws={'label': 'Gene Set Score'}, ax=ax)

        ax.set_title('Cluster-Level Gene Set Scores (Pseudo-Bulk)')
        ax.set_xlabel('Gene Sets')
        ax.set_ylabel('Clusters')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'cluster_pseudobulk_scores_heatmap.png'),
                   bbox_inches='tight', dpi=300)
        plt.close()

        print(f"  Saved: cluster_pseudobulk_scores_heatmap.png")

def main():
    """Main function to run the cluster pseudo-bulk gene set scoring pipeline."""

    print("Cluster Pseudo-Bulk Gene Set Scoring and Spatial Visualization Pipeline")
    print("=" * 80)

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

    # Check for cluster annotations
    cluster_key = 'leiden_new'
    if cluster_key not in adata.obs.columns:
        print(f"Error: Cluster key '{cluster_key}' not found in data")
        sys.exit(1)

    print(f"Using cluster key: {cluster_key}")
    print(f"Number of clusters: {len(adata.obs[cluster_key].unique())}")

    # Create pseudo-bulk profiles by cluster
    pseudobulk_adata = create_pseudobulk_by_cluster(adata, cluster_key=cluster_key)

    # Calculate gene set scores for clusters
    pseudobulk_adata = calculate_cluster_gene_set_scores(pseudobulk_adata, gene_sets_dict)

    # Assign cluster scores to individual cells
    adata_with_cluster_scores = assign_cluster_scores_to_cells(
        adata, pseudobulk_adata, gene_sets_dict, cluster_key=cluster_key
    )

    # Create spatial plots
    create_spatial_plots(adata_with_cluster_scores, gene_sets_dict, output_dir, cluster_key=cluster_key)

    # Create summary tables and heatmap
    create_cluster_score_summary(adata, pseudobulk_adata, gene_sets_dict, output_dir, cluster_key=cluster_key)

    print("\n" + "=" * 80)
    print("Pipeline completed successfully!")
    print(f"Results saved to: {output_dir}")
    print("  - Spatial plots with cluster-level scores")
    print("  - Cluster score summary CSV")
    print("  - Cluster score heatmap")
    print("=" * 80)

if __name__ == "__main__":
    main()
