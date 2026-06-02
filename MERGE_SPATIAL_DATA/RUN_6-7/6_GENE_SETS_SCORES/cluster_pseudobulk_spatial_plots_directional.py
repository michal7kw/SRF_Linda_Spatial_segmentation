#!/usr/bin/env python3
"""
Cluster-Level Pseudo-Bulk Gene Set Scoring with Directional Information

This script calculates gene set scores at the cluster level (pseudo-bulk approach)
using directional gene lists where each gene is marked as positively or negatively
correlated with a biological process. Three types of scores are calculated:
1. Positive scores (using only upregulated genes)
2. Negative scores (using only downregulated genes)
3. Aggregate scores (combining both with correct signs)

The workflow:
1. Load spatial transcriptomics data
2. Load gene sets with directional information
3. Calculate pseudo-bulk expression for each cluster
4. Score gene sets on pseudo-bulk cluster profiles (positive, negative, aggregate)
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
gene_sets_path = 'DATA/GENE_LISTS/directional_lists/gene_lists.csv'
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/cluster_pseudobulk_spatial_plots_directional'

# Create output directory
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

def load_directional_gene_sets(gene_sets_path):
    """
    Load gene sets from CSV file with directional information.

    The CSV format is:
    - Row 1: Gene set names (biological processes)
    - Row 2: Direction ('up' or 'down') indicating correlation with the process
    - Rows 3+: Gene symbols

    Parameters:
    -----------
    gene_sets_path : str
        Path to the gene sets CSV file

    Returns:
    --------
    dict
        Dictionary with structure:
        {
            'gene_set_name': {
                'up': [list of upregulated genes],
                'down': [list of downregulated genes],
                'all_genes': [all genes],
                'direction': 'up' or 'down' (from row 2)
            }
        }
    """
    print(f"Loading directional gene sets from: {gene_sets_path}")

    # Read the CSV file
    df = pd.read_csv(gene_sets_path, header=None)

    # First row contains gene set names, second row contains directions
    gene_set_names = df.iloc[0].values
    directions = df.iloc[1].values

    # Group columns by gene set name
    gene_sets_dict = {}

    for col_idx in range(len(gene_set_names)):
        gene_set_name = gene_set_names[col_idx]
        direction = directions[col_idx]

        # Skip if gene set name or direction is NaN
        if pd.isna(gene_set_name) or pd.isna(direction):
            continue

        # Extract genes (rows 2 onwards), excluding NaN values
        genes = df.iloc[2:, col_idx].dropna().tolist()

        # Remove empty strings
        genes = [g for g in genes if g and str(g).strip()]

        if not genes:
            continue

        # Initialize gene set entry if not exists
        if gene_set_name not in gene_sets_dict:
            gene_sets_dict[gene_set_name] = {
                'up': [],
                'down': [],
                'all_genes': []
            }

        # Add genes to appropriate list
        if direction == 'up':
            gene_sets_dict[gene_set_name]['up'].extend(genes)
        elif direction == 'down':
            gene_sets_dict[gene_set_name]['down'].extend(genes)

        gene_sets_dict[gene_set_name]['all_genes'].extend(genes)

    print(f"Loaded {len(gene_sets_dict)} gene sets:")
    for gene_set_name, gene_info in gene_sets_dict.items():
        n_up = len(gene_info['up'])
        n_down = len(gene_info['down'])
        n_total = len(gene_info['all_genes'])
        print(f"  - {gene_set_name}: {n_total} genes ({n_up} up, {n_down} down)")

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

def calculate_cluster_gene_set_scores_directional(pseudobulk_adata, gene_sets_dict):
    """
    Calculate directional gene set scores for each cluster using pseudo-bulk profiles.

    For each gene set, calculates:
    - Positive score (using only upregulated genes)
    - Negative score (using only downregulated genes)
    - Aggregate score (upregulated genes positively, downregulated genes negatively)

    Parameters:
    -----------
    pseudobulk_adata : AnnData
        Pseudo-bulk AnnData object with one observation per cluster
    gene_sets_dict : dict
        Dictionary with gene set information including directional genes

    Returns:
    --------
    AnnData
        Updated pseudo-bulk AnnData with directional scores
    """
    print("Calculating directional gene set scores for clusters...")

    for gene_set_name, gene_info in gene_sets_dict.items():
        print(f"  Scoring gene set: {gene_set_name}")

        # Process upregulated genes
        up_genes = gene_info['up']
        up_genes_in_data = [gene for gene in up_genes if gene in pseudobulk_adata.var.index]

        if len(up_genes_in_data) > 0:
            print(f"    Positive: Using {len(up_genes_in_data)}/{len(up_genes)} upregulated genes")
            sc.tl.score_genes(pseudobulk_adata, gene_list=up_genes_in_data,
                            score_name=f"{gene_set_name}_positive_score", use_raw=False)
        else:
            print(f"    Warning: No upregulated genes found in dataset")
            pseudobulk_adata.obs[f"{gene_set_name}_positive_score"] = 0.0

        # Process downregulated genes
        down_genes = gene_info['down']
        down_genes_in_data = [gene for gene in down_genes if gene in pseudobulk_adata.var.index]

        if len(down_genes_in_data) > 0:
            print(f"    Negative: Using {len(down_genes_in_data)}/{len(down_genes)} downregulated genes")
            sc.tl.score_genes(pseudobulk_adata, gene_list=down_genes_in_data,
                            score_name=f"{gene_set_name}_negative_score", use_raw=False)
        else:
            print(f"    Warning: No downregulated genes found in dataset")
            pseudobulk_adata.obs[f"{gene_set_name}_negative_score"] = 0.0

        # Calculate aggregate score
        # Aggregate = positive contribution from upregulated - positive contribution from downregulated
        # (high expression of up genes is good, high expression of down genes is bad)
        if f"{gene_set_name}_positive_score" in pseudobulk_adata.obs.columns and \
           f"{gene_set_name}_negative_score" in pseudobulk_adata.obs.columns:
            pseudobulk_adata.obs[f"{gene_set_name}_aggregate_score"] = \
                pseudobulk_adata.obs[f"{gene_set_name}_positive_score"] - \
                pseudobulk_adata.obs[f"{gene_set_name}_negative_score"]
            print(f"    Aggregate score calculated")

    return pseudobulk_adata

def assign_cluster_scores_to_cells(adata, pseudobulk_adata, gene_sets_dict, cluster_key='leiden_new'):
    """
    Assign cluster-level directional gene set scores to individual cells.

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

    # Create mapping from cluster to scores for all score types
    for gene_set_name in gene_sets_dict.keys():
        for score_type in ['positive', 'negative', 'aggregate']:
            score_column = f"{gene_set_name}_{score_type}_score"

            if score_column not in pseudobulk_adata.obs.columns:
                continue

            cluster_score_map = dict(zip(
                pseudobulk_adata.obs.index,
                pseudobulk_adata.obs[score_column]
            ))

            # Map scores to cells based on their cluster
            mapped_scores = adata_copy.obs[cluster_key].map(cluster_score_map)
            adata_copy.obs[f"cluster_{score_column}"] = pd.to_numeric(mapped_scores, errors='coerce')
            print(f"  Assigned {score_column} to cells")

    return adata_copy

def create_spatial_plots(adata, gene_sets_dict, output_dir, cluster_key='leiden_new'):
    """
    Create spatial plots showing cluster-level directional gene set scores with consistent color scales per time point.

    Creates separate plots for positive, negative, and aggregate scores.

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
    cmaps = {'positive': 'Reds', 'negative': 'Blues', 'aggregate': 'coolwarm'}

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
        print(f"Creating spatial plots for gene set: {gene_set_name}")

        # Clean gene set name for filename
        clean_name = gene_set_name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')

        # Create plots for each score type
        for score_type in ['positive', 'negative', 'aggregate']:
            score_column = f"cluster_{gene_set_name}_{score_type}_score"

            # Skip if score column doesn't exist
            if score_column not in adata.obs.columns:
                continue

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
                    plt.close(fig)
                    continue

                # Create scatter plot
                scatter = ax.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                                   c=scores, cmap=cmaps[score_type], s=5, alpha=0.8,
                                   vmin=vmin, vmax=vmax)

                # Add colorbar
                cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
                cbar.set_label(f'{gene_set_name} {score_type.capitalize()} Score\n(Cluster Level)',
                             rotation=270, labelpad=20)

                # Set labels and title
                ax.set_xlabel('X coordinate')
                ax.set_ylabel('Y coordinate')
                ax.set_title(f'{gene_set_name} - {score_type.capitalize()} Score\n'
                           f'Sample {sample_id} - Cluster-level pseudo-bulk scores\n'
                           f'Scale: {vmin:.3f} to {vmax:.3f} ({time_point} range)')
                ax.set_aspect('equal')

                # Save the plot
                plot_filename = f"{clean_name}_{score_type}_cluster_pseudobulk_sample_{sample_id}.png"
                fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight', dpi=300)
                plt.close(fig)

            print(f"  Saved {score_type} score plots")

def create_cluster_score_summary(adata, pseudobulk_adata, gene_sets_dict, output_dir, cluster_key='leiden_new'):
    """
    Create summary table showing cluster-level directional scores.

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

    # Create summary for each score type
    for score_type in ['positive', 'negative', 'aggregate']:
        summary_data = []

        for cluster in pseudobulk_adata.obs.index:
            n_cells = len(adata[adata.obs[cluster_key] == cluster])

            row = {
                'Cluster': cluster,
                'N_Cells': n_cells
            }

            for gene_set_name in gene_sets_dict.keys():
                score_column = f"{gene_set_name}_{score_type}_score"
                if score_column in pseudobulk_adata.obs.columns:
                    row[f"{gene_set_name}_{score_type}"] = pseudobulk_adata.obs.loc[cluster, score_column]

            summary_data.append(row)

        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(os.path.join(output_dir, f'cluster_pseudobulk_{score_type}_scores_summary.csv'),
                         index=False)
        print(f"  Saved: cluster_pseudobulk_{score_type}_scores_summary.csv")

        # Create a heatmap for this score type
        score_columns = [col for col in summary_df.columns
                        if col.endswith(f'_{score_type}') and col not in ['Cluster', 'N_Cells']]

        if len(score_columns) > 0:
            heatmap_data = summary_df.set_index('Cluster')[score_columns]

            # Rename columns to remove score type suffix for cleaner display
            heatmap_data.columns = [col.replace(f'_{score_type}', '') for col in heatmap_data.columns]

            # Create heatmap
            fig, ax = plt.subplots(figsize=(12, max(6, len(heatmap_data) * 0.3)))

            # Use different colormaps for different score types
            if score_type == 'positive':
                cmap = 'Reds'
            elif score_type == 'negative':
                cmap = 'Blues'
            else:  # aggregate
                cmap = 'coolwarm'

            center_val = 0 if score_type == 'aggregate' else None

            sns.heatmap(heatmap_data, cmap=cmap, center=center_val,
                       annot=True, fmt='.3f', linewidths=0.5,
                       cbar_kws={'label': f'{score_type.capitalize()} Gene Set Score'}, ax=ax)

            ax.set_title(f'Cluster-Level {score_type.capitalize()} Gene Set Scores (Pseudo-Bulk)')
            ax.set_xlabel('Gene Sets')
            ax.set_ylabel('Clusters')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'cluster_pseudobulk_{score_type}_scores_heatmap.png'),
                       bbox_inches='tight', dpi=300)
            plt.close()

            print(f"  Saved: cluster_pseudobulk_{score_type}_scores_heatmap.png")

def main():
    """Main function to run the directional cluster pseudo-bulk gene set scoring pipeline."""

    print("Directional Cluster Pseudo-Bulk Gene Set Scoring Pipeline")
    print("=" * 80)

    # Load the merged AnnData object
    if os.path.exists(merged_h5ad_path):
        print(f"Loading merged data from: {merged_h5ad_path}")
        adata = sc.read_h5ad(merged_h5ad_path)
        print(f"Loaded data: {adata.n_obs} cells x {adata.n_vars} genes")
    else:
        print(f"Error: Merged data file not found at {merged_h5ad_path}")
        sys.exit(1)

    # Load directional gene sets
    if os.path.exists(gene_sets_path):
        gene_sets_dict = load_directional_gene_sets(gene_sets_path)
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

    # Calculate directional gene set scores for clusters
    pseudobulk_adata = calculate_cluster_gene_set_scores_directional(pseudobulk_adata, gene_sets_dict)

    # Assign cluster scores to individual cells
    adata_with_cluster_scores = assign_cluster_scores_to_cells(
        adata, pseudobulk_adata, gene_sets_dict, cluster_key=cluster_key
    )

    # Create spatial plots
    create_spatial_plots(adata_with_cluster_scores, gene_sets_dict, output_dir, cluster_key=cluster_key)

    # Create summary tables and heatmaps
    create_cluster_score_summary(adata, pseudobulk_adata, gene_sets_dict, output_dir, cluster_key=cluster_key)

    print("\n" + "=" * 80)
    print("Pipeline completed successfully!")
    print(f"Results saved to: {output_dir}")
    print("  - Spatial plots with cluster-level directional scores (positive, negative, aggregate)")
    print("  - Cluster score summary CSVs for each score type")
    print("  - Cluster score heatmaps for each score type")
    print("=" * 80)

if __name__ == "__main__":
    main()
