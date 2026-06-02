#!/usr/bin/env python3
"""
Solution 1: Global Expression Filtering

This approach filters out genes that are not expressed anywhere in the dataset
before calculating gene set scores. This reduces noise from genes that are
absent from the data, but uses the SAME filtered gene list for all clusters.

Key features:
- Removes genes with mean expression < threshold across entire dataset
- All clusters scored using the same filtered gene list
- Good for removing technical noise and missing genes
- Does NOT account for cell-type-specific expression patterns

Use case: When you want to remove non-expressed genes but maintain uniform
          scoring across all clusters for easy comparison.
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
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_1_global_filter/results'

# Create output directory
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

# Parameters
GLOBAL_EXPR_THRESHOLD = 0.01  # Mean expression threshold across all cells

def load_gene_sets(gene_sets_path):
    """Load gene sets from CSV file."""
    print(f"Loading gene sets from: {gene_sets_path}")

    gene_sets_df = pd.read_csv(gene_sets_path, header=None,
                              names=['gene_set', 'gene', 'col3', 'col4'])

    gene_sets_dict = {}
    for gene_set, group in gene_sets_df.groupby('gene_set'):
        gene_sets_dict[gene_set] = group['gene'].tolist()

    print(f"Loaded {len(gene_sets_dict)} gene sets:")
    for gene_set, genes in gene_sets_dict.items():
        print(f"  - {gene_set}: {len(genes)} genes")

    return gene_sets_dict

def filter_genes_by_global_expression(adata, threshold=0.01):
    """
    Identify genes that are expressed above threshold across the entire dataset.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    threshold : float
        Minimum mean expression across all cells

    Returns:
    --------
    set
        Set of gene names that pass the expression threshold
    """
    print(f"\nFiltering genes by global expression (threshold={threshold})...")

    # Calculate mean expression across all cells
    if hasattr(adata.X, 'toarray'):
        mean_expr = np.array(adata.X.mean(axis=0)).flatten()
    else:
        mean_expr = adata.X.mean(axis=0)

    # Identify expressed genes
    expressed_genes = set(adata.var.index[mean_expr > threshold])

    print(f"  Total genes in dataset: {adata.n_vars}")
    print(f"  Genes passing threshold: {len(expressed_genes)}")
    print(f"  Genes filtered out: {adata.n_vars - len(expressed_genes)}")

    return expressed_genes

def create_pseudobulk_by_cluster(adata, cluster_key='leiden_new'):
    """Create pseudo-bulk expression profiles for each cluster."""
    print(f"\nCreating pseudo-bulk profiles for clusters...")

    clusters = sorted(adata.obs[cluster_key].unique())
    print(f"  Found {len(clusters)} clusters: {clusters}")

    pseudobulk_data = []
    cluster_names = []

    for cluster in clusters:
        cluster_cells = adata[adata.obs[cluster_key] == cluster]
        mean_expr = cluster_cells.X.mean(axis=0)
        if hasattr(mean_expr, 'A1'):
            mean_expr = mean_expr.A1
        pseudobulk_data.append(mean_expr)
        cluster_names.append(cluster)

    pseudobulk_data = np.array(pseudobulk_data)
    pseudobulk_adata = sc.AnnData(
        X=pseudobulk_data,
        obs=pd.DataFrame(index=cluster_names, data={cluster_key: cluster_names}),
        var=adata.var.copy()
    )

    print(f"  Created pseudo-bulk AnnData: {pseudobulk_adata.n_obs} clusters x {pseudobulk_adata.n_vars} genes")

    return pseudobulk_adata

def calculate_cluster_gene_set_scores_with_global_filter(pseudobulk_adata, gene_sets_dict,
                                                          expressed_genes):
    """
    Calculate gene set scores using globally filtered gene lists.

    Parameters:
    -----------
    pseudobulk_adata : AnnData
        Pseudo-bulk data
    gene_sets_dict : dict
        Dictionary of gene sets
    expressed_genes : set
        Set of genes that pass global expression threshold

    Returns:
    --------
    tuple
        (pseudobulk_adata, filtering_stats)
    """
    print("\nCalculating gene set scores with global filtering...")

    filtering_stats = []

    for gene_set_name, genes in gene_sets_dict.items():
        print(f"\n  Scoring gene set: {gene_set_name}")

        # Filter 1: Genes in dataset
        genes_in_data = [gene for gene in genes if gene in pseudobulk_adata.var.index]

        # Filter 2: Genes passing global expression threshold
        genes_globally_expressed = [gene for gene in genes_in_data if gene in expressed_genes]

        genes_filtered_by_expr = set(genes_in_data) - set(genes_globally_expressed)
        genes_missing = set(genes) - set(genes_in_data)

        print(f"    Total genes in set: {len(genes)}")
        print(f"    Genes in dataset: {len(genes_in_data)}")
        print(f"    Genes passing global expression filter: {len(genes_globally_expressed)}")
        print(f"    Genes filtered out (low/no expression): {len(genes_filtered_by_expr)}")

        # Store statistics
        filtering_stats.append({
            'Gene_Set': gene_set_name,
            'Total_Genes': len(genes),
            'Genes_In_Dataset': len(genes_in_data),
            'Genes_After_Global_Filter': len(genes_globally_expressed),
            'Genes_Filtered_Low_Expr': len(genes_filtered_by_expr),
            'Genes_Missing': len(genes_missing),
            'Filtering_Rate': len(genes_filtered_by_expr) / len(genes) if len(genes) > 0 else 0
        })

        if len(genes_globally_expressed) == 0:
            print(f"    WARNING: No genes passed global expression filter. Skipping.")
            continue

        # Calculate gene set score using filtered gene list
        sc.tl.score_genes(pseudobulk_adata, gene_list=genes_globally_expressed,
                         score_name=f"{gene_set_name}_score", use_raw=False)

    return pseudobulk_adata, filtering_stats

def assign_cluster_scores_to_cells(adata, pseudobulk_adata, gene_sets_dict, cluster_key='leiden_new'):
    """Assign cluster-level gene set scores to individual cells."""
    print("\nAssigning cluster scores to individual cells...")

    adata_copy = adata.copy()

    cluster_score_map = {}
    for gene_set_name in gene_sets_dict.keys():
        score_column = f"{gene_set_name}_score"
        if score_column not in pseudobulk_adata.obs.columns:
            continue

        cluster_score_map[score_column] = dict(zip(
            pseudobulk_adata.obs.index,
            pseudobulk_adata.obs[score_column]
        ))

    for score_column, cluster_scores in cluster_score_map.items():
        mapped_scores = adata_copy.obs[cluster_key].map(cluster_scores)
        adata_copy.obs[f"cluster_{score_column}"] = pd.to_numeric(mapped_scores, errors='coerce')
        print(f"  Assigned {score_column} to cells")

    return adata_copy

def create_spatial_plots(adata, gene_sets_dict, output_dir, cluster_key='leiden_new'):
    """Create spatial plots showing cluster-level gene set scores with consistent color scales per time point."""
    print("\nCreating spatial plots...")

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

        if score_column not in adata.obs.columns:
            print(f"Skipping {gene_set_name} - no cluster scores calculated")
            continue

        print(f"Creating spatial plots for gene set: {gene_set_name}")

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

        for sample_id in samples:
            sample_adata = adata[adata.obs['sample_id'] == sample_id].copy()

            if len(sample_adata) == 0:
                continue

            time_point = get_time_point(sample_id)
            if time_point not in time_point_scales:
                continue

            vmin, vmax = time_point_scales[time_point]

            fig, ax = plt.subplots(figsize=(12, 10))

            spatial_coords = sample_adata.obsm['spatial']
            scores = sample_adata.obs[score_column]

            valid_scores = scores.dropna()
            if len(valid_scores) == 0:
                plt.close(fig)
                continue

            scatter = ax.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                               c=scores, cmap='coolwarm', s=5, alpha=0.8,
                               vmin=vmin, vmax=vmax)

            cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
            cbar.set_label(f'{gene_set_name} Score\n(Global Filter)', rotation=270, labelpad=20)

            ax.set_xlabel('X coordinate')
            ax.set_ylabel('Y coordinate')
            ax.set_title(f'{gene_set_name} - Sample {sample_id}\n'
                       f'Global expression filtering (threshold={GLOBAL_EXPR_THRESHOLD})\n'
                       f'Scale: {vmin:.3f} to {vmax:.3f} ({time_point} range)')
            ax.set_aspect('equal')

            plot_filename = f"{clean_name}_global_filter_sample_{sample_id}.png"
            fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight', dpi=300)
            plt.close(fig)

            print(f"  Saved: {plot_filename}")

def create_cluster_score_heatmap(pseudobulk_adata, gene_sets_dict, output_dir, cluster_key='leiden_new'):
    """Create heatmap of cluster scores."""
    print("\nCreating cluster score heatmap...")

    score_columns = [f"{gene_set}_score" for gene_set in gene_sets_dict.keys()
                     if f"{gene_set}_score" in pseudobulk_adata.obs.columns]

    if len(score_columns) == 0:
        print("  No scores to plot")
        return

    heatmap_data = pseudobulk_adata.obs[score_columns].copy()
    heatmap_data.columns = [col.replace('_score', '') for col in heatmap_data.columns]

    fig, ax = plt.subplots(figsize=(12, max(6, len(heatmap_data) * 0.3)))
    sns.heatmap(heatmap_data, cmap='coolwarm', center=0,
               annot=True, fmt='.3f', linewidths=0.5,
               cbar_kws={'label': 'Gene Set Score'}, ax=ax)

    ax.set_title('Cluster-Level Gene Set Scores\n(Global Expression Filtering)')
    ax.set_xlabel('Gene Sets')
    ax.set_ylabel('Clusters')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cluster_scores_heatmap_global_filter.png'),
               bbox_inches='tight', dpi=300)
    plt.close()

    print(f"  Saved: cluster_scores_heatmap_global_filter.png")

def save_filtering_stats(filtering_stats, output_dir):
    """Save filtering statistics to CSV."""
    print("\nSaving filtering statistics...")

    stats_df = pd.DataFrame(filtering_stats)
    stats_df.to_csv(os.path.join(output_dir, 'global_filtering_statistics.csv'), index=False)

    print(f"  Saved: global_filtering_statistics.csv")
    print("\nFiltering Summary:")
    print(stats_df.to_string(index=False))

def main():
    """Main function."""

    print("=" * 80)
    print("SOLUTION 1: GLOBAL EXPRESSION FILTERING")
    print("=" * 80)

    # Load data
    print(f"\nLoading merged data from: {merged_h5ad_path}")
    adata = sc.read_h5ad(merged_h5ad_path)
    print(f"Loaded data: {adata.n_obs} cells x {adata.n_vars} genes")

    # Load gene sets
    gene_sets_dict = load_gene_sets(gene_sets_path)

    # Filter genes by global expression
    expressed_genes = filter_genes_by_global_expression(adata, threshold=GLOBAL_EXPR_THRESHOLD)

    # Create pseudo-bulk profiles
    cluster_key = 'leiden_new'
    pseudobulk_adata = create_pseudobulk_by_cluster(adata, cluster_key=cluster_key)

    # Calculate scores with global filtering
    pseudobulk_adata, filtering_stats = calculate_cluster_gene_set_scores_with_global_filter(
        pseudobulk_adata, gene_sets_dict, expressed_genes
    )

    # Save filtering statistics
    save_filtering_stats(filtering_stats, output_dir)

    # Assign scores to cells
    adata_with_scores = assign_cluster_scores_to_cells(
        adata, pseudobulk_adata, gene_sets_dict, cluster_key=cluster_key
    )

    # Create visualizations
    create_spatial_plots(adata_with_scores, gene_sets_dict, output_dir, cluster_key=cluster_key)
    create_cluster_score_heatmap(pseudobulk_adata, gene_sets_dict, output_dir, cluster_key=cluster_key)

    print("\n" + "=" * 80)
    print("SOLUTION 1 COMPLETED")
    print(f"Results saved to: {output_dir}")
    print("=" * 80)

if __name__ == "__main__":
    main()
