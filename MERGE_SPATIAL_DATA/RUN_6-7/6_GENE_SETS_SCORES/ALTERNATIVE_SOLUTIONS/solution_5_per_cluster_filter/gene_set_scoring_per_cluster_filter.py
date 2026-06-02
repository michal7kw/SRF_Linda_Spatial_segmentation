#!/usr/bin/env python3
"""
Solution 5: Per-Cluster Gene Filtering

This is THE SOLUTION for the dilution problem raised in your Reelin example!

For each cluster, genes are filtered based on whether they are expressed in THAT
specific cluster before scoring. Different clusters get scored on different gene
subsets from the same gene set.

Example: Reelin pathway (20 genes)
- Granule cells: Scored on 3 receptor genes (only ones expressed in GCs)
- Cajal-Retzius: Scored on 1 ligand gene (only one expressed in CR cells)
- Other clusters: Scored on their relevant subset

Key features:
- Cluster-specific gene filtering
- Each cluster scored on biologically relevant genes
- Directly addresses cell-type-specific pathway components
- Maintains standard scoring method (scanpy.tl.score_genes)

Use case: THE solution when pathways have modular, cell-type-specific components
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
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_5_per_cluster_filter/results'

# Create output directory
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

# Parameters
CLUSTER_EXPR_THRESHOLD = 0.5  # Expression in pseudobulk to consider gene expressed
MIN_CLUSTER_FRACTION = 0.1  # Minimum fraction of cells in cluster where gene must be expressed
MIN_GENES_FOR_SCORE = 2  # Minimum genes needed for valid score

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

def filter_genes_per_cluster(pseudobulk_adata, genes, expr_threshold=0.5):
    """
    For each cluster, identify which genes from the gene set are expressed.

    Returns dictionary mapping cluster -> list of expressed genes
    """
    genes_in_data = [g for g in genes if g in pseudobulk_adata.var.index]

    if len(genes_in_data) == 0:
        return {}

    gene_indices = [pseudobulk_adata.var.index.get_loc(g) for g in genes_in_data]

    cluster_gene_map = {}

    for i, cluster in enumerate(pseudobulk_adata.obs.index):
        cluster_expr = pseudobulk_adata.X[i, gene_indices]

        # Filter to expressed genes in this cluster
        expressed_mask = cluster_expr > expr_threshold
        expressed_genes = [genes_in_data[j] for j, is_expr in enumerate(expressed_mask) if is_expr]

        cluster_gene_map[cluster] = expressed_genes

    return cluster_gene_map

def calculate_per_cluster_filtered_scores(pseudobulk_adata, adata, gene_sets_dict,
                                          cluster_key='leiden_new',
                                          expr_threshold=0.5,
                                          min_genes=2):
    """
    Calculate gene set scores with per-cluster gene filtering.

    For each cluster:
    1. Identify which genes are expressed in that cluster
    2. Score only those genes for cells in that cluster
    3. Different clusters may be scored on different gene subsets

    Returns both pseudobulk scores and single-cell scores
    """
    print(f"\nCalculating per-cluster filtered gene set scores...")
    print(f"  Cluster expression threshold: {expr_threshold}")
    print(f"  Minimum genes for scoring: {min_genes}")

    filtering_stats = []
    cluster_gene_usage = {}

    # Create copy of adata for scoring
    adata_scored = adata.copy()

    for gene_set_name, genes in gene_sets_dict.items():
        print(f"\n  Processing gene set: {gene_set_name}")
        print(f"    Total genes in set: {len(genes)}")

        # Get cluster-specific gene lists
        cluster_gene_map = filter_genes_per_cluster(pseudobulk_adata, genes, expr_threshold)

        # Store for reporting
        cluster_gene_usage[gene_set_name] = cluster_gene_map

        # Score each cluster separately
        cluster_scores = {}
        cluster_n_genes = {}

        for cluster, cluster_genes in cluster_gene_map.items():
            n_genes = len(cluster_genes)

            print(f"      {cluster}: {n_genes}/{len([g for g in genes if g in pseudobulk_adata.var.index])} genes expressed")

            if n_genes < min_genes:
                print(f"        WARNING: Only {n_genes} genes - score will be NaN")
                cluster_scores[cluster] = np.nan
                cluster_n_genes[cluster] = n_genes
                continue

            # Get cells from this cluster
            cluster_cells = adata_scored[adata_scored.obs[cluster_key] == cluster].copy()

            if len(cluster_cells) == 0:
                cluster_scores[cluster] = np.nan
                cluster_n_genes[cluster] = 0
                continue

            # Score this cluster with its specific gene list
            sc.tl.score_genes(cluster_cells, gene_list=cluster_genes,
                            score_name=f"{gene_set_name}_score", use_raw=False)

            # Get mean score for this cluster (for pseudobulk)
            cluster_scores[cluster] = cluster_cells.obs[f"{gene_set_name}_score"].mean()
            cluster_n_genes[cluster] = n_genes

            # Assign scores back to main adata
            cell_indices = adata_scored.obs[cluster_key] == cluster
            adata_scored.obs.loc[cell_indices, f"{gene_set_name}_score"] = \
                cluster_cells.obs[f"{gene_set_name}_score"].values

            # Store filtering stats
            filtering_stats.append({
                'Gene_Set': gene_set_name,
                'Cluster': cluster,
                'N_Genes_Total': len([g for g in genes if g in pseudobulk_adata.var.index]),
                'N_Genes_Used': n_genes,
                'Frac_Genes_Used': n_genes / len([g for g in genes if g in pseudobulk_adata.var.index]) if len([g for g in genes if g in pseudobulk_adata.var.index]) > 0 else 0,
                'Genes_Used': ', '.join(cluster_genes),
                'Score': cluster_scores[cluster],
                'Expression_Threshold': expr_threshold
            })

        # Add cluster-level metrics to pseudobulk
        pseudobulk_adata.obs[f"{gene_set_name}_score"] = [cluster_scores.get(c, np.nan) for c in pseudobulk_adata.obs.index]
        pseudobulk_adata.obs[f"{gene_set_name}_n_genes"] = [cluster_n_genes.get(c, 0) for c in pseudobulk_adata.obs.index]

        # Add n_genes to single-cell data
        for cluster in pseudobulk_adata.obs.index:
            cell_indices = adata_scored.obs[cluster_key] == cluster
            adata_scored.obs.loc[cell_indices, f"{gene_set_name}_n_genes"] = cluster_n_genes.get(cluster, 0)

    return adata_scored, pseudobulk_adata, filtering_stats, cluster_gene_usage

def create_spatial_plots(adata, gene_sets_dict, output_dir):
    """Create spatial plots with per-cluster filtered scores and consistent color scales per time point."""
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
        score_column = f"{gene_set_name}_score"
        n_genes_column = f"{gene_set_name}_n_genes"

        if score_column not in adata.obs.columns:
            print(f"Skipping {gene_set_name} - no scores calculated")
            continue

        print(f"Creating spatial plots for gene set: {gene_set_name}")

        clean_name = gene_set_name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')

        # Calculate color scales per time point for scores
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

            vmin_score, vmax_score = time_point_scales[time_point]

            # Create figure with two subplots
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 10))

            spatial_coords = sample_adata.obsm['spatial']
            scores = sample_adata.obs[score_column]
            n_genes = sample_adata.obs[n_genes_column]

            # Plot 1: Per-cluster filtered score with consistent scale
            valid_scores = scores.dropna()
            if len(valid_scores) > 0:
                scatter1 = ax1.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                                      c=scores, cmap='coolwarm', s=5, alpha=0.8,
                                      vmin=vmin_score, vmax=vmax_score)

                cbar1 = plt.colorbar(scatter1, ax=ax1, shrink=0.8)
                cbar1.set_label(f'{gene_set_name} Score\n(Per-Cluster Filtered)', rotation=270, labelpad=20)

                ax1.set_xlabel('X coordinate')
                ax1.set_ylabel('Y coordinate')
                ax1.set_title(f'Per-Cluster Filtered Score\n{gene_set_name}\nScale: {vmin_score:.3f} to {vmax_score:.3f} ({time_point} range)')
                ax1.set_aspect('equal')

            # Plot 2: Number of genes used per cluster
            valid_n_genes = n_genes.dropna()
            if len(valid_n_genes) > 0:
                scatter2 = ax2.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                                      c=n_genes, cmap='viridis', s=5, alpha=0.8,
                                      vmin=0, vmax=valid_n_genes.max())

                cbar2 = plt.colorbar(scatter2, ax=ax2, shrink=0.8)
                cbar2.set_label('N Genes in Cluster', rotation=270, labelpad=20)

                ax2.set_xlabel('X coordinate')
                ax2.set_ylabel('Y coordinate')
                ax2.set_title(f'Genes Used Per Cluster\n{gene_set_name}\nRange: {int(valid_n_genes.min())} to {int(valid_n_genes.max())}')
                ax2.set_aspect('equal')

            fig.suptitle(f'Sample {sample_id}', fontsize=16, y=0.98)

            plot_filename = f"{clean_name}_per_cluster_filter_sample_{sample_id}.png"
            fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight', dpi=300)
            plt.close(fig)

            print(f"  Saved: {plot_filename}")

def create_gene_usage_visualization(cluster_gene_usage, output_dir):
    """Create visualization showing which genes are used in which clusters."""
    print("\nCreating gene usage visualization...")

    for gene_set_name, cluster_genes in cluster_gene_usage.items():
        # Create binary matrix: clusters x genes
        all_genes = set()
        for genes in cluster_genes.values():
            all_genes.update(genes)

        all_genes = sorted(all_genes)
        clusters = sorted(cluster_genes.keys())

        if len(all_genes) == 0:
            continue

        # Create binary matrix
        matrix = np.zeros((len(clusters), len(all_genes)))
        for i, cluster in enumerate(clusters):
            for j, gene in enumerate(all_genes):
                if gene in cluster_genes[cluster]:
                    matrix[i, j] = 1

        # Create heatmap
        fig, ax = plt.subplots(figsize=(max(12, len(all_genes) * 0.5), max(6, len(clusters) * 0.4)))

        sns.heatmap(matrix, cmap='RdYlGn', vmin=0, vmax=1,
                   xticklabels=all_genes, yticklabels=clusters,
                   cbar_kws={'label': 'Gene Used (1) / Not Used (0)'},
                   linewidths=0.5, ax=ax)

        ax.set_title(f'Gene Usage Per Cluster\n{gene_set_name}\n(Which genes are scored in which clusters)')
        ax.set_xlabel('Genes')
        ax.set_ylabel('Clusters')

        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()

        clean_name = gene_set_name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
        plt.savefig(os.path.join(output_dir, f'gene_usage_{clean_name}.png'),
                   bbox_inches='tight', dpi=300)
        plt.close()

        print(f"  Saved: gene_usage_{clean_name}.png")

def create_score_heatmaps(pseudobulk_adata, gene_sets_dict, output_dir):
    """Create heatmaps of scores and N genes used."""
    print("\nCreating score heatmaps...")

    # Heatmap 1: Scores
    score_columns = [f"{gs}_score" for gs in gene_sets_dict.keys()
                     if f"{gs}_score" in pseudobulk_adata.obs.columns]

    if score_columns:
        score_data = pseudobulk_adata.obs[score_columns].copy()
        score_data.columns = [c.replace('_score', '') for c in score_data.columns]

        fig, ax = plt.subplots(figsize=(12, max(6, len(score_data) * 0.3)))
        sns.heatmap(score_data, cmap='coolwarm', center=0,
                   annot=True, fmt='.2f', linewidths=0.5,
                   cbar_kws={'label': 'Gene Set Score'}, ax=ax)

        ax.set_title('Per-Cluster Filtered Gene Set Scores')
        ax.set_xlabel('Gene Sets')
        ax.set_ylabel('Clusters')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'per_cluster_scores_heatmap.png'),
                   bbox_inches='tight', dpi=300)
        plt.close()

        print(f"  Saved: per_cluster_scores_heatmap.png")

    # Heatmap 2: N genes
    n_genes_columns = [f"{gs}_n_genes" for gs in gene_sets_dict.keys()
                       if f"{gs}_n_genes" in pseudobulk_adata.obs.columns]

    if n_genes_columns:
        n_genes_data = pseudobulk_adata.obs[n_genes_columns].copy()
        n_genes_data.columns = [c.replace('_n_genes', '') for c in n_genes_data.columns]

        fig, ax = plt.subplots(figsize=(12, max(6, len(n_genes_data) * 0.3)))
        sns.heatmap(n_genes_data, cmap='YlGnBu',
                   annot=True, fmt='.0f', linewidths=0.5,
                   cbar_kws={'label': 'N Genes Used'}, ax=ax)

        ax.set_title('Number of Genes Used Per Cluster')
        ax.set_xlabel('Gene Sets')
        ax.set_ylabel('Clusters')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'per_cluster_n_genes_heatmap.png'),
                   bbox_inches='tight', dpi=300)
        plt.close()

        print(f"  Saved: per_cluster_n_genes_heatmap.png")

def save_filtering_stats(filtering_stats, output_dir):
    """Save filtering statistics."""
    print("\nSaving filtering statistics...")

    df = pd.DataFrame(filtering_stats)
    df.to_csv(os.path.join(output_dir, 'per_cluster_filtering_statistics.csv'), index=False)

    print(f"  Saved: per_cluster_filtering_statistics.csv")

    # Print summary
    print("\nPer-Cluster Filtering Summary:")
    df_valid = df[df['Score'].notna()]
    if not df_valid.empty:
        summary = df_valid.groupby('Gene_Set').agg({
            'N_Genes_Used': ['mean', 'std', 'min', 'max'],
            'Frac_Genes_Used': ['mean', 'min', 'max'],
            'Score': ['mean', 'std', 'min', 'max']
        }).round(3)
        print(summary)

def main():
    """Main function."""

    print("=" * 80)
    print("SOLUTION 5: PER-CLUSTER GENE FILTERING")
    print("=" * 80)
    print("\nTHIS IS THE SOLUTION FOR THE REELIN DILUTION PROBLEM!")
    print("Each cluster is scored on only the genes it actually expresses.")
    print("=" * 80)

    # Load data
    print(f"\nLoading merged data from: {merged_h5ad_path}")
    adata = sc.read_h5ad(merged_h5ad_path)
    print(f"Loaded data: {adata.n_obs} cells x {adata.n_vars} genes")

    # Load gene sets
    gene_sets_dict = load_gene_sets(gene_sets_path)

    # Create pseudo-bulk profiles
    cluster_key = 'leiden_new'
    pseudobulk_adata = create_pseudobulk_by_cluster(adata, cluster_key=cluster_key)

    # Calculate per-cluster filtered scores
    adata_scored, pseudobulk_scored, filtering_stats, cluster_gene_usage = \
        calculate_per_cluster_filtered_scores(
            pseudobulk_adata, adata, gene_sets_dict,
            cluster_key=cluster_key,
            expr_threshold=CLUSTER_EXPR_THRESHOLD,
            min_genes=MIN_GENES_FOR_SCORE
        )

    # Save statistics
    save_filtering_stats(filtering_stats, output_dir)

    # Create visualizations
    create_spatial_plots(adata_scored, gene_sets_dict, output_dir)
    create_gene_usage_visualization(cluster_gene_usage, output_dir)
    create_score_heatmaps(pseudobulk_scored, gene_sets_dict, output_dir)

    print("\n" + "=" * 80)
    print("SOLUTION 5 COMPLETED")
    print(f"Results saved to: {output_dir}")
    print("\nKey outputs:")
    print("  - per_cluster_filtering_statistics.csv: Which genes used per cluster")
    print("  - gene_usage_*.png: Visualize which genes scored where")
    print("  - Spatial plots: Scores based on cluster-relevant genes only")
    print("\nThis solution directly addresses the Reelin pathway concern!")
    print("  GCs scored on receptors, CR cells scored on ligands, etc.")
    print("=" * 80)

if __name__ == "__main__":
    main()
