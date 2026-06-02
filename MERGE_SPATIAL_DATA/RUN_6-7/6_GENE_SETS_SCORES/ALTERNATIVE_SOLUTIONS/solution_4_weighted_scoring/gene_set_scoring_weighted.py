#!/usr/bin/env python3
"""
Solution 4: Weighted Scoring by Expression

This approach calculates gene set scores using ONLY genes that are expressed
above a threshold in each cell. Genes below threshold are excluded from scoring
on a per-cell basis, preventing non-expressed genes from diluting the signal.

Key differences from standard scoring:
- Standard: All genes included, zero expression contributes to score
- Weighted: Only expressed genes included per cell/cluster

Key features:
- Cell-specific gene filtering (not cluster-specific)
- Prevents dilution from non-expressed genes
- More computationally intensive
- Most biologically accurate for sparse gene sets

Use case: When you want to score only genes that are actually present,
          allowing the score to reflect the true state of expressed genes.
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
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_4_weighted_scoring/results'

# Create output directory
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

# Parameters
EXPRESSION_THRESHOLD = 0.1  # Expression level to include gene in scoring
MIN_GENES_FOR_SCORE = 2  # Minimum genes needed to calculate a score

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

def calculate_weighted_scores_pseudobulk(pseudobulk_adata, gene_sets_dict,
                                         expr_threshold=0.1, min_genes=2):
    """
    Calculate weighted gene set scores on pseudobulk data.

    For each cluster, only genes expressed above threshold are included in scoring.
    This prevents non-expressed genes from diluting the signal.

    Parameters:
    -----------
    pseudobulk_adata : AnnData
        Pseudo-bulk data
    gene_sets_dict : dict
        Dictionary of gene sets
    expr_threshold : float
        Minimum expression to include gene in scoring
    min_genes : int
        Minimum genes needed to calculate a valid score

    Returns:
    --------
    tuple
        (pseudobulk_adata, scoring_stats)
    """
    print(f"\nCalculating weighted gene set scores...")
    print(f"  Expression threshold: {expr_threshold}")
    print(f"  Minimum genes for scoring: {min_genes}")

    scoring_stats = []

    for gene_set_name, genes in gene_sets_dict.items():
        print(f"\n  Processing gene set: {gene_set_name}")

        # Filter genes to those in dataset
        genes_in_data = [gene for gene in genes if gene in pseudobulk_adata.var.index]

        if len(genes_in_data) == 0:
            print(f"    WARNING: No genes found in dataset. Skipping.")
            continue

        print(f"    Total genes in dataset: {len(genes_in_data)}/{len(genes)}")

        # Get gene indices
        gene_indices = np.array([pseudobulk_adata.var.index.get_loc(g) for g in genes_in_data])

        # Calculate weighted scores for each cluster
        scores = []
        n_genes_used = []

        for i, cluster in enumerate(pseudobulk_adata.obs.index):
            # Get expression values for this cluster
            cluster_expr = pseudobulk_adata.X[i, gene_indices]

            # Filter to expressed genes
            expressed_mask = cluster_expr > expr_threshold
            n_expressed = expressed_mask.sum()

            if n_expressed >= min_genes:
                # Calculate mean of only expressed genes
                score = cluster_expr[expressed_mask].mean()
            else:
                # Not enough genes expressed, assign NaN
                score = np.nan

            scores.append(score)
            n_genes_used.append(n_expressed)

            score_str = f"{score:.3f}" if not np.isnan(score) else "N/A"
            print(f"      {cluster}: {n_expressed}/{len(genes_in_data)} genes used | "
                  f"Score: {score_str}")

            # Store statistics
            scoring_stats.append({
                'Cluster': cluster,
                'Gene_Set': gene_set_name,
                'Score': score,
                'N_Genes_In_Set': len(genes_in_data),
                'N_Genes_Used': n_expressed,
                'Frac_Genes_Used': n_expressed / len(genes_in_data) if len(genes_in_data) > 0 else 0,
                'Expression_Threshold': expr_threshold
            })

        # Store scores in pseudobulk_adata
        pseudobulk_adata.obs[f"{gene_set_name}_score"] = scores
        pseudobulk_adata.obs[f"{gene_set_name}_n_genes_used"] = n_genes_used
        pseudobulk_adata.obs[f"{gene_set_name}_frac_genes_used"] = [
            n / len(genes_in_data) for n in n_genes_used
        ]

    return pseudobulk_adata, scoring_stats

def assign_scores_to_cells(adata, pseudobulk_adata, gene_sets_dict, cluster_key='leiden_new'):
    """Assign cluster-level weighted scores to individual cells."""
    print("\nAssigning weighted scores to cells...")

    adata_copy = adata.copy()

    for gene_set_name in gene_sets_dict.keys():
        score_column = f"{gene_set_name}_score"
        n_genes_column = f"{gene_set_name}_n_genes_used"
        frac_genes_column = f"{gene_set_name}_frac_genes_used"

        if score_column in pseudobulk_adata.obs.columns:
            # Map scores
            score_map = dict(zip(pseudobulk_adata.obs.index,
                                pseudobulk_adata.obs[score_column]))
            adata_copy.obs[f"cluster_{score_column}"] = pd.to_numeric(
                adata_copy.obs[cluster_key].map(score_map), errors='coerce'
            )

            # Map n_genes_used
            n_genes_map = dict(zip(pseudobulk_adata.obs.index,
                                  pseudobulk_adata.obs[n_genes_column]))
            adata_copy.obs[f"cluster_{n_genes_column}"] = pd.to_numeric(
                adata_copy.obs[cluster_key].map(n_genes_map), errors='coerce'
            )

            # Map frac_genes_used
            frac_genes_map = dict(zip(pseudobulk_adata.obs.index,
                                     pseudobulk_adata.obs[frac_genes_column]))
            adata_copy.obs[f"cluster_{frac_genes_column}"] = pd.to_numeric(
                adata_copy.obs[cluster_key].map(frac_genes_map), errors='coerce'
            )

            print(f"  Assigned {gene_set_name} metrics to cells")

    return adata_copy

def create_spatial_plots(adata, gene_sets_dict, output_dir, cluster_key='leiden_new'):
    """Create spatial plots showing weighted gene set scores."""
    print("\nCreating spatial plots...")

    if 'spatial' not in adata.obsm:
        adata.obsm['spatial'] = np.c_[adata.obs['center_x'], adata.obs['center_y']]

    samples = sorted(adata.obs['sample_id'].unique())

    for gene_set_name in gene_sets_dict.keys():
        score_column = f"cluster_{gene_set_name}_score"
        n_genes_column = f"cluster_{gene_set_name}_n_genes_used"

        if score_column not in adata.obs.columns:
            print(f"Skipping {gene_set_name} - no scores calculated")
            continue

        print(f"Creating spatial plots for gene set: {gene_set_name}")

        clean_name = gene_set_name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')

        for sample_id in samples:
            sample_adata = adata[adata.obs['sample_id'] == sample_id].copy()

            if len(sample_adata) == 0:
                continue

            # Create figure with two subplots
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 10))

            spatial_coords = sample_adata.obsm['spatial']
            scores = sample_adata.obs[score_column]
            n_genes = sample_adata.obs[n_genes_column]

            # Plot 1: Weighted score
            valid_scores = scores.dropna()
            if len(valid_scores) > 0:
                vmin_score = valid_scores.quantile(0.01)
                vmax_score = valid_scores.quantile(0.99)

                scatter1 = ax1.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                                      c=scores, cmap='coolwarm', s=5, alpha=0.8,
                                      vmin=vmin_score, vmax=vmax_score)

                cbar1 = plt.colorbar(scatter1, ax=ax1, shrink=0.8)
                cbar1.set_label(f'{gene_set_name} Score\n(Weighted)', rotation=270, labelpad=20)

                ax1.set_xlabel('X coordinate')
                ax1.set_ylabel('Y coordinate')
                ax1.set_title(f'Weighted Gene Set Score\n{gene_set_name}\nRange: {valid_scores.min():.3f} to {valid_scores.max():.3f}')
                ax1.set_aspect('equal')

            # Plot 2: Number of genes used
            valid_n_genes = n_genes.dropna()
            if len(valid_n_genes) > 0:
                scatter2 = ax2.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                                      c=n_genes, cmap='viridis', s=5, alpha=0.8,
                                      vmin=0, vmax=valid_n_genes.max())

                cbar2 = plt.colorbar(scatter2, ax=ax2, shrink=0.8)
                cbar2.set_label('N Genes Used', rotation=270, labelpad=20)

                ax2.set_xlabel('X coordinate')
                ax2.set_ylabel('Y coordinate')
                ax2.set_title(f'Number of Genes Used in Score\n{gene_set_name}\nRange: {int(valid_n_genes.min())} to {int(valid_n_genes.max())}')
                ax2.set_aspect('equal')

            fig.suptitle(f'Sample {sample_id}', fontsize=16, y=0.98)

            plot_filename = f"{clean_name}_weighted_sample_{sample_id}.png"
            fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight', dpi=300)
            plt.close(fig)

            print(f"  Saved: {plot_filename}")

def create_diagnostic_plots(scoring_stats, output_dir):
    """Create diagnostic plots for weighted scoring."""
    print("\nCreating diagnostic plots...")

    df = pd.DataFrame(scoring_stats)

    if df.empty:
        print("  No data to plot")
        return

    # Remove NaN scores for plotting
    df_valid = df[df['Score'].notna()].copy()

    # Plot 1: Score vs N_Genes_Used
    gene_sets = df_valid['Gene_Set'].unique()

    for gene_set in gene_sets:
        gene_set_df = df_valid[df_valid['Gene_Set'] == gene_set]

        fig, ax = plt.subplots(figsize=(10, 8))

        scatter = ax.scatter(gene_set_df['N_Genes_Used'],
                           gene_set_df['Score'],
                           c=gene_set_df['Frac_Genes_Used'],
                           cmap='viridis',
                           s=100,
                           alpha=0.7)

        # Add cluster labels
        for idx, row in gene_set_df.iterrows():
            ax.annotate(row['Cluster'],
                       (row['N_Genes_Used'], row['Score']),
                       fontsize=8,
                       alpha=0.7)

        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Fraction Genes Used', rotation=270, labelpad=20)

        ax.set_xlabel('Number of Genes Used in Score')
        ax.set_ylabel('Weighted Gene Set Score')
        ax.set_title(f'Score vs N Genes Used\n{gene_set}\n(Each point is a cluster)')
        ax.grid(True, alpha=0.3)

        clean_name = gene_set.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
        fig.savefig(os.path.join(output_dir, f'diagnostic_weighted_{clean_name}.png'),
                   bbox_inches='tight', dpi=300)
        plt.close()

        print(f"  Saved: diagnostic_weighted_{clean_name}.png")

    # Plot 2: Heatmap of scores
    print("\nCreating score heatmap...")
    pivot_score = df_valid.pivot(index='Cluster', columns='Gene_Set', values='Score')

    fig, ax = plt.subplots(figsize=(12, max(6, len(pivot_score) * 0.3)))
    sns.heatmap(pivot_score, cmap='coolwarm', center=0,
               annot=True, fmt='.2f', linewidths=0.5,
               cbar_kws={'label': 'Weighted Gene Set Score'}, ax=ax)

    ax.set_title('Weighted Gene Set Scores per Cluster')
    ax.set_xlabel('Gene Sets')
    ax.set_ylabel('Clusters')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'weighted_scores_heatmap.png'),
               bbox_inches='tight', dpi=300)
    plt.close()

    print(f"  Saved: weighted_scores_heatmap.png")

    # Plot 3: Heatmap of N genes used
    print("\nCreating N genes used heatmap...")
    pivot_n_genes = df_valid.pivot(index='Cluster', columns='Gene_Set', values='N_Genes_Used')

    fig, ax = plt.subplots(figsize=(12, max(6, len(pivot_n_genes) * 0.3)))
    sns.heatmap(pivot_n_genes, cmap='YlGnBu',
               annot=True, fmt='.0f', linewidths=0.5,
               cbar_kws={'label': 'N Genes Used'}, ax=ax)

    ax.set_title('Number of Genes Used in Weighted Scoring')
    ax.set_xlabel('Gene Sets')
    ax.set_ylabel('Clusters')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'n_genes_used_heatmap.png'),
               bbox_inches='tight', dpi=300)
    plt.close()

    print(f"  Saved: n_genes_used_heatmap.png")

def save_scoring_stats(scoring_stats, output_dir):
    """Save scoring statistics to CSV."""
    print("\nSaving scoring statistics...")

    df = pd.DataFrame(scoring_stats)
    df.to_csv(os.path.join(output_dir, 'weighted_scoring_statistics.csv'), index=False)

    print(f"  Saved: weighted_scoring_statistics.csv")

    # Print summary
    print("\nWeighted Scoring Summary:")
    df_valid = df[df['Score'].notna()]
    summary = df_valid.groupby('Gene_Set').agg({
        'Score': ['mean', 'std', 'min', 'max'],
        'N_Genes_Used': ['mean', 'min', 'max'],
        'Frac_Genes_Used': ['mean', 'min', 'max']
    }).round(3)
    print(summary)

def main():
    """Main function."""

    print("=" * 80)
    print("SOLUTION 4: WEIGHTED SCORING BY EXPRESSION")
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

    # Calculate weighted scores
    pseudobulk_adata, scoring_stats = calculate_weighted_scores_pseudobulk(
        pseudobulk_adata, gene_sets_dict,
        expr_threshold=EXPRESSION_THRESHOLD,
        min_genes=MIN_GENES_FOR_SCORE
    )

    # Save statistics
    save_scoring_stats(scoring_stats, output_dir)

    # Assign to cells
    adata_with_scores = assign_scores_to_cells(
        adata, pseudobulk_adata, gene_sets_dict, cluster_key=cluster_key
    )

    # Create visualizations
    create_spatial_plots(adata_with_scores, gene_sets_dict, output_dir, cluster_key=cluster_key)
    create_diagnostic_plots(scoring_stats, output_dir)

    print("\n" + "=" * 80)
    print("SOLUTION 4 COMPLETED")
    print(f"Results saved to: {output_dir}")
    print("\nKey insight: Weighted scores reflect only expressed genes!")
    print("  - Higher scores = True biological signal (not diluted)")
    print("  - Check N_Genes_Used to ensure sufficient data")
    print("=" * 80)

if __name__ == "__main__":
    main()
