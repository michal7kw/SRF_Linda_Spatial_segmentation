#!/usr/bin/env python3
"""
Solution 2: Expression Fraction Metrics

This approach calculates BOTH the gene set score AND the fraction of genes
in the set that are actually expressed. This helps diagnose whether low scores
are due to true downregulation or simply because few genes are expressed.

Key features:
- Calculates standard gene set scores (same as original)
- Adds fraction_expressed metric per cluster
- Adds number of expressed genes per cluster
- Helps distinguish signal from dilution effects

Use case: Diagnostic tool to understand if dilution is affecting your results.
          Low score + low fraction = dilution problem
          Low score + high fraction = true biological signal
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
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_2_expression_fraction/results'

# Create output directory
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

# Parameters
EXPRESSION_THRESHOLD = 0.5  # Expression level to consider gene "expressed" in pseudobulk

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

def calculate_scores_with_expression_fraction(pseudobulk_adata, gene_sets_dict,
                                               expr_threshold=0.5):
    """
    Calculate gene set scores AND expression fraction metrics.

    For each gene set and cluster, calculates:
    - Standard gene set score
    - Fraction of genes expressed above threshold
    - Number of genes expressed above threshold

    Parameters:
    -----------
    pseudobulk_adata : AnnData
        Pseudo-bulk data
    gene_sets_dict : dict
        Dictionary of gene sets
    expr_threshold : float
        Expression level to consider gene "expressed"

    Returns:
    --------
    tuple
        (pseudobulk_adata, expression_fraction_data)
    """
    print(f"\nCalculating gene set scores with expression fraction metrics...")
    print(f"  Expression threshold: {expr_threshold}")

    expression_fraction_data = []

    for gene_set_name, genes in gene_sets_dict.items():
        print(f"\n  Processing gene set: {gene_set_name}")

        # Filter genes to those in dataset
        genes_in_data = [gene for gene in genes if gene in pseudobulk_adata.var.index]

        if len(genes_in_data) == 0:
            print(f"    WARNING: No genes found in dataset. Skipping.")
            continue

        print(f"    Using {len(genes_in_data)}/{len(genes)} genes")

        # Get indices of genes in data
        gene_indices = [pseudobulk_adata.var.index.get_loc(g) for g in genes_in_data]

        # Calculate standard gene set score
        sc.tl.score_genes(pseudobulk_adata, gene_list=genes_in_data,
                         score_name=f"{gene_set_name}_score", use_raw=False)

        # Calculate expression fraction for each cluster
        for i, cluster in enumerate(pseudobulk_adata.obs.index):
            cluster_expr = pseudobulk_adata.X[i, gene_indices]

            # Count genes expressed above threshold
            n_expressed = (cluster_expr > expr_threshold).sum()
            frac_expressed = n_expressed / len(genes_in_data) if len(genes_in_data) > 0 else 0

            # Store in obs
            pseudobulk_adata.obs.loc[cluster, f"{gene_set_name}_frac_expressed"] = frac_expressed
            pseudobulk_adata.obs.loc[cluster, f"{gene_set_name}_n_expressed"] = n_expressed

            # Store for reporting
            expression_fraction_data.append({
                'Cluster': cluster,
                'Gene_Set': gene_set_name,
                'Score': pseudobulk_adata.obs.loc[cluster, f"{gene_set_name}_score"],
                'Fraction_Expressed': frac_expressed,
                'N_Expressed': n_expressed,
                'N_Total_In_Set': len(genes_in_data),
                'Expression_Threshold': expr_threshold
            })

            print(f"      {cluster}: {n_expressed}/{len(genes_in_data)} genes expressed "
                  f"({frac_expressed:.2%}) | Score: {pseudobulk_adata.obs.loc[cluster, f'{gene_set_name}_score']:.3f}")

    return pseudobulk_adata, expression_fraction_data

def assign_scores_and_fractions_to_cells(adata, pseudobulk_adata, gene_sets_dict,
                                          cluster_key='leiden_new'):
    """Assign both scores and expression fractions to individual cells."""
    print("\nAssigning scores and expression fractions to cells...")

    adata_copy = adata.copy()

    for gene_set_name in gene_sets_dict.keys():
        score_column = f"{gene_set_name}_score"
        frac_column = f"{gene_set_name}_frac_expressed"
        n_expr_column = f"{gene_set_name}_n_expressed"

        if score_column in pseudobulk_adata.obs.columns:
            # Map scores
            score_map = dict(zip(pseudobulk_adata.obs.index,
                                pseudobulk_adata.obs[score_column]))
            adata_copy.obs[f"cluster_{score_column}"] = pd.to_numeric(
                adata_copy.obs[cluster_key].map(score_map), errors='coerce'
            )

            # Map fractions
            frac_map = dict(zip(pseudobulk_adata.obs.index,
                               pseudobulk_adata.obs[frac_column]))
            adata_copy.obs[f"cluster_{frac_column}"] = pd.to_numeric(
                adata_copy.obs[cluster_key].map(frac_map), errors='coerce'
            )

            # Map n_expressed
            n_expr_map = dict(zip(pseudobulk_adata.obs.index,
                                 pseudobulk_adata.obs[n_expr_column]))
            adata_copy.obs[f"cluster_{n_expr_column}"] = pd.to_numeric(
                adata_copy.obs[cluster_key].map(n_expr_map), errors='coerce'
            )

            print(f"  Assigned {gene_set_name} metrics to cells")

    return adata_copy

def create_spatial_plots_with_fraction(adata, gene_sets_dict, output_dir, cluster_key='leiden_new'):
    """Create dual spatial plots showing both score and expression fraction with consistent color scales per time point."""
    print("\nCreating spatial plots with expression fraction...")

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
        frac_column = f"cluster_{gene_set_name}_frac_expressed"

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
            fractions = sample_adata.obs[frac_column]

            # Plot 1: Gene set score with consistent scale
            valid_scores = scores.dropna()
            if len(valid_scores) > 0:
                scatter1 = ax1.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                                      c=scores, cmap='coolwarm', s=5, alpha=0.8,
                                      vmin=vmin_score, vmax=vmax_score)

                cbar1 = plt.colorbar(scatter1, ax=ax1, shrink=0.8)
                cbar1.set_label(f'{gene_set_name} Score', rotation=270, labelpad=20)

                ax1.set_xlabel('X coordinate')
                ax1.set_ylabel('Y coordinate')
                ax1.set_title(f'Gene Set Score\n{gene_set_name}\nScale: {vmin_score:.3f} to {vmax_score:.3f} ({time_point} range)')
                ax1.set_aspect('equal')

            # Plot 2: Expression fraction (always 0-1 scale)
            valid_fractions = fractions.dropna()
            if len(valid_fractions) > 0:
                scatter2 = ax2.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                                      c=fractions, cmap='viridis', s=5, alpha=0.8,
                                      vmin=0, vmax=1)

                cbar2 = plt.colorbar(scatter2, ax=ax2, shrink=0.8)
                cbar2.set_label('Fraction Expressed', rotation=270, labelpad=20)

                ax2.set_xlabel('X coordinate')
                ax2.set_ylabel('Y coordinate')
                ax2.set_title(f'Expression Fraction\n{gene_set_name}\nRange: {valid_fractions.min():.2%} to {valid_fractions.max():.2%}')
                ax2.set_aspect('equal')

            fig.suptitle(f'Sample {sample_id}', fontsize=16, y=0.98)

            plot_filename = f"{clean_name}_with_fraction_sample_{sample_id}.png"
            fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight', dpi=300)
            plt.close(fig)

            print(f"  Saved: {plot_filename}")

def create_diagnostic_plots(expression_fraction_data, output_dir):
    """Create diagnostic plots showing relationship between score and expression fraction."""
    print("\nCreating diagnostic plots...")

    df = pd.DataFrame(expression_fraction_data)

    if df.empty:
        print("  No data to plot")
        return

    # Plot 1: Score vs Fraction Expressed (scatter plot per gene set)
    gene_sets = df['Gene_Set'].unique()

    for gene_set in gene_sets:
        gene_set_df = df[df['Gene_Set'] == gene_set]

        fig, ax = plt.subplots(figsize=(10, 8))

        scatter = ax.scatter(gene_set_df['Fraction_Expressed'],
                           gene_set_df['Score'],
                           c=gene_set_df['N_Expressed'],
                           cmap='viridis',
                           s=100,
                           alpha=0.7)

        # Add cluster labels
        for idx, row in gene_set_df.iterrows():
            ax.annotate(row['Cluster'],
                       (row['Fraction_Expressed'], row['Score']),
                       fontsize=8,
                       alpha=0.7)

        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('N Genes Expressed', rotation=270, labelpad=20)

        ax.set_xlabel('Fraction of Gene Set Expressed')
        ax.set_ylabel('Gene Set Score')
        ax.set_title(f'Score vs Expression Fraction\n{gene_set}\n'
                    f'(Each point is a cluster)')
        ax.grid(True, alpha=0.3)

        clean_name = gene_set.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
        fig.savefig(os.path.join(output_dir, f'diagnostic_scatter_{clean_name}.png'),
                   bbox_inches='tight', dpi=300)
        plt.close()

        print(f"  Saved: diagnostic_scatter_{clean_name}.png")

    # Plot 2: Heatmap of expression fractions
    print("\nCreating expression fraction heatmap...")

    pivot_frac = df.pivot(index='Cluster', columns='Gene_Set', values='Fraction_Expressed')

    fig, ax = plt.subplots(figsize=(12, max(6, len(pivot_frac) * 0.3)))
    sns.heatmap(pivot_frac, cmap='YlGnBu', vmin=0, vmax=1,
               annot=True, fmt='.2f', linewidths=0.5,
               cbar_kws={'label': 'Fraction Expressed'}, ax=ax)

    ax.set_title('Fraction of Gene Set Expressed per Cluster')
    ax.set_xlabel('Gene Sets')
    ax.set_ylabel('Clusters')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'expression_fraction_heatmap.png'),
               bbox_inches='tight', dpi=300)
    plt.close()

    print(f"  Saved: expression_fraction_heatmap.png")

    # Plot 3: Combined score and fraction heatmap
    print("\nCreating combined score/fraction heatmap...")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, max(6, len(pivot_frac) * 0.3)))

    # Score heatmap
    pivot_score = df.pivot(index='Cluster', columns='Gene_Set', values='Score')
    sns.heatmap(pivot_score, cmap='coolwarm', center=0,
               annot=True, fmt='.2f', linewidths=0.5,
               cbar_kws={'label': 'Gene Set Score'}, ax=ax1)
    ax1.set_title('Gene Set Scores')
    ax1.set_xlabel('Gene Sets')
    ax1.set_ylabel('Clusters')

    # Fraction heatmap
    sns.heatmap(pivot_frac, cmap='YlGnBu', vmin=0, vmax=1,
               annot=True, fmt='.2f', linewidths=0.5,
               cbar_kws={'label': 'Fraction Expressed'}, ax=ax2)
    ax2.set_title('Fraction Expressed')
    ax2.set_xlabel('Gene Sets')
    ax2.set_ylabel('Clusters')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'combined_score_fraction_heatmap.png'),
               bbox_inches='tight', dpi=300)
    plt.close()

    print(f"  Saved: combined_score_fraction_heatmap.png")

def save_expression_fraction_data(expression_fraction_data, output_dir):
    """Save expression fraction data to CSV."""
    print("\nSaving expression fraction data...")

    df = pd.DataFrame(expression_fraction_data)
    df.to_csv(os.path.join(output_dir, 'expression_fraction_metrics.csv'), index=False)

    print(f"  Saved: expression_fraction_metrics.csv")

    # Print summary
    print("\nExpression Fraction Summary:")
    summary = df.groupby('Gene_Set').agg({
        'Fraction_Expressed': ['mean', 'min', 'max'],
        'N_Expressed': ['mean', 'min', 'max'],
        'Score': ['mean', 'min', 'max']
    }).round(3)
    print(summary)

def main():
    """Main function."""

    print("=" * 80)
    print("SOLUTION 2: EXPRESSION FRACTION METRICS")
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

    # Calculate scores with expression fractions
    pseudobulk_adata, expression_fraction_data = calculate_scores_with_expression_fraction(
        pseudobulk_adata, gene_sets_dict, expr_threshold=EXPRESSION_THRESHOLD
    )

    # Save expression fraction data
    save_expression_fraction_data(expression_fraction_data, output_dir)

    # Assign to cells
    adata_with_metrics = assign_scores_and_fractions_to_cells(
        adata, pseudobulk_adata, gene_sets_dict, cluster_key=cluster_key
    )

    # Create visualizations
    create_spatial_plots_with_fraction(adata_with_metrics, gene_sets_dict, output_dir,
                                       cluster_key=cluster_key)
    create_diagnostic_plots(expression_fraction_data, output_dir)

    print("\n" + "=" * 80)
    print("SOLUTION 2 COMPLETED")
    print(f"Results saved to: {output_dir}")
    print("\nKey insight: Compare score vs fraction_expressed to identify dilution effects!")
    print("  - Low score + Low fraction = Dilution problem")
    print("  - Low score + High fraction = True biological downregulation")
    print("=" * 80)

if __name__ == "__main__":
    main()
