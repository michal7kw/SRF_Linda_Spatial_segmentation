#!/usr/bin/env python3
"""
Solution 1: Global Expression Filtering with Directional Gene Sets

This approach filters out genes that are not expressed anywhere in the dataset
before calculating gene set scores, while properly handling directional information
(upregulated vs downregulated genes).

Three types of scores:
1. Positive scores (upregulated genes only)
2. Negative scores (downregulated genes only)
3. Aggregate scores (positive - negative)

Key features:
- Removes genes with mean expression < threshold across entire dataset
- All clusters scored using the same filtered gene list per direction
- Good for removing technical noise and missing genes
- Maintains directional information from gene lists
"""

import scanpy as sc
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# Set project root
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Define paths
merged_h5ad_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
gene_sets_path = '/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation/DATA/GENE_LISTS/directional_lists/gene_lists.csv'
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_1_global_filter/results'

os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

# Parameters
GLOBAL_EXPR_THRESHOLD = 0.01

def load_directional_gene_sets(gene_sets_path):
    """Load gene sets with directional information."""
    print(f"Loading directional gene sets from: {gene_sets_path}")

    df = pd.read_csv(gene_sets_path, header=None)
    gene_set_names = df.iloc[0].values
    directions = df.iloc[1].values

    gene_sets_dict = {}
    for col_idx in range(len(gene_set_names)):
        gene_set_name = gene_set_names[col_idx]
        direction = directions[col_idx]

        if pd.isna(gene_set_name) or pd.isna(direction):
            continue

        genes = df.iloc[2:, col_idx].dropna().tolist()
        genes = [g for g in genes if g and str(g).strip()]

        if not genes:
            continue

        if gene_set_name not in gene_sets_dict:
            gene_sets_dict[gene_set_name] = {'up': [], 'down': [], 'all_genes': []}

        if direction == 'up':
            gene_sets_dict[gene_set_name]['up'].extend(genes)
        elif direction == 'down':
            gene_sets_dict[gene_set_name]['down'].extend(genes)

        gene_sets_dict[gene_set_name]['all_genes'].extend(genes)

    print(f"Loaded {len(gene_sets_dict)} gene sets")
    return gene_sets_dict

def filter_genes_by_global_expression(adata, threshold=0.01):
    """Identify genes expressed above threshold across entire dataset."""
    print(f"\nFiltering genes by global expression (threshold={threshold})...")

    if hasattr(adata.X, 'toarray'):
        mean_expr = np.array(adata.X.mean(axis=0)).flatten()
    else:
        mean_expr = adata.X.mean(axis=0)

    expressed_genes = set(adata.var.index[mean_expr > threshold])

    print(f"  Total genes: {adata.n_vars}")
    print(f"  Genes passing threshold: {len(expressed_genes)}")

    return expressed_genes

def create_pseudobulk_by_cluster(adata, cluster_key='leiden_new'):
    """Create pseudo-bulk expression profiles for each cluster."""
    print(f"\nCreating pseudo-bulk profiles...")

    clusters = sorted(adata.obs[cluster_key].unique())
    pseudobulk_data = []

    for cluster in clusters:
        cluster_cells = adata[adata.obs[cluster_key] == cluster]
        mean_expr = cluster_cells.X.mean(axis=0)
        if hasattr(mean_expr, 'A1'):
            mean_expr = mean_expr.A1
        pseudobulk_data.append(mean_expr)

    pseudobulk_adata = sc.AnnData(
        X=np.array(pseudobulk_data),
        obs=pd.DataFrame(index=clusters, data={cluster_key: clusters}),
        var=adata.var.copy()
    )

    print(f"  Created: {pseudobulk_adata.n_obs} clusters x {pseudobulk_adata.n_vars} genes")
    return pseudobulk_adata

def calculate_directional_scores_with_global_filter(pseudobulk_adata, gene_sets_dict, expressed_genes):
    """Calculate directional gene set scores with global filtering."""
    print("\nCalculating directional scores with global filtering...")

    filtering_stats = []

    for gene_set_name, gene_info in gene_sets_dict.items():
        print(f"\n  {gene_set_name}")

        for direction in ['up', 'down']:
            genes = gene_info[direction]
            if not genes:
                continue

            # Filter: in dataset AND globally expressed
            genes_in_data = [g for g in genes if g in pseudobulk_adata.var.index]
            genes_expressed = [g for g in genes_in_data if g in expressed_genes]

            print(f"    {direction}: {len(genes_expressed)}/{len(genes)} genes after filtering")

            if len(genes_expressed) == 0:
                print(f"      WARNING: No genes passed filter")
                continue

            score_name = f"{gene_set_name}_{direction}_score"
            sc.tl.score_genes(pseudobulk_adata, gene_list=genes_expressed,
                            score_name=score_name, use_raw=False)

            filtering_stats.append({
                'Gene_Set': gene_set_name,
                'Direction': direction,
                'Total_Genes': len(genes),
                'Genes_After_Filter': len(genes_expressed),
                'Filtering_Rate': 1 - (len(genes_expressed) / len(genes)) if genes else 0
            })

        # Calculate aggregate score
        up_col = f"{gene_set_name}_up_score"
        down_col = f"{gene_set_name}_down_score"

        if up_col in pseudobulk_adata.obs.columns:
            if down_col not in pseudobulk_adata.obs.columns:
                pseudobulk_adata.obs[down_col] = 0.0
        elif down_col in pseudobulk_adata.obs.columns:
            pseudobulk_adata.obs[up_col] = 0.0

        if up_col in pseudobulk_adata.obs.columns and down_col in pseudobulk_adata.obs.columns:
            pseudobulk_adata.obs[f"{gene_set_name}_aggregate_score"] = \
                pseudobulk_adata.obs[up_col] - pseudobulk_adata.obs[down_col]

    return pseudobulk_adata, filtering_stats

def assign_cluster_scores_to_cells(adata, pseudobulk_adata, gene_sets_dict, cluster_key='leiden_new'):
    """Assign cluster scores to cells."""
    print("\nAssigning scores to cells...")

    adata_copy = adata.copy()

    for gene_set_name in gene_sets_dict.keys():
        for score_type in ['up', 'down', 'aggregate']:
            score_column = f"{gene_set_name}_{score_type}_score"

            if score_column in pseudobulk_adata.obs.columns:
                score_map = dict(zip(pseudobulk_adata.obs.index,
                                    pseudobulk_adata.obs[score_column]))
                adata_copy.obs[f"cluster_{score_column}"] = pd.to_numeric(
                    adata_copy.obs[cluster_key].map(score_map), errors='coerce'
                )

    return adata_copy

def create_spatial_plots(adata, gene_sets_dict, output_dir, score_type='aggregate'):
    """Create spatial plots for specified score type with consistent color scales per time point."""
    print(f"\nCreating spatial plots for {score_type} scores...")

    if 'spatial' not in adata.obsm:
        adata.obsm['spatial'] = np.c_[adata.obs['center_x'], adata.obs['center_y']]

    samples = sorted(adata.obs['sample_id'].unique())

    # Extract time points from sample IDs (e.g., 'e16', 'p0', 'p7', 'p30')
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

    # Colormap selection
    cmaps = {'up': 'Reds', 'down': 'Blues', 'aggregate': 'coolwarm'}
    cmap = cmaps.get(score_type, 'coolwarm')

    for gene_set_name in gene_sets_dict.keys():
        score_column = f"cluster_{gene_set_name}_{score_type}_score"

        if score_column not in adata.obs.columns:
            continue

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

        # Create plots with consistent scales per time point
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
                               c=scores, cmap=cmap, s=5, alpha=0.8,
                               vmin=vmin, vmax=vmax)

            cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
            cbar.set_label(f'{score_type.capitalize()} Score', rotation=270, labelpad=20)

            ax.set_xlabel('X coordinate')
            ax.set_ylabel('Y coordinate')
            ax.set_title(f'{gene_set_name} - {score_type.capitalize()}\n'
                       f'Sample {sample_id} (Global Filter)\n'
                       f'Scale: {vmin:.3f} to {vmax:.3f} ({time_point} range)')
            ax.set_aspect('equal')

            plot_filename = f"{clean_name}_{score_type}_global_filter_sample_{sample_id}.png"
            fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight', dpi=300)
            plt.close(fig)

def create_heatmaps(pseudobulk_adata, gene_sets_dict, output_dir):
    """Create heatmaps for all score types."""
    print("\nCreating heatmaps...")

    for score_type in ['up', 'down', 'aggregate']:
        score_columns = [f"{gs}_{score_type}_score" for gs in gene_sets_dict.keys()
                        if f"{gs}_{score_type}_score" in pseudobulk_adata.obs.columns]

        if len(score_columns) == 0:
            continue

        heatmap_data = pseudobulk_adata.obs[score_columns].copy()
        heatmap_data.columns = [col.replace(f'_{score_type}_score', '') for col in heatmap_data.columns]

        fig, ax = plt.subplots(figsize=(12, max(6, len(heatmap_data) * 0.3)))

        cmap = 'Reds' if score_type == 'up' else ('Blues' if score_type == 'down' else 'coolwarm')
        center_val = 0 if score_type == 'aggregate' else None

        sns.heatmap(heatmap_data, cmap=cmap, center=center_val,
                   annot=True, fmt='.3f', linewidths=0.5,
                   cbar_kws={'label': f'{score_type.capitalize()} Score'}, ax=ax)

        ax.set_title(f'Cluster {score_type.capitalize()} Scores (Global Filter)')
        ax.set_xlabel('Gene Sets')
        ax.set_ylabel('Clusters')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'heatmap_{score_type}_global_filter.png'),
                   bbox_inches='tight', dpi=300)
        plt.close()

def main():
    print("=" * 80)
    print("SOLUTION 1: GLOBAL FILTERING WITH DIRECTIONAL GENE SETS")
    print("=" * 80)

    # Load data
    adata = sc.read_h5ad(merged_h5ad_path)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # Load directional gene sets
    gene_sets_dict = load_directional_gene_sets(gene_sets_path)

    # Filter genes globally
    expressed_genes = filter_genes_by_global_expression(adata, GLOBAL_EXPR_THRESHOLD)

    # Create pseudobulk
    pseudobulk_adata = create_pseudobulk_by_cluster(adata)

    # Calculate scores
    pseudobulk_adata, filtering_stats = calculate_directional_scores_with_global_filter(
        pseudobulk_adata, gene_sets_dict, expressed_genes
    )

    # Save stats
    pd.DataFrame(filtering_stats).to_csv(
        os.path.join(output_dir, 'filtering_statistics.csv'), index=False
    )

    # Assign to cells
    adata_with_scores = assign_cluster_scores_to_cells(adata, pseudobulk_adata, gene_sets_dict)

    # Create visualizations
    for score_type in ['up', 'down', 'aggregate']:
        create_spatial_plots(adata_with_scores, gene_sets_dict, output_dir, score_type)

    create_heatmaps(pseudobulk_adata, gene_sets_dict, output_dir)

    print("\n" + "=" * 80)
    print("SOLUTION 1 COMPLETED")
    print(f"Results saved to: {output_dir}")
    print("=" * 80)

if __name__ == "__main__":
    main()
