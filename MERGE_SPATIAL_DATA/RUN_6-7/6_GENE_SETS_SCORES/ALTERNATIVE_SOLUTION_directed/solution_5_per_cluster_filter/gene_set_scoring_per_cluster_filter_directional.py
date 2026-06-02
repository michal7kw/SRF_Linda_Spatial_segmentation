#!/usr/bin/env python3
"""
Solution 5: Per-Cluster Gene Filtering with Directional Gene Sets

THE SOLUTION for the dilution problem with directional gene sets!

For each cluster, genes are filtered based on whether they are expressed in THAT
specific cluster before scoring, separately for up and down genes.

Example: Migration pathway
- Cluster A: Scored on 5 upregulated + 2 downregulated genes (only ones expressed)
- Cluster B: Scored on 3 upregulated + 4 downregulated genes
- Aggregate = (up score) - (down score) for each cluster

Key features:
- Cluster-specific gene filtering per direction
- Each cluster scored on biologically relevant gene subsets
- Directly addresses cell-type-specific pathway components
- Maintains directional information throughout
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

project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

merged_h5ad_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
gene_sets_path = '/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation/DATA/GENE_LISTS/directional_lists/gene_lists.csv'
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_5_per_cluster_filter/results'

os.makedirs(output_dir, exist_ok=True)

CLUSTER_EXPR_THRESHOLD = 0.5
MIN_GENES_FOR_SCORE = 2

def load_directional_gene_sets(gene_sets_path):
    df = pd.read_csv(gene_sets_path, header=None)
    gene_set_names, directions = df.iloc[0].values, df.iloc[1].values
    gene_sets_dict = {}

    for col_idx in range(len(gene_set_names)):
        if pd.isna(gene_set_names[col_idx]) or pd.isna(directions[col_idx]):
            continue
        gene_set_name, direction = gene_set_names[col_idx], directions[col_idx]
        genes = [g for g in df.iloc[2:, col_idx].dropna().tolist() if g and str(g).strip()]

        if gene_set_name not in gene_sets_dict:
            gene_sets_dict[gene_set_name] = {'up': [], 'down': []}

        if direction in ['up', 'down']:
            gene_sets_dict[gene_set_name][direction].extend(genes)

    return gene_sets_dict

def create_pseudobulk(adata, cluster_key='leiden_new'):
    clusters = sorted(adata.obs[cluster_key].unique())
    pseudobulk_data = []
    for cluster in clusters:
        mean_expr = adata[adata.obs[cluster_key] == cluster].X.mean(axis=0)
        pseudobulk_data.append(mean_expr.A1 if hasattr(mean_expr, 'A1') else mean_expr)

    return sc.AnnData(
        X=np.array(pseudobulk_data),
        obs=pd.DataFrame(index=clusters, data={cluster_key: clusters}),
        var=adata.var.copy()
    )

def filter_genes_per_cluster(pseudobulk_adata, genes, expr_threshold=0.5):
    """For each cluster, identify which genes are expressed."""
    genes_in_data = [g for g in genes if g in pseudobulk_adata.var.index]
    if not genes_in_data:
        return {}

    gene_indices = [pseudobulk_adata.var.index.get_loc(g) for g in genes_in_data]
    cluster_gene_map = {}

    for i, cluster in enumerate(pseudobulk_adata.obs.index):
        cluster_expr = pseudobulk_adata.X[i, gene_indices]
        expressed_mask = cluster_expr > expr_threshold
        expressed_genes = [genes_in_data[j] for j, is_expr in enumerate(expressed_mask) if is_expr]
        cluster_gene_map[cluster] = expressed_genes

    return cluster_gene_map

def calculate_per_cluster_directional_scores(pseudobulk_adata, adata, gene_sets_dict,
                                             cluster_key='leiden_new',
                                             expr_threshold=0.5,
                                             min_genes=2):
    """Calculate directional gene set scores with per-cluster filtering."""
    print("\nCalculating per-cluster filtered directional scores...")

    filtering_stats = []
    adata_scored = adata.copy()

    for gene_set_name, gene_info in gene_sets_dict.items():
        print(f"\n  {gene_set_name}")

        cluster_scores = {direction: {} for direction in ['up', 'down']}
        cluster_n_genes = {direction: {} for direction in ['up', 'down']}

        for direction in ['up', 'down']:
            genes = gene_info[direction]
            if not genes:
                continue

            # Get cluster-specific gene lists
            cluster_gene_map = filter_genes_per_cluster(pseudobulk_adata, genes, expr_threshold)

            # Score each cluster separately
            for cluster, cluster_genes in cluster_gene_map.items():
                n_genes = len(cluster_genes)

                print(f"    {cluster} ({direction}): {n_genes} genes")

                if n_genes < min_genes:
                    cluster_scores[direction][cluster] = np.nan
                    cluster_n_genes[direction][cluster] = n_genes
                    continue

                # Get cells from this cluster
                cluster_cells = adata_scored[adata_scored.obs[cluster_key] == cluster].copy()

                if len(cluster_cells) == 0:
                    cluster_scores[direction][cluster] = np.nan
                    cluster_n_genes[direction][cluster] = 0
                    continue

                # Score this cluster with its specific gene list
                score_col = f"{gene_set_name}_{direction}_score"
                sc.tl.score_genes(cluster_cells, gene_list=cluster_genes,
                                score_name=score_col, use_raw=False)

                # Store mean score
                cluster_scores[direction][cluster] = cluster_cells.obs[score_col].mean()
                cluster_n_genes[direction][cluster] = n_genes

                # Assign scores to cells
                cell_indices = adata_scored.obs[cluster_key] == cluster
                adata_scored.obs.loc[cell_indices, score_col] = cluster_cells.obs[score_col].values

                # Store stats
                filtering_stats.append({
                    'Gene_Set': gene_set_name,
                    'Cluster': cluster,
                    'Direction': direction,
                    'N_Genes_Total': len([g for g in genes if g in pseudobulk_adata.var.index]),
                    'N_Genes_Used': n_genes,
                    'Frac_Genes_Used': n_genes / len([g for g in genes if g in pseudobulk_adata.var.index]) if genes else 0,
                    'Score': cluster_scores[direction][cluster]
                })

        # Calculate aggregate scores per cluster
        for cluster in pseudobulk_adata.obs.index:
            up_score = cluster_scores['up'].get(cluster, 0.0)
            down_score = cluster_scores['down'].get(cluster, 0.0)

            if not np.isnan(up_score) or not np.isnan(down_score):
                up_score = 0.0 if np.isnan(up_score) else up_score
                down_score = 0.0 if np.isnan(down_score) else down_score
                aggregate_score = up_score - down_score

                # Assign to pseudobulk
                pseudobulk_adata.obs.loc[cluster, f"{gene_set_name}_aggregate_score"] = aggregate_score

                # Assign to cells
                cell_indices = adata_scored.obs[cluster_key] == cluster
                adata_scored.obs.loc[cell_indices, f"{gene_set_name}_aggregate_score"] = aggregate_score

        # Add to pseudobulk
        for direction in ['up', 'down']:
            score_col = f"{gene_set_name}_{direction}_score"
            n_genes_col = f"{gene_set_name}_{direction}_n_genes"

            pseudobulk_adata.obs[score_col] = [cluster_scores[direction].get(c, np.nan)
                                                for c in pseudobulk_adata.obs.index]
            pseudobulk_adata.obs[n_genes_col] = [cluster_n_genes[direction].get(c, 0)
                                                  for c in pseudobulk_adata.obs.index]

    return adata_scored, pseudobulk_adata, filtering_stats

def create_spatial_plots(adata, gene_sets_dict, output_dir):
    """Create spatial plots for all score types with consistent color scales per time point."""
    print("\nCreating spatial plots...")

    if 'spatial' not in adata.obsm:
        adata.obsm['spatial'] = np.c_[adata.obs['center_x'], adata.obs['center_y']]

    samples = sorted(adata.obs['sample_id'].unique())
    cmaps = {'up': 'Reds', 'down': 'Blues', 'aggregate': 'coolwarm'}

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
        clean_name = gene_set_name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')

        for score_type in ['up', 'down', 'aggregate']:
            score_column = f"{gene_set_name}_{score_type}_score"

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
                                   c=scores, cmap=cmaps[score_type], s=5, alpha=0.8,
                                   vmin=vmin, vmax=vmax)

                cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
                cbar.set_label(f'{score_type.capitalize()} Score', rotation=270, labelpad=20)

                ax.set_xlabel('X coordinate')
                ax.set_ylabel('Y coordinate')
                ax.set_title(f'{gene_set_name} - {score_type.capitalize()}\n'
                           f'Sample {sample_id} (Per-Cluster Filter)\n'
                           f'Scale: {vmin:.3f} to {vmax:.3f} ({time_point} range)')
                ax.set_aspect('equal')

                plot_filename = f"{clean_name}_{score_type}_per_cluster_sample_{sample_id}.png"
                fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight', dpi=300)
                plt.close(fig)

def create_heatmaps(pseudobulk_adata, gene_sets_dict, output_dir):
    """Create heatmaps for all score types."""
    print("\nCreating heatmaps...")

    for score_type in ['up', 'down', 'aggregate']:
        score_columns = [f"{gs}_{score_type}_score" for gs in gene_sets_dict.keys()
                        if f"{gs}_{score_type}_score" in pseudobulk_adata.obs.columns]

        if not score_columns:
            continue

        heatmap_data = pseudobulk_adata.obs[score_columns].copy()
        heatmap_data.columns = [col.replace(f'_{score_type}_score', '') for col in heatmap_data.columns]

        fig, ax = plt.subplots(figsize=(12, max(6, len(heatmap_data) * 0.3)))

        cmap = 'Reds' if score_type == 'up' else ('Blues' if score_type == 'down' else 'coolwarm')
        center_val = 0 if score_type == 'aggregate' else None

        sns.heatmap(heatmap_data, cmap=cmap, center=center_val,
                   annot=True, fmt='.3f', linewidths=0.5,
                   cbar_kws={'label': f'{score_type.capitalize()} Score'}, ax=ax)

        ax.set_title(f'Cluster {score_type.capitalize()} Scores (Per-Cluster Filter)')
        ax.set_xlabel('Gene Sets')
        ax.set_ylabel('Clusters')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'heatmap_{score_type}_per_cluster.png'),
                   bbox_inches='tight', dpi=300)
        plt.close()

def main():
    print("=" * 80)
    print("SOLUTION 5: PER-CLUSTER FILTERING WITH DIRECTIONAL GENE SETS")
    print("=" * 80)

    adata = sc.read_h5ad(merged_h5ad_path)
    gene_sets_dict = load_directional_gene_sets(gene_sets_path)
    pseudobulk_adata = create_pseudobulk(adata)

    adata_scored, pseudobulk_adata, filtering_stats = calculate_per_cluster_directional_scores(
        pseudobulk_adata, adata, gene_sets_dict,
        expr_threshold=CLUSTER_EXPR_THRESHOLD,
        min_genes=MIN_GENES_FOR_SCORE
    )

    pd.DataFrame(filtering_stats).to_csv(
        os.path.join(output_dir, 'per_cluster_filtering_stats.csv'), index=False
    )

    create_spatial_plots(adata_scored, gene_sets_dict, output_dir)
    create_heatmaps(pseudobulk_adata, gene_sets_dict, output_dir)

    print("\n" + "=" * 80)
    print("SOLUTION 5 COMPLETED - THE BEST SOLUTION FOR DILUTION PROBLEMS!")
    print(f"Results: {output_dir}")
    print("=" * 80)

if __name__ == "__main__":
    main()
