#!/usr/bin/env python3
"""
Solution 2: Expression Fraction Metrics with Directional Gene Sets

Calculates gene set scores AND fraction of genes expressed, for directional gene sets.
Helps diagnose whether low scores are due to true downregulation or dilution.

Three score types (up, down, aggregate) with expression fraction metrics for each.
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
gene_sets_path = 'DATA/GENE_LISTS/directional_lists/gene_lists.csv'
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_2_expression_fraction/results'

os.makedirs(output_dir, exist_ok=True)
EXPRESSION_THRESHOLD = 0.5

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

def calculate_scores_with_fraction(pseudobulk_adata, gene_sets_dict, expr_threshold=0.5):
    """Calculate directional scores with expression fraction metrics."""
    print("\nCalculating scores with expression fractions...")

    fraction_data = []

    for gene_set_name, gene_info in gene_sets_dict.items():
        print(f"\n  {gene_set_name}")

        for direction in ['up', 'down']:
            genes = gene_info[direction]
            if not genes:
                continue

            genes_in_data = [g for g in genes if g in pseudobulk_adata.var.index]
            if not genes_in_data:
                continue

            # Calculate scores
            sc.tl.score_genes(pseudobulk_adata, gene_list=genes_in_data,
                            score_name=f"{gene_set_name}_{direction}_score", use_raw=False)

            # Calculate expression fractions per cluster
            gene_indices = [pseudobulk_adata.var.index.get_loc(g) for g in genes_in_data]

            for i, cluster in enumerate(pseudobulk_adata.obs.index):
                cluster_expr = pseudobulk_adata.X[i, gene_indices]
                n_expressed = (cluster_expr > expr_threshold).sum()
                frac_expressed = n_expressed / len(genes_in_data)

                pseudobulk_adata.obs.loc[cluster, f"{gene_set_name}_{direction}_frac_expressed"] = frac_expressed
                pseudobulk_adata.obs.loc[cluster, f"{gene_set_name}_{direction}_n_expressed"] = n_expressed

                fraction_data.append({
                    'Cluster': cluster,
                    'Gene_Set': gene_set_name,
                    'Direction': direction,
                    'Score': pseudobulk_adata.obs.loc[cluster, f"{gene_set_name}_{direction}_score"],
                    'Fraction_Expressed': frac_expressed,
                    'N_Expressed': n_expressed,
                    'N_Total': len(genes_in_data)
                })

        # Calculate aggregate scores and fractions
        up_score_col = f"{gene_set_name}_up_score"
        down_score_col = f"{gene_set_name}_down_score"
        up_frac_col = f"{gene_set_name}_up_frac_expressed"
        down_frac_col = f"{gene_set_name}_down_frac_expressed"
        up_n_col = f"{gene_set_name}_up_n_expressed"
        down_n_col = f"{gene_set_name}_down_n_expressed"

        if up_score_col in pseudobulk_adata.obs.columns:
            if down_score_col not in pseudobulk_adata.obs.columns:
                pseudobulk_adata.obs[down_score_col] = 0.0
        elif down_score_col in pseudobulk_adata.obs.columns:
            pseudobulk_adata.obs[up_score_col] = 0.0

        if up_score_col in pseudobulk_adata.obs.columns and down_score_col in pseudobulk_adata.obs.columns:
            # Aggregate score
            pseudobulk_adata.obs[f"{gene_set_name}_aggregate_score"] = \
                pseudobulk_adata.obs[up_score_col] - pseudobulk_adata.obs[down_score_col]

            # Aggregate fraction (combined fraction of all genes - up and down)
            for cluster in pseudobulk_adata.obs.index:
                n_up = pseudobulk_adata.obs.loc[cluster, up_n_col] if up_n_col in pseudobulk_adata.obs.columns else 0
                n_down = pseudobulk_adata.obs.loc[cluster, down_n_col] if down_n_col in pseudobulk_adata.obs.columns else 0

                n_up_total = len(gene_info['up']) if gene_info['up'] else 0
                n_down_total = len(gene_info['down']) if gene_info['down'] else 0
                n_total_genes = n_up_total + n_down_total

                if n_total_genes > 0:
                    frac_aggregate = (n_up + n_down) / n_total_genes
                else:
                    frac_aggregate = 0.0

                pseudobulk_adata.obs.loc[cluster, f"{gene_set_name}_aggregate_frac_expressed"] = frac_aggregate
                pseudobulk_adata.obs.loc[cluster, f"{gene_set_name}_aggregate_n_expressed"] = n_up + n_down

    return pseudobulk_adata, fraction_data

def assign_metrics_to_cells(adata, pseudobulk_adata, gene_sets_dict, cluster_key='leiden_new'):
    adata_copy = adata.copy()

    for gene_set_name in gene_sets_dict.keys():
        for direction in ['up', 'down', 'aggregate']:
            score_col = f"{gene_set_name}_{direction}_score"
            frac_col = f"{gene_set_name}_{direction}_frac_expressed"

            if score_col in pseudobulk_adata.obs.columns:
                score_map = dict(zip(pseudobulk_adata.obs.index, pseudobulk_adata.obs[score_col]))
                adata_copy.obs[f"cluster_{score_col}"] = pd.to_numeric(
                    adata_copy.obs[cluster_key].map(score_map), errors='coerce'
                )

                if frac_col in pseudobulk_adata.obs.columns:
                    frac_map = dict(zip(pseudobulk_adata.obs.index, pseudobulk_adata.obs[frac_col]))
                    adata_copy.obs[f"cluster_{frac_col}"] = pd.to_numeric(
                        adata_copy.obs[cluster_key].map(frac_map), errors='coerce'
                    )

    return adata_copy

def create_spatial_plots_with_fraction(adata, gene_sets_dict, output_dir):
    """Create dual spatial plots showing both score and expression fraction."""
    print("\nCreating spatial plots with expression fraction...")

    if 'spatial' not in adata.obsm:
        adata.obsm['spatial'] = np.c_[adata.obs['center_x'], adata.obs['center_y']]

    samples = sorted(adata.obs['sample_id'].unique())
    cmaps = {'up': 'Reds', 'down': 'Blues', 'aggregate': 'coolwarm'}

    for gene_set_name in gene_sets_dict.keys():
        clean_name = gene_set_name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')

        for direction in ['up', 'down', 'aggregate']:
            score_column = f"cluster_{gene_set_name}_{direction}_score"
            frac_column = f"cluster_{gene_set_name}_{direction}_frac_expressed"

            if score_column not in adata.obs.columns:
                continue

            print(f"Creating spatial plots for: {gene_set_name} ({direction})")

            for sample_id in samples:
                sample_adata = adata[adata.obs['sample_id'] == sample_id].copy()

                if len(sample_adata) == 0:
                    continue

                # Create figure with two subplots
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 10))

                spatial_coords = sample_adata.obsm['spatial']
                scores = sample_adata.obs[score_column]
                fractions = sample_adata.obs[frac_column] if frac_column in sample_adata.obs.columns else None

                # Plot 1: Gene set score
                valid_scores = scores.dropna()
                if len(valid_scores) > 0:
                    vmin_score = valid_scores.quantile(0.01)
                    vmax_score = valid_scores.quantile(0.99)

                    scatter1 = ax1.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                                          c=scores, cmap=cmaps[direction], s=5, alpha=0.8,
                                          vmin=vmin_score, vmax=vmax_score)

                    cbar1 = plt.colorbar(scatter1, ax=ax1, shrink=0.8)
                    cbar1.set_label(f'{direction.capitalize()} Score', rotation=270, labelpad=20)

                    ax1.set_xlabel('X coordinate')
                    ax1.set_ylabel('Y coordinate')
                    ax1.set_title(f'{direction.capitalize()} Score\n{gene_set_name}\n'
                                 f'Range: {valid_scores.min():.3f} to {valid_scores.max():.3f}')
                    ax1.set_aspect('equal')

                # Plot 2: Expression fraction
                if fractions is not None:
                    valid_fractions = fractions.dropna()
                    if len(valid_fractions) > 0:
                        scatter2 = ax2.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                                              c=fractions, cmap='viridis', s=5, alpha=0.8,
                                              vmin=0, vmax=1)

                        cbar2 = plt.colorbar(scatter2, ax=ax2, shrink=0.8)
                        cbar2.set_label('Fraction Expressed', rotation=270, labelpad=20)

                        ax2.set_xlabel('X coordinate')
                        ax2.set_ylabel('Y coordinate')
                        ax2.set_title(f'Expression Fraction ({direction})\n{gene_set_name}\n'
                                     f'Range: {valid_fractions.min():.2%} to {valid_fractions.max():.2%}')
                        ax2.set_aspect('equal')

                fig.suptitle(f'Sample {sample_id}', fontsize=16, y=0.98)

                plot_filename = f"{clean_name}_{direction}_with_fraction_sample_{sample_id}.png"
                fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight', dpi=300)
                plt.close(fig)

            print(f"  Saved {direction} spatial plots")

def create_diagnostic_plots(fraction_data, output_dir):
    """Create diagnostic plots for each direction."""
    print("\nCreating diagnostic plots...")

    df = pd.DataFrame(fraction_data)
    if df.empty:
        return

    for direction in ['up', 'down']:
        dir_df = df[df['Direction'] == direction]
        if dir_df.empty:
            continue

        for gene_set in dir_df['Gene_Set'].unique():
            gs_df = dir_df[dir_df['Gene_Set'] == gene_set]

            fig, ax = plt.subplots(figsize=(10, 8))
            scatter = ax.scatter(gs_df['Fraction_Expressed'], gs_df['Score'],
                               c=gs_df['N_Expressed'], cmap='viridis', s=100, alpha=0.7)

            for _, row in gs_df.iterrows():
                ax.annotate(row['Cluster'], (row['Fraction_Expressed'], row['Score']),
                          fontsize=8, alpha=0.7)

            plt.colorbar(scatter, ax=ax, label='N Genes Expressed')
            ax.set_xlabel('Fraction Expressed')
            ax.set_ylabel('Score')
            ax.set_title(f'{gene_set} - {direction.upper()}\nScore vs Expression Fraction')
            ax.grid(True, alpha=0.3)

            clean_name = gene_set.replace(' ', '_').replace('/', '_')
            plt.savefig(os.path.join(output_dir, f'diagnostic_{clean_name}_{direction}.png'),
                       bbox_inches='tight', dpi=300)
            plt.close()

def main():
    print("=" * 80)
    print("SOLUTION 2: EXPRESSION FRACTION WITH DIRECTIONAL GENE SETS")
    print("=" * 80)

    adata = sc.read_h5ad(merged_h5ad_path)
    gene_sets_dict = load_directional_gene_sets(gene_sets_path)
    pseudobulk_adata = create_pseudobulk(adata)

    pseudobulk_adata, fraction_data = calculate_scores_with_fraction(
        pseudobulk_adata, gene_sets_dict, EXPRESSION_THRESHOLD
    )

    pd.DataFrame(fraction_data).to_csv(
        os.path.join(output_dir, 'expression_fraction_metrics.csv'), index=False
    )

    adata_with_metrics = assign_metrics_to_cells(adata, pseudobulk_adata, gene_sets_dict)

    # Create spatial plots with fraction
    create_spatial_plots_with_fraction(adata_with_metrics, gene_sets_dict, output_dir)

    # Create diagnostic plots
    create_diagnostic_plots(fraction_data, output_dir)

    print("\n" + "=" * 80)
    print("SOLUTION 2 COMPLETED")
    print(f"Results: {output_dir}")
    print("\nGenerated outputs:")
    print("  - Spatial plots: *_{up/down/aggregate}_with_fraction_sample_*.png")
    print("  - Diagnostic plots: diagnostic_*_{up/down}.png")
    print("  - Metrics CSV: expression_fraction_metrics.csv")
    print("=" * 80)

if __name__ == "__main__":
    main()
