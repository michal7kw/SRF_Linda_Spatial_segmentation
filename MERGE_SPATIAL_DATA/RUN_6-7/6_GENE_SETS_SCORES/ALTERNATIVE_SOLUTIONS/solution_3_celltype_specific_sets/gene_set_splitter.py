#!/usr/bin/env python3
"""
Solution 3: Cell-Type Specific Gene Sets

This is a HELPER TOOL to split your gene sets into cell-type-specific subsets
based on expression patterns. Instead of scoring "Reelin pathway" with all genes
in all cell types, you create separate gene sets like:
- Reelin_receptors
- Reelin_ligands
- Reelin_signaling_GC (granule cell specific)

This tool analyzes your data and suggests how to split gene sets based on
which genes are co-expressed in which clusters.

Key features:
- Analyzes expression patterns across clusters
- Identifies which genes cluster together by cell type
- Suggests cell-type-specific gene set splits
- Generates new gene set files ready for scoring

Use case: When pathways have modular, cell-type-specific components that
          should be scored separately (ligands vs receptors, etc.)
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
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
import warnings
warnings.filterwarnings('ignore')

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Define paths
merged_h5ad_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
gene_sets_path = 'DATA/GENE_LISTS/output_combined/import_ready.csv'
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_3_celltype_specific_sets/results'

# Create output directory
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

# Parameters
EXPRESSION_THRESHOLD = 0.5  # Threshold for considering gene expressed
MIN_CLUSTER_SPECIFICITY = 0.3  # Minimum fraction of clusters where gene should be expressed

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

def analyze_gene_expression_patterns(pseudobulk_adata, genes, gene_set_name, expr_threshold=0.5):
    """
    Analyze which genes are expressed in which clusters.

    Returns expression patterns and suggests potential splits.
    """
    print(f"\nAnalyzing expression patterns for: {gene_set_name}")

    genes_in_data = [g for g in genes if g in pseudobulk_adata.var.index]

    if len(genes_in_data) == 0:
        print(f"  No genes found in dataset")
        return None

    print(f"  Analyzing {len(genes_in_data)}/{len(genes)} genes")

    # Get expression matrix for these genes across clusters
    gene_indices = [pseudobulk_adata.var.index.get_loc(g) for g in genes_in_data]
    expr_matrix = pseudobulk_adata.X[:, gene_indices]  # clusters x genes

    # Create binary expression matrix (expressed or not)
    binary_expr = (expr_matrix > expr_threshold).astype(int)

    # Create DataFrame
    expr_df = pd.DataFrame(
        expr_matrix,
        index=pseudobulk_adata.obs.index,
        columns=genes_in_data
    )

    binary_df = pd.DataFrame(
        binary_expr,
        index=pseudobulk_adata.obs.index,
        columns=genes_in_data
    )

    # Calculate gene statistics
    gene_stats = []
    for gene in genes_in_data:
        n_clusters_expressed = binary_df[gene].sum()
        frac_clusters_expressed = n_clusters_expressed / len(binary_df)
        mean_expr = expr_df[gene].mean()
        max_expr = expr_df[gene].max()

        # Find which clusters express this gene
        expressing_clusters = binary_df[binary_df[gene] == 1].index.tolist()

        gene_stats.append({
            'Gene': gene,
            'N_Clusters_Expressed': n_clusters_expressed,
            'Frac_Clusters_Expressed': frac_clusters_expressed,
            'Mean_Expression': mean_expr,
            'Max_Expression': max_expr,
            'Expressing_Clusters': ', '.join(expressing_clusters)
        })

    gene_stats_df = pd.DataFrame(gene_stats)

    return expr_df, binary_df, gene_stats_df

def cluster_genes_by_expression_pattern(binary_df, gene_stats_df, gene_set_name, output_dir):
    """
    Cluster genes based on their expression patterns across clusters.
    Genes with similar expression patterns will cluster together.
    """
    print(f"\nClustering genes by expression pattern for: {gene_set_name}")

    if len(binary_df.columns) < 2:
        print(f"  Not enough genes to cluster")
        return None

    # Transpose: genes x clusters
    gene_expr_patterns = binary_df.T

    # Calculate correlation between genes
    gene_corr = gene_expr_patterns.T.corr()

    # Hierarchical clustering
    linkage_matrix = linkage(gene_expr_patterns, method='ward')

    # Create dendrogram
    fig, ax = plt.subplots(figsize=(12, 8))
    dendrogram(linkage_matrix, labels=gene_expr_patterns.index, ax=ax, orientation='right')
    ax.set_xlabel('Distance')
    ax.set_ylabel('Genes')
    ax.set_title(f'Gene Clustering by Expression Pattern\n{gene_set_name}')
    plt.tight_layout()

    clean_name = gene_set_name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
    plt.savefig(os.path.join(output_dir, f'gene_clustering_{clean_name}.png'),
               bbox_inches='tight', dpi=300)
    plt.close()

    print(f"  Saved dendrogram: gene_clustering_{clean_name}.png")

    # Perform clustering (find optimal number of clusters)
    n_genes = len(gene_expr_patterns)
    max_clusters = min(5, n_genes - 1)  # At most 5 sub-gene-sets

    if max_clusters < 2:
        print(f"  Not enough genes for meaningful clustering")
        return None

    clustering = AgglomerativeClustering(n_clusters=max_clusters, linkage='ward')
    cluster_labels = clustering.fit_predict(gene_expr_patterns)

    # Add cluster labels to gene stats
    gene_stats_with_clusters = gene_stats_df.copy()
    gene_stats_with_clusters['Expression_Cluster'] = cluster_labels

    return gene_stats_with_clusters, cluster_labels

def create_expression_heatmap(expr_df, binary_df, gene_stats_with_clusters, gene_set_name, output_dir):
    """Create heatmap showing gene expression across clusters."""
    print(f"\nCreating expression heatmap for: {gene_set_name}")

    # Sort genes by expression cluster if available
    if 'Expression_Cluster' in gene_stats_with_clusters.columns:
        gene_order = gene_stats_with_clusters.sort_values('Expression_Cluster')['Gene'].tolist()
        expr_df_sorted = expr_df[gene_order]
        binary_df_sorted = binary_df[gene_order]
    else:
        expr_df_sorted = expr_df
        binary_df_sorted = binary_df

    # Create figure with two heatmaps
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))

    # Heatmap 1: Continuous expression
    sns.heatmap(expr_df_sorted.T, cmap='viridis', ax=ax1,
               cbar_kws={'label': 'Expression Level'},
               yticklabels=True, xticklabels=True)
    ax1.set_title(f'Gene Expression Levels\n{gene_set_name}')
    ax1.set_xlabel('Clusters')
    ax1.set_ylabel('Genes')

    # Heatmap 2: Binary expression
    sns.heatmap(binary_df_sorted.T, cmap='RdYlGn', vmin=0, vmax=1, ax=ax2,
               cbar_kws={'label': 'Expressed (1) / Not Expressed (0)'},
               yticklabels=True, xticklabels=True)
    ax2.set_title(f'Gene Expression Patterns (Binary)\n{gene_set_name}')
    ax2.set_xlabel('Clusters')
    ax2.set_ylabel('Genes')

    # Add expression cluster annotations if available
    if 'Expression_Cluster' in gene_stats_with_clusters.columns:
        for ax in [ax1, ax2]:
            # Add vertical lines to separate expression clusters
            gene_order_clusters = gene_stats_with_clusters.sort_values('Expression_Cluster')['Expression_Cluster'].tolist()
            cluster_changes = [i for i in range(1, len(gene_order_clusters))
                             if gene_order_clusters[i] != gene_order_clusters[i-1]]
            for pos in cluster_changes:
                ax.axhline(y=pos, color='red', linewidth=2, linestyle='--')

    plt.tight_layout()

    clean_name = gene_set_name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
    plt.savefig(os.path.join(output_dir, f'expression_heatmap_{clean_name}.png'),
               bbox_inches='tight', dpi=300)
    plt.close()

    print(f"  Saved heatmap: expression_heatmap_{clean_name}.png")

def suggest_gene_set_splits(gene_stats_with_clusters, gene_set_name):
    """
    Suggest how to split the gene set based on expression clustering.
    """
    print(f"\nSuggested gene set splits for: {gene_set_name}")
    print("=" * 60)

    if 'Expression_Cluster' not in gene_stats_with_clusters.columns:
        print("  No clustering performed (not enough genes)")
        return []

    splits = []

    for cluster_id in sorted(gene_stats_with_clusters['Expression_Cluster'].unique()):
        cluster_genes = gene_stats_with_clusters[
            gene_stats_with_clusters['Expression_Cluster'] == cluster_id
        ]

        # Identify predominant clusters where these genes are expressed
        all_expressing_clusters = []
        for clusters_str in cluster_genes['Expressing_Clusters']:
            all_expressing_clusters.extend(clusters_str.split(', '))

        from collections import Counter
        cluster_counts = Counter(all_expressing_clusters)
        top_clusters = cluster_counts.most_common(3)

        # Create suggested name
        suggested_name = f"{gene_set_name}_subset{cluster_id+1}"
        if top_clusters:
            top_cluster_names = '_'.join([c[0] for c in top_clusters[:2]])
            suggested_name = f"{gene_set_name}_{top_cluster_names}"

        split_info = {
            'Original_Gene_Set': gene_set_name,
            'Suggested_Name': suggested_name,
            'Genes': cluster_genes['Gene'].tolist(),
            'N_Genes': len(cluster_genes),
            'Primary_Clusters': ', '.join([f"{c[0]} ({c[1]})" for c in top_clusters[:3]])
        }

        splits.append(split_info)

        print(f"\nSubset {cluster_id + 1}: {suggested_name}")
        print(f"  Genes ({len(cluster_genes)}): {', '.join(cluster_genes['Gene'].tolist())}")
        print(f"  Primarily expressed in: {split_info['Primary_Clusters']}")

    print("\n" + "=" * 60)

    return splits

def save_split_gene_sets(all_splits, output_dir):
    """
    Save suggested gene set splits to CSV files.
    Creates files in the same format as input gene sets.
    """
    print("\nSaving suggested gene set splits...")

    # Create DataFrame in import_ready format
    split_records = []
    for split in all_splits:
        for gene in split['Genes']:
            split_records.append({
                'gene_set': split['Suggested_Name'],
                'gene': gene,
                'col3': '',
                'col4': ''
            })

    if split_records:
        split_df = pd.DataFrame(split_records)
        output_file = os.path.join(output_dir, 'suggested_gene_set_splits_import_ready.csv')
        split_df.to_csv(output_file, index=False, header=False)
        print(f"  Saved: suggested_gene_set_splits_import_ready.csv")
        print(f"  This file can be used directly as input for gene set scoring!")

    # Also save split summary
    split_summary = []
    for split in all_splits:
        split_summary.append({
            'Original_Gene_Set': split['Original_Gene_Set'],
            'Suggested_Name': split['Suggested_Name'],
            'N_Genes': split['N_Genes'],
            'Genes': ', '.join(split['Genes']),
            'Primary_Clusters': split['Primary_Clusters']
        })

    if split_summary:
        summary_df = pd.DataFrame(split_summary)
        summary_file = os.path.join(output_dir, 'gene_set_split_summary.csv')
        summary_df.to_csv(summary_file, index=False)
        print(f"  Saved: gene_set_split_summary.csv")

def main():
    """Main function."""

    print("=" * 80)
    print("SOLUTION 3: CELL-TYPE SPECIFIC GENE SETS")
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

    # Analyze each gene set
    all_splits = []
    all_gene_stats = []

    for gene_set_name, genes in gene_sets_dict.items():
        print("\n" + "=" * 80)

        # Analyze expression patterns
        result = analyze_gene_expression_patterns(
            pseudobulk_adata, genes, gene_set_name, expr_threshold=EXPRESSION_THRESHOLD
        )

        if result is None:
            continue

        expr_df, binary_df, gene_stats_df = result

        # Save gene statistics
        gene_stats_df['Gene_Set'] = gene_set_name
        all_gene_stats.append(gene_stats_df)

        # Cluster genes by expression pattern
        clustering_result = cluster_genes_by_expression_pattern(
            binary_df, gene_stats_df, gene_set_name, output_dir
        )

        if clustering_result is not None:
            gene_stats_with_clusters, cluster_labels = clustering_result

            # Create heatmap
            create_expression_heatmap(
                expr_df, binary_df, gene_stats_with_clusters, gene_set_name, output_dir
            )

            # Suggest splits
            splits = suggest_gene_set_splits(gene_stats_with_clusters, gene_set_name)
            all_splits.extend(splits)
        else:
            # No clustering, but still create heatmap
            create_expression_heatmap(
                expr_df, binary_df, gene_stats_df, gene_set_name, output_dir
            )

    # Save all gene statistics
    if all_gene_stats:
        combined_stats = pd.concat(all_gene_stats, ignore_index=True)
        combined_stats.to_csv(os.path.join(output_dir, 'all_gene_expression_statistics.csv'),
                            index=False)
        print(f"\nSaved: all_gene_expression_statistics.csv")

    # Save suggested splits
    if all_splits:
        save_split_gene_sets(all_splits, output_dir)

    print("\n" + "=" * 80)
    print("SOLUTION 3 COMPLETED")
    print(f"Results saved to: {output_dir}")
    print("\nNext steps:")
    print("  1. Review the suggested gene set splits in: gene_set_split_summary.csv")
    print("  2. Examine expression heatmaps to understand patterns")
    print("  3. Manually refine splits based on biological knowledge")
    print("  4. Use suggested_gene_set_splits_import_ready.csv for scoring")
    print("=" * 80)

if __name__ == "__main__":
    main()
