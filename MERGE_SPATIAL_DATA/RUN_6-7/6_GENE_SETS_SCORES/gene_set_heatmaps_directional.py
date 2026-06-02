#!/usr/bin/env python3
"""
Directional Gene Set Scoring and Heatmap Visualization

This script performs directional gene set scoring at single cell level and creates heatmap
visualizations showing gene set scores across different clusters, conditions, regions, and
time points. Each gene set includes genes marked as positively or negatively correlated
with a biological process.

Three types of scores are calculated:
1. Positive scores (using only upregulated genes)
2. Negative scores (using only downregulated genes)
3. Aggregate scores (combining both with correct signs)

Additionally, it creates temporal heatmaps showing how gene scores change across time points
(e16, p0, p7, p30) for different conditions, regions, and clusters.
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
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/gene_set_heatmaps_directional'

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
                'all_genes': [all genes]
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

def calculate_directional_gene_set_scores(adata, gene_sets_dict):
    """
    Calculate directional gene set scores using scanpy's score_genes function.

    For each gene set, calculates:
    - Positive score (using only upregulated genes)
    - Negative score (using only downregulated genes)
    - Aggregate score (upregulated genes positively, downregulated genes negatively)

    Parameters:
    -----------
    adata : AnnData
        Annotated data object containing the gene expression data
    gene_sets_dict : dict
        Dictionary with gene set information including directional genes

    Returns:
    --------
    AnnData
        Updated AnnData object with gene set scores in .obs
    """
    print("Calculating directional gene set scores...")

    # Make a copy to avoid modifying the original
    adata_copy = adata.copy()

    # Calculate scores for each gene set
    for gene_set_name, gene_info in gene_sets_dict.items():
        print(f"  Scoring gene set: {gene_set_name}")

        # Process upregulated genes
        up_genes = gene_info['up']
        up_genes_in_data = [gene for gene in up_genes if gene in adata_copy.var.index]
        up_genes_missing = [gene for gene in up_genes if gene not in adata_copy.var.index]

        if len(up_genes_missing) > 0:
            print(f"    Warning: {len(up_genes_missing)}/{len(up_genes)} upregulated genes not found")

        if len(up_genes_in_data) > 0:
            print(f"    Positive: Using {len(up_genes_in_data)}/{len(up_genes)} upregulated genes")
            sc.tl.score_genes(adata_copy, gene_list=up_genes_in_data,
                            score_name=f"{gene_set_name}_positive_score", use_raw=False)
        else:
            print(f"    Warning: No upregulated genes found in dataset")
            adata_copy.obs[f"{gene_set_name}_positive_score"] = 0.0

        # Process downregulated genes
        down_genes = gene_info['down']
        down_genes_in_data = [gene for gene in down_genes if gene in adata_copy.var.index]
        down_genes_missing = [gene for gene in down_genes if gene not in adata_copy.var.index]

        if len(down_genes_missing) > 0:
            print(f"    Warning: {len(down_genes_missing)}/{len(down_genes)} downregulated genes not found")

        if len(down_genes_in_data) > 0:
            print(f"    Negative: Using {len(down_genes_in_data)}/{len(down_genes)} downregulated genes")
            sc.tl.score_genes(adata_copy, gene_list=down_genes_in_data,
                            score_name=f"{gene_set_name}_negative_score", use_raw=False)
        else:
            print(f"    Warning: No downregulated genes found in dataset")
            adata_copy.obs[f"{gene_set_name}_negative_score"] = 0.0

        # Calculate aggregate score
        # Aggregate = positive contribution from upregulated - positive contribution from downregulated
        # (high expression of up genes is good, high expression of down genes is bad)
        if f"{gene_set_name}_positive_score" in adata_copy.obs.columns and \
           f"{gene_set_name}_negative_score" in adata_copy.obs.columns:
            adata_copy.obs[f"{gene_set_name}_aggregate_score"] = \
                adata_copy.obs[f"{gene_set_name}_positive_score"] - \
                adata_copy.obs[f"{gene_set_name}_negative_score"]
            print(f"    Aggregate score calculated")

    return adata_copy

def create_cluster_heatmap(adata, gene_sets_dict, time_point, score_type, output_dir):
    """
    Create heatmap showing gene set scores for leiden_new clusters within a given time point,
    across different conditions (ctrl and mut) and regions.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names
    time_point : str
        Time point to analyze (e.g., 'e16', 'p0', 'p7', 'p30')
    score_type : str
        Type of score ('positive', 'negative', or 'aggregate')
    output_dir : str
        Directory to save the plots
    """
    print(f"Creating cluster heatmap for time point: {time_point}, score type: {score_type}")

    # Filter data for the specified time point
    time_adata = adata[adata.obs['time_point'] == time_point].copy()

    if len(time_adata) == 0:
        print(f"  No data found for time point: {time_point}")
        return

    # Get unique conditions and regions
    conditions = time_adata.obs['type'].unique()
    regions = time_adata.obs['region'].unique()

    # Collect all data to calculate min/max for consistent color scaling
    all_scores_for_timepoint = []
    all_heatmap_data = {}

    for region in regions:
        heatmap_data = []
        region_samples = [s for s in time_adata.obs['sample_id'].unique() if region in str(s)]

        for condition in conditions:
            condition_samples = [s for s in region_samples if condition in str(s)]

            if not condition_samples:
                continue

            condition_adata = time_adata[time_adata.obs['sample_id'].isin(condition_samples)]

            for cluster in sorted(condition_adata.obs['leiden_new'].unique()):
                cluster_data = condition_adata[condition_adata.obs['leiden_new'] == cluster]

                if len(cluster_data) == 0:
                    continue

                row_data = {
                    'Condition': condition,
                    'Region': region,
                    'Cluster': cluster
                }

                for gene_set_name in gene_sets_dict.keys():
                    score_column = f"{gene_set_name}_{score_type}_score"
                    if score_column in cluster_data.obs.columns:
                        score = cluster_data.obs[score_column].mean()
                        row_data[gene_set_name] = score
                        all_scores_for_timepoint.append(score)

                heatmap_data.append(row_data)

        all_heatmap_data[region] = heatmap_data

    # Calculate min and max for this time point
    if all_scores_for_timepoint:
        file_min = min(all_scores_for_timepoint)
        file_max = max(all_scores_for_timepoint)
        print(f"  Using color scale range: [{file_min:.3f}, {file_max:.3f}]")
    else:
        file_min = -1
        file_max = 1

    # Create figure with subplots for each region
    fig, axes = plt.subplots(len(regions), 1, figsize=(12, 8 * len(regions)))
    if len(regions) == 1:
        axes = [axes]

    # Choose colormap based on score type
    if score_type == 'positive':
        cmap = 'Reds'
        center_val = None
    elif score_type == 'negative':
        cmap = 'Blues'
        center_val = None
    else:  # aggregate
        cmap = 'coolwarm'
        center_val = 0

    for region_idx, region in enumerate(regions):
        ax = axes[region_idx]
        heatmap_data = all_heatmap_data[region]

        if not heatmap_data:
            continue

        heatmap_df = pd.DataFrame(heatmap_data)
        plot_data = heatmap_df.pivot_table(
            index=['Condition', 'Region', 'Cluster'],
            values=[col for col in heatmap_df.columns if col not in ['Condition', 'Region', 'Cluster']],
            aggfunc='mean'
        )

        sns.heatmap(plot_data, ax=ax, cmap=cmap, center=center_val,
                   vmin=file_min, vmax=file_max,
                   cbar_kws={'label': f'{score_type.capitalize()} Gene Set Score'})

        ax.set_title(f'{score_type.capitalize()} Gene Set Scores by Cluster - Region {region}\n'
                    f'Time Point: {time_point}')
        ax.set_xlabel('Gene Sets')
        ax.set_ylabel('Condition, Region, Cluster')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'cluster_heatmap_{score_type}_timepoint_{time_point}.png'),
               bbox_inches='tight', dpi=300)
    plt.close()

    print(f"  Saved: cluster_heatmap_{score_type}_timepoint_{time_point}.png")

def create_celltype_heatmap(adata, gene_sets_dict, cluster, score_type, output_dir):
    """
    Create heatmap showing gene set scores for a given leiden_new cluster across different samples,
    conditions, regions, and time points.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names
    cluster : str
        Leiden cluster to analyze
    score_type : str
        Type of score ('positive', 'negative', or 'aggregate')
    output_dir : str
        Directory to save the plots
    """
    print(f"Creating cluster heatmap for: {cluster}, score type: {score_type}")

    # Filter data for the specified cluster
    cluster_adata = adata[adata.obs['leiden_new'] == cluster].copy()

    if len(cluster_adata) == 0:
        print(f"  No data found for cluster: {cluster}")
        return

    # Define the desired time point order
    time_point_order = ["e16", "p0", "p7", "p30"]

    # Get unique conditions and regions
    conditions = cluster_adata.obs['type'].unique()
    regions = cluster_adata.obs['region'].unique()

    # Prepare data for heatmap
    heatmap_data = []
    all_scores_for_cluster = []

    for time_point in time_point_order:
        if time_point not in cluster_adata.obs['time_point'].unique():
            continue

        for condition in conditions:
            for region in regions:
                samples = [s for s in cluster_adata.obs['sample_id'].unique()
                          if time_point in str(s) and condition in str(s) and region in str(s)]

                if not samples:
                    continue

                combo_adata = cluster_adata[cluster_adata.obs['sample_id'].isin(samples)]

                if len(combo_adata) == 0:
                    continue

                row_data = {
                    'Time_Point': time_point,
                    'Condition': condition,
                    'Region': region
                }

                for gene_set_name in gene_sets_dict.keys():
                    score_column = f"{gene_set_name}_{score_type}_score"
                    if score_column in combo_adata.obs.columns:
                        score = combo_adata.obs[score_column].mean()
                        row_data[gene_set_name] = score
                        all_scores_for_cluster.append(score)

                heatmap_data.append(row_data)

    # Calculate min and max for this cluster
    if all_scores_for_cluster:
        file_min = min(all_scores_for_cluster)
        file_max = max(all_scores_for_cluster)
    else:
        file_min = -1
        file_max = 1

    # Convert to DataFrame
    heatmap_df = pd.DataFrame(heatmap_data)

    if heatmap_df.empty:
        print(f"  No data available for cluster: {cluster}")
        return

    # Convert Time_Point to categorical with the correct order
    heatmap_df['Time_Point'] = pd.Categorical(
        heatmap_df['Time_Point'],
        categories=time_point_order,
        ordered=True
    )

    # Set up the data for heatmap
    plot_data = heatmap_df.pivot_table(
        index=['Time_Point', 'Condition', 'Region'],
        values=[col for col in heatmap_df.columns if col not in ['Time_Point', 'Condition', 'Region']],
        aggfunc='mean'
    )

    # Choose colormap based on score type
    if score_type == 'positive':
        cmap = 'Reds'
        center_val = None
    elif score_type == 'negative':
        cmap = 'Blues'
        center_val = None
    else:  # aggregate
        cmap = 'coolwarm'
        center_val = 0

    # Create heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(plot_data, cmap=cmap, center=center_val,
               vmin=file_min, vmax=file_max,
               cbar_kws={'label': f'{score_type.capitalize()} Gene Set Score'})

    plt.title(f'{score_type.capitalize()} Gene Set Scores for Cluster {cluster}\n'
             f'Across Time Points, Conditions, and Regions')
    plt.xlabel('Gene Sets')
    plt.ylabel('Time Point, Condition, Region')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'cluster_heatmap_{cluster}_{score_type}.png'),
               bbox_inches='tight', dpi=300)
    plt.close()

    print(f"  Saved: cluster_heatmap_{cluster}_{score_type}.png")

def create_temporal_heatmap(adata, gene_sets_dict, score_type, output_dir):
    """
    Create heatmap showing gene set scores with time points, conditions, and regions on x-axis
    and clusters on y-axis. One plot per gene set, combining all regions.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names
    score_type : str
        Type of score ('positive', 'negative', or 'aggregate')
    output_dir : str
        Directory to save the plots
    """
    print(f"Creating temporal heatmaps for {score_type} scores...")

    # Define the desired time point order
    time_point_order = ["e16", "p0", "p7", "p30"]

    # Get unique regions and clusters
    regions = sorted(adata.obs['region'].unique())
    clusters = sorted(adata.obs['leiden_new'].unique())
    conditions = sorted(adata.obs['type'].unique())

    # First pass: calculate global min and max for consistent color scaling
    all_scores = []
    for gene_set_name in gene_sets_dict.keys():
        score_column = f"{gene_set_name}_{score_type}_score"

        if score_column not in adata.obs.columns:
            continue

        gene_set_scores = adata.obs[score_column].dropna()
        all_scores.extend(gene_set_scores.tolist())

    # Calculate global min and max
    if all_scores:
        global_min = min(all_scores)
        global_max = max(all_scores)
        print(f"  Global score range: [{global_min:.3f}, {global_max:.3f}]")
    else:
        global_min = -1
        global_max = 1

    # Choose colormap based on score type
    if score_type == 'positive':
        cmap = 'Reds'
        center_val = None
    elif score_type == 'negative':
        cmap = 'Blues'
        center_val = None
    else:  # aggregate
        cmap = 'coolwarm'
        center_val = 0

    # Create a separate heatmap for each gene set
    for gene_set_name in gene_sets_dict.keys():
        score_column = f"{gene_set_name}_{score_type}_score"

        if score_column not in adata.obs.columns:
            continue

        print(f"  Processing gene set: {gene_set_name}")

        # Prepare data for heatmap
        heatmap_data = []

        for cluster in clusters:
            cluster_adata = adata[adata.obs['leiden_new'] == cluster]

            if len(cluster_adata) == 0:
                continue

            for time_point in time_point_order:
                for condition in conditions:
                    for region in regions:
                        combo_adata = cluster_adata[
                            (cluster_adata.obs['time_point'] == time_point) &
                            (cluster_adata.obs['type'] == condition) &
                            (cluster_adata.obs['region'] == region)
                        ]

                        if len(combo_adata) == 0:
                            heatmap_data.append({
                                'Cluster': cluster,
                                'Time_Condition_Region': f"{time_point}_{condition}_{region}",
                                'Score': np.nan
                            })
                        else:
                            mean_score = combo_adata.obs[score_column].mean()
                            heatmap_data.append({
                                'Cluster': cluster,
                                'Time_Condition_Region': f"{time_point}_{condition}_{region}",
                                'Score': mean_score
                            })

        # Convert to DataFrame
        heatmap_df = pd.DataFrame(heatmap_data)

        if heatmap_df.empty:
            continue

        # Pivot to create matrix for heatmap
        plot_data = heatmap_df.pivot(
            index='Cluster',
            columns='Time_Condition_Region',
            values='Score'
        )

        # Reorder columns
        ordered_columns = [f"{tp}_{cond}_{region}" for tp in time_point_order
                          for cond in conditions for region in regions]
        ordered_columns = [col for col in ordered_columns if col in plot_data.columns]
        plot_data = plot_data[ordered_columns]

        # Create heatmap
        plt.figure(figsize=(24, 12))
        sns.heatmap(plot_data, cmap=cmap, center=center_val,
                   vmin=global_min, vmax=global_max,
                   cbar_kws={'label': f'{score_type.capitalize()} Gene Set Score'},
                   annot=False, fmt='.2f')

        plt.title(f'{score_type.capitalize()} Gene Set Scores: {gene_set_name}\nAll Regions Combined')
        plt.xlabel('Time Point & Condition & Region')
        plt.ylabel('Cluster')
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        plt.tight_layout()

        safe_gene_set = gene_set_name.replace('/', '_').replace('\\', '_').replace(' ', '_')
        filename = f'temporal_heatmap_{score_type}_all_regions_{safe_gene_set}.png'
        plt.savefig(os.path.join(output_dir, filename), bbox_inches='tight', dpi=300)
        plt.close()

        print(f"  Saved: {filename}")

def create_summary_report(adata, gene_sets_dict, output_dir):
    """
    Create a summary report of directional gene set scores.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names
    output_dir : str
        Directory to save the report
    """
    print("Creating summary report...")

    # Create summary for each score type
    for score_type in ['positive', 'negative', 'aggregate']:
        summary_data = []

        for gene_set_name in gene_sets_dict.keys():
            score_column = f"{gene_set_name}_{score_type}_score"

            if score_column not in adata.obs.columns:
                continue

            scores = adata.obs[score_column]

            # Overall statistics
            overall_stats = {
                'Gene_Set': gene_set_name,
                'Score_Type': score_type,
                'Sample': 'All',
                'N_Cells': len(scores),
                'Mean_Score': scores.mean(),
                'Std_Score': scores.std(),
                'Min_Score': scores.min(),
                'Max_Score': scores.max(),
                'N_Genes_Up': len(gene_sets_dict[gene_set_name]['up']),
                'N_Genes_Down': len(gene_sets_dict[gene_set_name]['down']),
                'N_Genes_Up_Found': len([g for g in gene_sets_dict[gene_set_name]['up'] if g in adata.var.index]),
                'N_Genes_Down_Found': len([g for g in gene_sets_dict[gene_set_name]['down'] if g in adata.var.index])
            }
            summary_data.append(overall_stats)

            # Per-sample statistics
            for sample_id in adata.obs['sample_id'].unique():
                sample_scores = adata.obs[adata.obs['sample_id'] == sample_id][score_column]

                sample_stats = {
                    'Gene_Set': gene_set_name,
                    'Score_Type': score_type,
                    'Sample': sample_id,
                    'N_Cells': len(sample_scores),
                    'Mean_Score': sample_scores.mean(),
                    'Std_Score': sample_scores.std(),
                    'Min_Score': sample_scores.min(),
                    'Max_Score': sample_scores.max(),
                    'N_Genes_Up': len(gene_sets_dict[gene_set_name]['up']),
                    'N_Genes_Down': len(gene_sets_dict[gene_set_name]['down']),
                    'N_Genes_Up_Found': len([g for g in gene_sets_dict[gene_set_name]['up'] if g in adata.var.index]),
                    'N_Genes_Down_Found': len([g for g in gene_sets_dict[gene_set_name]['down'] if g in adata.var.index])
                }
                summary_data.append(sample_stats)

        # Convert to DataFrame and save
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(os.path.join(output_dir, f'gene_set_{score_type}_scores_summary.csv'), index=False)
        print(f"  Saved: gene_set_{score_type}_scores_summary.csv")

    # Create a formatted text report
    with open(os.path.join(output_dir, 'gene_set_scores_report.txt'), 'w') as f:
        f.write("Directional Gene Set Scoring Summary Report\n")
        f.write("=" * 50 + "\n\n")

        f.write(f"Total gene sets analyzed: {len(gene_sets_dict)}\n")
        f.write(f"Total cells: {len(adata)}\n")
        f.write(f"Total samples: {len(adata.obs['sample_id'].unique())}\n\n")

        f.write("Gene Sets Overview:\n")
        f.write("-" * 20 + "\n")
        for gene_set_name, gene_info in gene_sets_dict.items():
            n_up = len(gene_info['up'])
            n_down = len(gene_info['down'])
            n_up_found = len([g for g in gene_info['up'] if g in adata.var.index])
            n_down_found = len([g for g in gene_info['down'] if g in adata.var.index])
            f.write(f"{gene_set_name}:\n")
            f.write(f"  Upregulated: {n_up_found}/{n_up} genes found\n")
            f.write(f"  Downregulated: {n_down_found}/{n_down} genes found\n")

        f.write(f"\nDetailed statistics saved to: gene_set_*_scores_summary.csv\n")
        f.write(f"Heatmaps saved to: {output_dir}\n")

def main():
    """Main function to run the directional gene set scoring and heatmap visualization pipeline."""

    print("Directional Gene Set Scoring and Heatmap Visualization Pipeline")
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

    # Calculate directional gene set scores
    adata_scored = calculate_directional_gene_set_scores(adata, gene_sets_dict)

    # Define score types
    score_types = ['positive', 'negative', 'aggregate']

    # Create cluster heatmaps for each time point and score type
    time_point_order = ["e16", "p0", "p7", "p30"]
    time_points = [tp for tp in time_point_order if tp in adata_scored.obs['time_point'].unique()]

    for score_type in score_types:
        print(f"\n{'=' * 80}")
        print(f"Creating heatmaps for {score_type} scores")
        print(f"{'=' * 80}")

        for time_point in time_points:
            create_cluster_heatmap(adata_scored, gene_sets_dict, time_point, score_type, output_dir)

        # Create cell type heatmaps for each cluster
        unique_clusters = sorted(adata_scored.obs['leiden_new'].unique())
        for cluster in unique_clusters:
            create_celltype_heatmap(adata_scored, gene_sets_dict, cluster, score_type, output_dir)

        # Create temporal heatmaps
        # create_temporal_heatmap(adata_scored, gene_sets_dict, score_type, output_dir)

    # Create summary report
    create_summary_report(adata_scored, gene_sets_dict, output_dir)

    print("\n" + "=" * 80)
    print("Pipeline completed successfully!")
    print(f"Results saved to: {output_dir}")
    print("  - Heatmaps for positive, negative, and aggregate scores")
    print("  - Summary reports for each score type")
    print("=" * 80)

if __name__ == "__main__":
    main()
