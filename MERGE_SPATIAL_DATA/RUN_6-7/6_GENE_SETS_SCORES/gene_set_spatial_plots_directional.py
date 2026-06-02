#!/usr/bin/env python3
"""
Directional Gene Set Scoring and Spatial Visualization

This script performs directional gene set scoring at single cell level and creates spatial
visualizations. Each gene set includes genes marked as positively or negatively correlated
with a biological process.

Three types of scores are calculated:
1. Positive scores (using only upregulated genes)
2. Negative scores (using only downregulated genes)
3. Aggregate scores (combining both with correct signs)

The workflow:
1. Load spatial transcriptomics data
2. Load gene sets with directional information
3. Calculate gene set scores for each cell (positive, negative, aggregate)
4. Create spatial plots showing the distribution of scores
"""

import scanpy as sc
import squidpy as sq
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/gene_set_spatial_plots_directional'

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

    # Calculate scores for each gene set
    for gene_set_name, gene_info in gene_sets_dict.items():
        print(f"  Scoring gene set: {gene_set_name}")

        # Process upregulated genes
        up_genes = gene_info['up']
        up_genes_in_data = [gene for gene in up_genes if gene in adata.var.index]
        up_genes_missing = [gene for gene in up_genes if gene not in adata.var.index]

        if len(up_genes_missing) > 0:
            print(f"    Warning: {len(up_genes_missing)}/{len(up_genes)} upregulated genes not found")

        if len(up_genes_in_data) > 0:
            print(f"    Positive: Using {len(up_genes_in_data)}/{len(up_genes)} upregulated genes")
            sc.tl.score_genes(adata, gene_list=up_genes_in_data,
                            score_name=f"{gene_set_name}_positive_score", use_raw=False)
        else:
            print(f"    Warning: No upregulated genes found in dataset")
            adata.obs[f"{gene_set_name}_positive_score"] = 0.0

        # Process downregulated genes
        down_genes = gene_info['down']
        down_genes_in_data = [gene for gene in down_genes if gene in adata.var.index]
        down_genes_missing = [gene for gene in down_genes if gene not in adata.var.index]

        if len(down_genes_missing) > 0:
            print(f"    Warning: {len(down_genes_missing)}/{len(down_genes)} downregulated genes not found")

        if len(down_genes_in_data) > 0:
            print(f"    Negative: Using {len(down_genes_in_data)}/{len(down_genes)} downregulated genes")
            sc.tl.score_genes(adata, gene_list=down_genes_in_data,
                            score_name=f"{gene_set_name}_negative_score", use_raw=False)
        else:
            print(f"    Warning: No downregulated genes found in dataset")
            adata.obs[f"{gene_set_name}_negative_score"] = 0.0

        # Calculate aggregate score
        # Aggregate = positive contribution from upregulated - positive contribution from downregulated
        # (high expression of up genes is good, high expression of down genes is bad)
        if f"{gene_set_name}_positive_score" in adata.obs.columns and \
           f"{gene_set_name}_negative_score" in adata.obs.columns:
            adata.obs[f"{gene_set_name}_aggregate_score"] = \
                adata.obs[f"{gene_set_name}_positive_score"] - \
                adata.obs[f"{gene_set_name}_negative_score"]
            print(f"    Aggregate score calculated")

    return adata

def create_spatial_plots(adata, gene_sets_dict, output_dir):
    """
    Create spatial plots for directional gene set scores with consistent color scales per time point.

    Creates separate plots for positive, negative, and aggregate scores.

    Parameters:
    -----------
    adata : AnnData
        AnnData object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names
    output_dir : str
        Directory to save plots
    """
    print("Creating spatial plots...")

    # Ensure spatial coordinates exist
    if 'spatial' not in adata.obsm:
        adata.obsm['spatial'] = np.c_[adata.obs['center_x'], adata.obs['center_y']]

    # Create spatial metadata in .uns (required by squidpy)
    if 'spatial' not in adata.uns:
        adata.uns['spatial'] = {}
        for sample_id in adata.obs['sample_id'].unique():
            adata.uns['spatial'][sample_id] = {
                'images': {},
                'scalefactors': {'spot_diameter_fullres': 1.0, 'hires': 1.0}
            }

    samples = sorted(adata.obs['sample_id'].unique())

    # Define score types and their visualization properties
    score_configs = {
        'positive': {'cmap': 'Reds', 'label': 'Positive Score'},
        'negative': {'cmap': 'Blues', 'label': 'Negative Score'},
        'aggregate': {'cmap': 'coolwarm', 'label': 'Aggregate Score'}
    }

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
        for score_type, config in score_configs.items():
            score_column = f"{gene_set_name}_{score_type}_score"

            # Skip if score column doesn't exist
            if score_column not in adata.obs.columns:
                print(f"  Skipping {score_type} - no scores calculated")
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

                # Calculate color scale based on quantiles
                valid_scores = scores.dropna()
                if len(valid_scores) == 0:
                    plt.close(fig)
                    continue

                # Create scatter plot
                scatter = ax.scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                                   c=scores, cmap=config['cmap'], s=5, alpha=0.8,
                                   vmin=vmin, vmax=vmax)

                # Add colorbar
                cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
                cbar.set_label(f'{gene_set_name}\n{config["label"]}',
                             rotation=270, labelpad=20)

                # Set labels and title
                ax.set_xlabel('X coordinate')
                ax.set_ylabel('Y coordinate')
                ax.set_title(f'{gene_set_name} - {score_type.capitalize()} Score\n'
                           f'Sample {sample_id}\n'
                           f'Scale: {vmin:.3f} to {vmax:.3f} ({time_point} range)')
                ax.set_aspect('equal')

                # Save the plot
                plot_filename = f"{clean_name}_{score_type}_sample_{sample_id}.png"
                fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight', dpi=300)
                plt.close(fig)

            print(f"  Saved {score_type} score plots")

def create_summary_statistics(adata, gene_sets_dict, output_dir):
    """
    Create summary statistics for directional gene set scores.

    Parameters:
    -----------
    adata : AnnData
        AnnData object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set information
    output_dir : str
        Directory to save summary
    """
    print("\nCreating summary statistics...")

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
                'N_Genes_Up_Found': len([g for g in gene_sets_dict[gene_set_name]['up']
                                        if g in adata.var.index]),
                'N_Genes_Down_Found': len([g for g in gene_sets_dict[gene_set_name]['down']
                                          if g in adata.var.index])
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
                    'N_Genes_Up_Found': len([g for g in gene_sets_dict[gene_set_name]['up']
                                            if g in adata.var.index]),
                    'N_Genes_Down_Found': len([g for g in gene_sets_dict[gene_set_name]['down']
                                              if g in adata.var.index])
                }
                summary_data.append(sample_stats)

        # Convert to DataFrame and save
        summary_df = pd.DataFrame(summary_data)
        summary_filename = f'gene_set_{score_type}_scores_summary.csv'
        summary_df.to_csv(os.path.join(output_dir, summary_filename), index=False)
        print(f"  Saved: {summary_filename}")

    # Create a formatted text report
    print("\nCreating summary report...")
    with open(os.path.join(output_dir, 'gene_set_scores_report.txt'), 'w') as f:
        f.write("Directional Gene Set Scoring Summary Report\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"Total gene sets analyzed: {len(gene_sets_dict)}\n")
        f.write(f"Total cells: {len(adata)}\n")
        f.write(f"Total samples: {len(adata.obs['sample_id'].unique())}\n\n")

        f.write("Gene Sets Overview:\n")
        f.write("-" * 30 + "\n")
        for gene_set_name, gene_info in gene_sets_dict.items():
            n_up = len(gene_info['up'])
            n_down = len(gene_info['down'])
            n_up_found = len([g for g in gene_info['up'] if g in adata.var.index])
            n_down_found = len([g for g in gene_info['down'] if g in adata.var.index])
            f.write(f"\n{gene_set_name}:\n")
            f.write(f"  Total genes: {n_up + n_down}\n")
            f.write(f"  Upregulated: {n_up_found}/{n_up} genes found\n")
            f.write(f"  Downregulated: {n_down_found}/{n_down} genes found\n")

        f.write(f"\nScore Types Generated:\n")
        f.write("-" * 30 + "\n")
        f.write("  - Positive scores: Based on upregulated genes only\n")
        f.write("  - Negative scores: Based on downregulated genes only\n")
        f.write("  - Aggregate scores: Combined score (positive - negative)\n")

        f.write(f"\nOutput Files:\n")
        f.write("-" * 30 + "\n")
        f.write(f"  - Spatial plots: {output_dir}/\n")
        f.write(f"  - Summary statistics: gene_set_*_scores_summary.csv\n")

    print(f"  Saved: gene_set_scores_report.txt")

def main():
    """Main function to run the directional gene set scoring and spatial visualization pipeline."""

    print("Directional Gene Set Scoring and Spatial Visualization Pipeline")
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
    adata = calculate_directional_gene_set_scores(adata, gene_sets_dict)

    # Create spatial plots
    create_spatial_plots(adata, gene_sets_dict, output_dir)

    # Create summary statistics
    create_summary_statistics(adata, gene_sets_dict, output_dir)

    print("\n" + "=" * 80)
    print("DIRECTIONAL GENE SET SCORING AND SPATIAL VISUALIZATION COMPLETED")
    print("=" * 80)
    print(f"Results saved to: {output_dir}")
    print("Generated outputs:")
    print("  - Spatial plots for each gene set and sample (positive, negative, aggregate)")
    print("  - Summary statistics: gene_set_*_scores_summary.csv")
    print("  - Summary report: gene_set_scores_report.txt")
    print("\nScore Types:")
    print("  - Positive: Based on upregulated genes (colormap: Reds)")
    print("  - Negative: Based on downregulated genes (colormap: Blues)")
    print("  - Aggregate: Combined score = positive - negative (colormap: coolwarm)")
    print("=" * 80)

if __name__ == "__main__":
    main()
