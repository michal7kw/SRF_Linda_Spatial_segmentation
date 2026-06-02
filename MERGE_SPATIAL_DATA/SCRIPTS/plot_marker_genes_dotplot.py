#!/usr/bin/env python3
"""
Script to generate marker gene expression dotplots for each sample
Creates plots similar to the reference image showing marker genes across cell types
"""

import os
import sys
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Configuration
MERGED_H5AD_PATH = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
MARKER_GENES_DIR = 'MERGE_SPATIAL_DATA/3_MARKER_GENES/marker_genes_per_sample_inclusive_new_clusters'
OUTPUT_DIR = 'MERGE_SPATIAL_DATA/8_ADDITIONAL_PLOTS'
CLUSTER_COLUMN = 'leiden_new'
TOP_N_GENES = 5  # Number of top marker genes per cluster to display

# Sample list (as they appear in the h5ad file, without "sample_" prefix)
SAMPLES = [
    'e16_R1_mut',
    'e16_R2_ctrl',
    'p0_R5p_mut',
    'p0_R6_ctrl',
    'p7_R1_mut',
    'p7_R2_ctrl',
    'p30_R3_mut',
    'p30_R4p_ctrl'
]

def load_marker_genes(sample_name, marker_genes_dir, top_n=5):
    """
    Load top marker genes for a specific sample from CSV file

    Args:
        sample_name: Sample identifier (without "sample_" prefix)
        marker_genes_dir: Directory containing marker gene files
        top_n: Number of top genes per cluster to extract

    Returns:
        Dictionary mapping cluster to list of top marker genes
    """
    # Add "sample_" prefix for folder lookup
    folder_name = f"sample_{sample_name}"
    csv_path = os.path.join(marker_genes_dir, folder_name, 'top_30_marker_genes.csv')

    if not os.path.exists(csv_path):
        print(f"Warning: Could not find {csv_path}")
        return None

    # Read CSV file
    df = pd.read_csv(csv_path, header=None)

    # First row contains cluster IDs
    cluster_ids = df.iloc[0].values

    # Subsequent rows contain marker genes
    marker_genes_by_cluster = {}

    for col_idx, cluster_id in enumerate(cluster_ids):
        # Get top N genes for this cluster (skip first row which is cluster ID)
        genes = df.iloc[1:top_n+1, col_idx].values
        # Remove any NaN values and convert to list
        genes = [g for g in genes if pd.notna(g)]
        marker_genes_by_cluster[str(cluster_id)] = genes

    return cluster_ids, marker_genes_by_cluster


def create_dotplot_for_sample(adata, sample_name, marker_genes_dict, cluster_ids, output_dir, cluster_column='leiden_new'):
    """
    Create dotplot for a specific sample showing marker gene expression across clusters

    Args:
        adata: AnnData object containing the merged spatial data
        sample_name: Name of the sample
        marker_genes_dict: Dictionary mapping cluster to marker genes
        cluster_ids: Array of cluster IDs in order
        output_dir: Directory to save output plots
        cluster_column: Column name containing cluster assignments
    """
    # Filter data for this sample
    sample_adata = adata[adata.obs['sample_id'] == sample_name].copy()

    if sample_adata.n_obs == 0:
        print(f"Warning: No cells found for sample {sample_name}")
        return

    print(f"\nProcessing {sample_name}:")
    print(f"  Cells: {sample_adata.n_obs}")
    print(f"  Clusters: {sorted(sample_adata.obs[cluster_column].unique())}")

    # Collect all marker genes in order
    all_marker_genes = []
    for cluster_id in cluster_ids:
        cluster_id_str = str(cluster_id)
        if cluster_id_str in marker_genes_dict:
            genes = marker_genes_dict[cluster_id_str]
            all_marker_genes.extend(genes)

    # Remove duplicates while preserving order
    unique_marker_genes = []
    seen = set()
    for gene in all_marker_genes:
        if gene not in seen:
            unique_marker_genes.append(gene)
            seen.add(gene)

    # Filter genes that exist in the dataset
    available_genes = [g for g in unique_marker_genes if g in sample_adata.var_names]

    print(f"  Total marker genes: {len(unique_marker_genes)}")
    print(f"  Available in data: {len(available_genes)}")

    if len(available_genes) == 0:
        print(f"Warning: No marker genes found in dataset for {sample_name}")
        return

    # Create dotplot with proper sizing and spacing
    n_clusters = len(sample_adata.obs[cluster_column].unique())
    n_genes = len(available_genes)

    # Calculate figure size: tight spacing for dots, extra width for legend
    # Each cluster needs ~0.4 inches, each gene needs ~0.35 inches
    plot_width = n_clusters * 0.4
    plot_height = n_genes * 0.35
    legend_width = 6  # Reserve 6 inches for legend

    fig_width = max(25, plot_width) + legend_width
    fig_height = max(18, plot_height)

    try:
        # Set large font sizes globally
        plt.rcParams.update({
            'font.size': 20,
            'axes.labelsize': 22,
            'axes.titlesize': 24,
            'xtick.labelsize': 20,
            'ytick.labelsize': 20,
            'legend.fontsize': 20,
            'legend.title_fontsize': 22,
            'figure.autolayout': False
        })

        # Create the dotplot with swap_axes=True (genes on y-axis, clusters on x-axis)
        # Dendrogram will be at the top (same side as cluster labels on x-axis)
        sc.pl.dotplot(
            sample_adata,
            var_names=available_genes,
            groupby=cluster_column,
            dendrogram=True,
            standard_scale='var',
            cmap='RdBu_r',
            figsize=(fig_width, fig_height),
            show=False,
            swap_axes=True,  # Clusters on x-axis, genes on y-axis
            dot_max=1.0,  # Maximum dot size - increased for bigger dots
            dot_min=0.0,
            size_title='Fraction of cells\nin group (%)',
            vmin=-2,
            vmax=2,
            var_group_rotation=0  # Keep gene labels horizontal if they were grouped
        )

        # Get the current figure
        fig = plt.gcf()
        fig.set_size_inches(fig_width, fig_height)

        # Customize all axes for better readability
        for ax in fig.axes:
            # Increase tick label sizes significantly
            ax.tick_params(axis='both', which='major', labelsize=20, length=8, width=2.5)

            # X-axis labels (cluster numbers) - rotated 90 degrees
            for label in ax.get_xticklabels():
                label.set_fontsize(20)
                label.set_rotation(90)
                label.set_ha('center')
                label.set_va('top')
                label.set_fontweight('bold')

            # Y-axis labels (gene names) - horizontal
            for label in ax.get_yticklabels():
                label.set_fontsize(20)
                label.set_ha('right')
                label.set_va('center')

            # Axis labels
            if ax.get_xlabel():
                ax.set_xlabel(ax.get_xlabel(), fontsize=24, fontweight='bold', labelpad=15)
            if ax.get_ylabel():
                ax.set_ylabel(ax.get_ylabel(), fontsize=24, fontweight='bold', labelpad=15)

            # Title
            if ax.get_title():
                ax.set_title(ax.get_title(), fontsize=26, fontweight='bold', pad=20)

        # Adjust layout: more space on right for legend, bottom for x-labels
        # Calculate the proportion for the plot area vs legend area
        plot_proportion = plot_width / fig_width
        plt.subplots_adjust(left=0.15, right=plot_proportion + 0.02, top=0.92, bottom=0.15)

        # Save with high DPI
        output_path = os.path.join(output_dir, f'{sample_name}_marker_genes_dotplot.png')
        plt.savefig(output_path, dpi=300, bbox_inches=None)  # Don't use tight to keep our spacing
        plt.close()

        print(f"  Saved: {output_path}")

    except Exception as e:
        print(f"Error creating dotplot for {sample_name}: {e}")
        plt.close()


def main():
    """Main execution function"""

    print("="*80)
    print("MARKER GENE DOTPLOT GENERATION")
    print("="*80)

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"\nOutput directory: {OUTPUT_DIR}")

    # Load merged spatial data
    print(f"\nLoading merged data from: {MERGED_H5AD_PATH}")
    adata = sc.read_h5ad(MERGED_H5AD_PATH)

    print(f"Loaded data shape: {adata.shape}")
    print(f"Samples in data: {adata.obs['sample_id'].unique()}")
    print(f"Clusters in data: {sorted(adata.obs[CLUSTER_COLUMN].unique())}")

    # Process each sample
    for sample_name in SAMPLES:
        print(f"\n{'='*80}")
        print(f"Processing: {sample_name}")
        print(f"{'='*80}")

        # Load marker genes for this sample
        result = load_marker_genes(sample_name, MARKER_GENES_DIR, top_n=TOP_N_GENES)

        if result is None:
            continue

        cluster_ids, marker_genes_dict = result

        print(f"Loaded marker genes for {len(marker_genes_dict)} clusters")

        # Create dotplot for this sample
        create_dotplot_for_sample(
            adata=adata,
            sample_name=sample_name,
            marker_genes_dict=marker_genes_dict,
            cluster_ids=cluster_ids,
            output_dir=OUTPUT_DIR,
            cluster_column=CLUSTER_COLUMN
        )

    print("\n" + "="*80)
    print("COMPLETED")
    print("="*80)
    print(f"All plots saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
