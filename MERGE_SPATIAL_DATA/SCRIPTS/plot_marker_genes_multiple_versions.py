#!/usr/bin/env python3
"""
Script to generate comprehensive marker gene visualizations with different numbers of markers
Creates separate folders for different marker counts to avoid overcrowding
"""

import os
import sys
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import warnings
warnings.filterwarnings('ignore')

# Set project root
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Configuration
MERGED_H5AD_PATH = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
MARKER_GENES_DIR = 'MERGE_SPATIAL_DATA/3_MARKER_GENES/marker_genes_per_sample_inclusive_new_clusters'
OUTPUT_BASE_DIR = 'MERGE_SPATIAL_DATA/8_ADDITIONAL_PLOTS'
CLUSTER_COLUMN = 'leiden_new'

# Different marker counts to generate
MARKER_CONFIGS = [
    {'top_n': 3, 'folder': 'top3_markers'},
    {'top_n': 5, 'folder': 'top5_markers'},
    {'top_n': 10, 'folder': 'top10_markers'}
]

# Sample list
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

# Set scanpy figure parameters
sc.set_figure_params(dpi=100, fontsize=14)


def load_marker_genes(sample_name, marker_genes_dir, top_n=5):
    """Load top marker genes for a specific sample from CSV file"""
    folder_name = f"sample_{sample_name}"
    csv_path = os.path.join(marker_genes_dir, folder_name, 'top_30_marker_genes.csv')

    if not os.path.exists(csv_path):
        print(f"Warning: Could not find {csv_path}")
        return None

    # Read CSV file without headers
    df = pd.read_csv(csv_path, header=None)

    # First row contains cluster IDs
    cluster_ids = df.iloc[0].values

    # Collect all marker genes (top N per cluster)
    all_marker_genes = []
    for col_idx in range(len(cluster_ids)):
        if len(df) > 1:
            genes = df.iloc[1:top_n+1, col_idx].values
            genes = [g for g in genes if pd.notna(g)]
            all_marker_genes.extend(genes)

    # Remove duplicates while preserving order
    unique_marker_genes = []
    seen = set()
    for gene in all_marker_genes:
        if gene not in seen:
            unique_marker_genes.append(gene)
            seen.add(gene)

    return cluster_ids, unique_marker_genes


def create_dotplot(sample_adata, available_genes, sample_name, output_dir, cluster_column='leiden_new'):
    """Create dotplot visualization"""
    try:
        print("    Creating dotplot...")

        n_clusters = len(sample_adata.obs[cluster_column].unique())
        n_genes = len(available_genes)

        # Adjust figure size based on number of genes/clusters
        fig_width = max(10, n_clusters * 0.6 + 3)
        fig_height = max(6, n_genes * 0.35 + 2)

        sc.pl.dotplot(
            sample_adata,
            var_names=available_genes,
            groupby=cluster_column,
            dendrogram=True,
            standard_scale='var',
            cmap='RdBu_r',
            figsize=(fig_width, fig_height),
            show=False,
            smallest_dot=50,
            dot_max=0.7,
            swap_axes=True
        )

        output_path = os.path.join(output_dir, f'{sample_name}_dotplot.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"      Saved: {os.path.basename(output_path)}")
        return True
    except Exception as e:
        print(f"      Error: {e}")
        plt.close()
        return False


def create_matrixplot(sample_adata, available_genes, sample_name, output_dir, cluster_column='leiden_new'):
    """Create matrixplot (heatmap of mean expression)"""
    try:
        print("    Creating matrixplot...")

        n_clusters = len(sample_adata.obs[cluster_column].unique())
        n_genes = len(available_genes)

        fig_width = max(10, n_clusters * 0.6 + 2)
        fig_height = max(6, n_genes * 0.35 + 2)

        sc.pl.matrixplot(
            sample_adata,
            var_names=available_genes,
            groupby=cluster_column,
            dendrogram=True,
            standard_scale='var',
            cmap='RdBu_r',
            figsize=(fig_width, fig_height),
            show=False,
            swap_axes=True
        )

        output_path = os.path.join(output_dir, f'{sample_name}_matrixplot.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"      Saved: {os.path.basename(output_path)}")
        return True
    except Exception as e:
        print(f"      Error: {e}")
        plt.close()
        return False


def create_stacked_violin(sample_adata, available_genes, sample_name, output_dir, cluster_column='leiden_new'):
    """Create stacked violin plot"""
    try:
        print("    Creating stacked violin...")

        n_clusters = len(sample_adata.obs[cluster_column].unique())
        n_genes = len(available_genes)

        fig_width = max(10, n_genes * 0.5 + 2)
        fig_height = max(6, n_clusters * 0.45 + 2)

        sc.pl.stacked_violin(
            sample_adata,
            var_names=available_genes,
            groupby=cluster_column,
            dendrogram=True,
            figsize=(fig_width, fig_height),
            show=False,
            swap_axes=True
        )

        output_path = os.path.join(output_dir, f'{sample_name}_stacked_violin.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"      Saved: {os.path.basename(output_path)}")
        return True
    except Exception as e:
        print(f"      Error: {e}")
        plt.close()
        return False


def create_heatmap(sample_adata, available_genes, sample_name, output_dir, cluster_column='leiden_new'):
    """Create heatmap (shows individual cells)"""
    try:
        print("    Creating heatmap...")

        n_genes = len(available_genes)
        fig_height = max(8, n_genes * 0.35 + 2)

        sc.pl.heatmap(
            sample_adata,
            var_names=available_genes,
            groupby=cluster_column,
            figsize=(12, fig_height),
            show=False,
            swap_axes=False,
            standard_scale='var',
            cmap='RdBu_r'
        )

        output_path = os.path.join(output_dir, f'{sample_name}_heatmap.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"      Saved: {os.path.basename(output_path)}")
        return True
    except Exception as e:
        print(f"      Error: {e}")
        plt.close()
        return False


def create_tracksplot(sample_adata, available_genes, sample_name, output_dir, cluster_column='leiden_new'):
    """Create tracksplot (expression as height)"""
    try:
        print("    Creating tracksplot...")

        n_genes = len(available_genes)
        fig_height = max(8, n_genes * 0.4 + 2)

        sc.pl.tracksplot(
            sample_adata,
            var_names=available_genes,
            groupby=cluster_column,
            figsize=(12, fig_height),
            show=False
        )

        output_path = os.path.join(output_dir, f'{sample_name}_tracksplot.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"      Saved: {os.path.basename(output_path)}")
        return True
    except Exception as e:
        print(f"      Error: {e}")
        plt.close()
        return False


def main():
    """Main execution function"""

    print("="*80)
    print("MARKER GENE VISUALIZATION - MULTIPLE VERSIONS")
    print("="*80)
    print()

    # Load merged spatial data
    print(f"Loading merged data from: {MERGED_H5AD_PATH}")
    adata = sc.read_h5ad(MERGED_H5AD_PATH)

    print(f"Loaded data shape: {adata.shape}")
    print(f"Samples: {sorted(adata.obs['sample_id'].unique())}")
    print(f"Clusters: {sorted(adata.obs[CLUSTER_COLUMN].unique())}\n")

    # Process each marker configuration
    for config in MARKER_CONFIGS:
        top_n = config['top_n']
        folder_name = config['folder']

        print("="*80)
        print(f"GENERATING PLOTS WITH TOP {top_n} MARKERS PER CLUSTER")
        print("="*80)
        print()

        # Create output directory for this configuration
        output_dir = os.path.join(OUTPUT_BASE_DIR, folder_name)
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory: {output_dir}\n")

        # Process each sample
        for sample_name in SAMPLES:
            print(f"  Processing: {sample_name}")

            # Load marker genes
            result = load_marker_genes(sample_name, MARKER_GENES_DIR, top_n=top_n)

            if result is None:
                continue

            cluster_ids, marker_genes = result
            print(f"    Loaded {len(marker_genes)} unique marker genes")

            # Subset data for this sample
            sample_adata = adata[adata.obs['sample_id'] == sample_name].copy()
            print(f"    Sample cells: {sample_adata.n_obs}")

            # Filter genes that exist in the dataset
            available_genes = [g for g in marker_genes if g in sample_adata.var_names]
            print(f"    Available genes: {len(available_genes)}")

            if len(available_genes) == 0:
                print("    Warning: No marker genes found")
                continue

            # Create all visualization types
            create_dotplot(sample_adata, available_genes, sample_name, output_dir, CLUSTER_COLUMN)
            create_matrixplot(sample_adata, available_genes, sample_name, output_dir, CLUSTER_COLUMN)
            create_stacked_violin(sample_adata, available_genes, sample_name, output_dir, CLUSTER_COLUMN)
            create_heatmap(sample_adata, available_genes, sample_name, output_dir, CLUSTER_COLUMN)
            create_tracksplot(sample_adata, available_genes, sample_name, output_dir, CLUSTER_COLUMN)

            print()

        print()

    print("="*80)
    print("COMPLETED")
    print("="*80)
    print(f"\nGenerated plots in {OUTPUT_BASE_DIR}/")
    print("Folders created:")
    for config in MARKER_CONFIGS:
        print(f"  - {config['folder']} (top {config['top_n']} markers per cluster)")
    print()


if __name__ == "__main__":
    main()
