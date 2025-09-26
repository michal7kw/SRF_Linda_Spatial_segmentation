# %%
import scanpy as sc
import squidpy as sq
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Define the path to the merged AnnData file
merged_h5ad_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'

# %%
# Load the merged AnnData object
if os.path.exists(merged_h5ad_path):
    print(f"Loading merged data from: {merged_h5ad_path}")
    merged_adata = sc.read_h5ad(merged_h5ad_path)
    print("Merged data loaded successfully.")
else:
    print(f"Error: Merged data file not found at {merged_h5ad_path}")
    sys.exit()

# %%
merged_adata

# %%
# Create spatial coordinates
merged_adata.obsm['spatial'] = np.c_[merged_adata.obs['center_x'], merged_adata.obs['center_y']]

# %%
# Create spatial metadata in .uns
merged_adata.uns['spatial'] = {}
for sample_id in merged_adata.obs['sample_id'].unique():
    merged_adata.uns['spatial'][sample_id] = {'images': {}, 'scalefactors': {'spot_diameter_fullres': 1.0, 'hires': 1.0}}

# %%
# Define output directory for plots
output_dir = 'MERGE_SPATIAL_DATA/plots_hippo_cluster_highlight_inclusive'
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

# %%
# Analyze cluster sizes per sample for reporting
print("Analyzing cluster sizes per sample...")
cluster_summary = {}

for sample_id in merged_adata.obs['sample_id'].unique():
    sample_data = merged_adata[merged_adata.obs['sample_id'] == sample_id]
    cluster_counts = sample_data.obs['leiden'].value_counts().sort_index()
    
    small_clusters = cluster_counts[cluster_counts < 50].index.tolist()
    large_clusters = cluster_counts[cluster_counts >= 50].index.tolist()
    
    cluster_summary[sample_id] = {
        'small_clusters': small_clusters,
        'large_clusters': large_clusters,
        'total_clusters': len(cluster_counts)
    }
    
    print(f"\nSample {sample_id}:")
    print(f"  Total clusters: {len(cluster_counts)}")
    print(f"  Large clusters (≥50 cells): {len(large_clusters)} - {large_clusters}")
    print(f"  Small clusters (<50 cells): {len(small_clusters)} - {small_clusters} (will be included in analysis)")

# %%
# Plot leiden clusters spatially, faceted by sample
sq.pl.spatial_scatter(merged_adata, color='leiden', shape=None,
                      library_key='sample_id',
                      save='spatial_leiden_by_sample.png')

# %%
# Plot leiden clusters spatially for each sample
for sample_id in merged_adata.obs['sample_id'].unique():
    sample_adata = merged_adata[merged_adata.obs['sample_id'] == sample_id].copy()
    
    title_text = f"Spatial Leiden Clusters - Sample: {sample_id}"
    
    print(f"Plotting spatial scatter for sample: {sample_id}")
    sq.pl.spatial_scatter(sample_adata, color='leiden', shape=None,
                          title=title_text,
                          save=f'spatial_leiden_by_sample_{sample_id}.png')

# %%
# Get unique samples and clusters
samples = merged_adata.obs['sample_id'].unique()
clusters = merged_adata.obs['leiden'].unique()

# %%
# Generate and save plots for each sample and cluster (no filtering)
for sample in samples:
    sample_adata = merged_adata[merged_adata.obs['sample_id'] == sample]
    
    for cluster in sorted(sample_adata.obs['leiden'].unique()):
        cluster_mask = sample_adata.obs['leiden'] == cluster
        cluster_cell_count = cluster_mask.sum()
        
        print(f"Plotting sample {sample}, cluster {cluster} ({cluster_cell_count} cells)")

        # Create a figure with two subplots
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))
        
        # Get the color for this cluster from leiden colors
        if 'leiden_colors' in merged_adata.uns:
            cluster_colors = dict(zip(merged_adata.obs['leiden'].cat.categories, merged_adata.uns['leiden_colors']))
            cluster_color = cluster_colors.get(cluster, 'blue')
        else:
            cluster_color = 'blue'
        
        # UMAP plot
        umap_coords = merged_adata.obsm['X_umap']
        sample_mask = merged_adata.obs['sample_id'] == sample
        full_cluster_mask = (merged_adata.obs['leiden'] == cluster) & sample_mask
        
        # Plot all cells in gray on the first subplot
        axes[0].scatter(umap_coords[:, 0], umap_coords[:, 1], c='lightgray', s=1, alpha=0.3)
        
        # Plot highlighted sample-cluster combination in color on the first subplot
        if full_cluster_mask.sum() > 0:
            axes[0].scatter(umap_coords[full_cluster_mask, 0], umap_coords[full_cluster_mask, 1],
                          c=cluster_color, s=5, alpha=0.8)
        
        axes[0].set_xlabel('UMAP1')
        axes[0].set_ylabel('UMAP2')
        axes[0].set_title(f"UMAP: {sample} - Cluster {cluster}", fontsize=10)
        
        # Spatial plot using matplotlib directly on the second subplot
        spatial_coords = sample_adata.obsm['spatial']
        
        # Plot all cells in gray on the second subplot
        axes[1].scatter(spatial_coords[:, 0], spatial_coords[:, 1], c='lightgray', s=1, alpha=0.5)
        
        # Plot highlighted cluster in the same color as UMAP on the second subplot
        axes[1].scatter(spatial_coords[cluster_mask, 0], spatial_coords[cluster_mask, 1],
                      c=cluster_color, s=3, alpha=1.0)
        
        axes[1].set_xlabel('X coordinate')
        axes[1].set_ylabel('Y coordinate')
        axes[1].set_title(f"Spatial: {sample} - Cluster {cluster}\n({cluster_cell_count} cells)", fontsize=10)
        axes[1].set_aspect('equal')
        
        # Save the figure
        plot_filename = f"sample_{sample}_cluster_{cluster}.png"
        fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight')
        plt.close(fig)

# %%
# Create an inclusive analysis summary report
print("\n" + "="*60)
print("INCLUSIVE ANALYSIS SUMMARY (All Clusters Included)")
print("="*60)

total_small_clusters = 0
total_large_clusters = 0
total_clusters_analyzed = 0

for sample_id, info in cluster_summary.items():
    print(f"\nSample: {sample_id}")
    print(f"  Total clusters found: {info['total_clusters']}")
    print(f"  Large clusters (≥50 cells): {len(info['large_clusters'])} - {info['large_clusters']}")
    print(f"  Small clusters (<50 cells): {len(info['small_clusters'])} - {info['small_clusters']}")
    
    total_small_clusters += len(info['small_clusters'])
    total_large_clusters += len(info['large_clusters'])
    total_clusters_analyzed += info['total_clusters']

print(f"\nOVERALL SUMMARY:")
print(f"  Total clusters analyzed across all samples: {total_clusters_analyzed}")
print(f"  Total large clusters: {total_large_clusters}")
print(f"  Total small clusters included in analysis: {total_small_clusters}")

# Save inclusive analysis summary to file
summary_file = os.path.join(output_dir, 'inclusive_analysis_summary.txt')
with open(summary_file, 'w') as f:
    f.write("INCLUSIVE ANALYSIS SUMMARY\n")
    f.write("=" * 60 + "\n\n")
    f.write("This file shows results from inclusive analysis where all clusters\n")
    f.write("are analyzed regardless of cell count.\n\n")
    
    for sample_id, info in cluster_summary.items():
        f.write(f"Sample: {sample_id}\n")
        f.write(f"  Total clusters found: {info['total_clusters']}\n")
        f.write(f"  Large clusters analyzed (≥50 cells): {len(info['large_clusters'])} - {info['large_clusters']}\n")
        f.write(f"  Small clusters analyzed (<50 cells): {len(info['small_clusters'])} - {info['small_clusters']}\n\n")
    
    f.write(f"OVERALL SUMMARY:\n")
    f.write(f"  Total clusters analyzed: {total_clusters_analyzed}\n")
    f.write(f"  Total large clusters: {total_large_clusters}\n")
    f.write(f"  Total small clusters included: {total_small_clusters}\n")

print(f"\nInclusive analysis summary saved to: {summary_file}")
print("Script finished successfully.")
# %%