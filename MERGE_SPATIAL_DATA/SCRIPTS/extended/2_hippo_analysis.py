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
merged_h5ad_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_extended.h5ad'

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
output_dir = 'MERGE_SPATIAL_DATA/plots_hippo_cluster_highlight_extended'
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

# %%
# Analyze cluster sizes per sample and identify clusters that would be extended
print("Analyzing cluster sizes per sample...")
extension_summary = {}

for sample_id in merged_adata.obs['sample_id'].unique():
    sample_adata = merged_adata[merged_adata.obs['sample_id'] == sample_id]
    cluster_counts = sample_adata.obs['leiden'].value_counts()
    small_clusters = cluster_counts[(cluster_counts > 0) & (cluster_counts < 50)].index.tolist()
    large_clusters = cluster_counts[cluster_counts >= 50].index.tolist()
    
    extension_summary[sample_id] = {
        'small_clusters': small_clusters,
        'large_clusters': large_clusters,
        'total_clusters': len(cluster_counts)
    }
    
    print(f"Sample {sample_id}:")
    print(f"  - Clusters with <50 cells (would be extended): {small_clusters}")
    print(f"  - Clusters with ≥50 cells: {large_clusters}")


# %%
# Plot leiden clusters spatially, faceted by sample
sq.pl.spatial_scatter(merged_adata, color='leiden', shape=None,
                      library_key='sample_id',
                      save='spatial_leiden_by_sample.png')

# %%
# Plot leiden clusters spatially for each sample with extension information
for sample_id in merged_adata.obs['sample_id'].unique():
    sample_adata = merged_adata[merged_adata.obs['sample_id'] == sample_id].copy()
    small_clusters = extension_summary[sample_id]['small_clusters']
    
    title_text = f"Spatial Leiden Clusters - Sample: {sample_id}"
    if small_clusters:
        title_text += f"\n(Clusters {small_clusters} would be extended in analysis)"
    
    print(f"Plotting spatial scatter for sample: {sample_id}")
    sq.pl.spatial_scatter(sample_adata, color='leiden', shape=None,
                          title=title_text,
                          save=f'spatial_leiden_by_sample_{sample_id}.png')

# %%
# Get unique samples and clusters
samples = merged_adata.obs['sample_id'].unique()
clusters = merged_adata.obs['leiden'].unique()

# %%
# Generate and save plots for each sample and cluster
for sample in samples:
    sample_adata = merged_adata[merged_adata.obs['sample_id'] == sample]
    sample_clusters = sample_adata.obs['leiden'].unique()
    small_clusters = extension_summary[sample]['small_clusters']
    
    for cluster in sample_clusters:
        # Check if cluster has more than 50 cells
        cluster_mask = sample_adata.obs['leiden'] == cluster
        cluster_cell_count = cluster_mask.sum()
        
        # Determine if this cluster would be extended
        is_extended_cluster = str(cluster) in [str(c) for c in small_clusters]
        extension_text = " (EXTENDED)" if is_extended_cluster else ""
        
        print(f"Plotting sample {sample}, cluster {cluster} ({cluster_cell_count} cells){extension_text}")

        # Create a figure with two subplots
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))
        
        # Get the color for this cluster from leiden colors
        if 'leiden_colors' in merged_adata.uns:
            cluster_colors = dict(zip(merged_adata.obs['leiden'].cat.categories, merged_adata.uns['leiden_colors']))
            cluster_color = cluster_colors.get(cluster, 'blue')
        else:
            cluster_color = 'blue'
        
        # UMAP plot - show extension concept
        umap_coords = merged_adata.obsm['X_umap']
        sample_mask = merged_adata.obs['sample_id'] == sample
        full_cluster_mask = (merged_adata.obs['leiden'] == cluster) & sample_mask
        
        # For extended clusters, also show all cells from this cluster across all samples
        if is_extended_cluster:
            all_cluster_mask = merged_adata.obs['leiden'] == cluster
            
            # Plot all cells in gray
            axes[0].scatter(umap_coords[:, 0], umap_coords[:, 1], c='lightgray', s=1, alpha=0.3)
            
            # Plot all cells from this cluster across all samples (lighter)
            axes[0].scatter(umap_coords[all_cluster_mask, 0], umap_coords[all_cluster_mask, 1],
                          c=cluster_color, s=3, alpha=0.4, label=f'All Cluster {cluster} cells')
            
            # Plot cells from current sample (darker)
            if full_cluster_mask.sum() > 0:
                axes[0].scatter(umap_coords[full_cluster_mask, 0], umap_coords[full_cluster_mask, 1],
                              c=cluster_color, s=8, alpha=1.0, edgecolors='black', linewidths=0.5,
                              label=f'Sample {sample} cells')
            
            axes[0].legend(loc='upper right', fontsize=8)
            umap_title = f"UMAP: {sample} - Cluster {cluster} (EXTENDED)\nShowing all cluster cells vs. sample-specific"
            
        else:
            # Plot all cells in gray
            axes[0].scatter(umap_coords[:, 0], umap_coords[:, 1], c='lightgray', s=1, alpha=0.3)
            
            # Plot highlighted sample-cluster combination in color
            if full_cluster_mask.sum() > 0:
                axes[0].scatter(umap_coords[full_cluster_mask, 0], umap_coords[full_cluster_mask, 1],
                              c=cluster_color, s=5, alpha=0.8)
            
            umap_title = f"UMAP: {sample} - Cluster {cluster}"
        
        axes[0].set_xlabel('UMAP1')
        axes[0].set_ylabel('UMAP2')
        axes[0].set_title(umap_title, fontsize=10)
        
        # Spatial plot using matplotlib directly
        spatial_coords = sample_adata.obsm['spatial']
        
        # Plot all cells in gray
        axes[1].scatter(spatial_coords[:, 0], spatial_coords[:, 1], c='lightgray', s=1, alpha=0.5)
        
        # Plot highlighted cluster in the same color as UMAP
        axes[1].scatter(spatial_coords[cluster_mask, 0], spatial_coords[cluster_mask, 1],
                      c=cluster_color, s=3, alpha=1.0)
        
        axes[1].set_xlabel('X coordinate')
        axes[1].set_ylabel('Y coordinate')
        spatial_title = f"Spatial: {sample} - Cluster {cluster}"
        if is_extended_cluster:
            spatial_title += f"\n({cluster_cell_count} cells, extended with cells from other samples)"
        axes[1].set_title(spatial_title, fontsize=10)
        axes[1].set_aspect('equal')
        
        # Save the figure
        plot_filename = f"sample_{sample}_cluster_{cluster}.png"
        fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight')
        plt.close(fig)

# %%
# Create an extension summary report
print("\n" + "="*60)
print("CLUSTER EXTENSION SUMMARY")
print("="*60)

total_small_clusters = 0
total_large_clusters = 0

for sample_id, info in extension_summary.items():
    print(f"\nSample: {sample_id}")
    print(f"  Total clusters: {info['total_clusters']}")
    print(f"  Large clusters (≥50 cells): {len(info['large_clusters'])} - {info['large_clusters']}")
    print(f"  Small clusters (<50 cells, EXTENDED): {len(info['small_clusters'])} - {info['small_clusters']}")
    
    total_small_clusters += len(info['small_clusters'])
    total_large_clusters += len(info['large_clusters'])

print(f"\nOVERALL SUMMARY:")
print(f"  Total large clusters across all samples: {total_large_clusters}")
print(f"  Total small clusters that will be extended: {total_small_clusters}")

# Save extension summary to file
summary_file = os.path.join(output_dir, 'cluster_extension_summary.txt')
with open(summary_file, 'w') as f:
    f.write("CLUSTER EXTENSION SUMMARY\n")
    f.write("=" * 60 + "\n\n")
    f.write("This file shows which clusters would be 'extended' in marker gene analysis.\n")
    f.write("Extended clusters are those with <50 cells in a sample's ROI.\n")
    f.write("For these clusters, marker gene analysis includes ALL cells from the\n")
    f.write("same cluster across the entire dataset for better statistical power.\n\n")
    
    for sample_id, info in extension_summary.items():
        f.write(f"Sample: {sample_id}\n")
        f.write(f"  Total clusters: {info['total_clusters']}\n")
        f.write(f"  Large clusters (≥50 cells): {len(info['large_clusters'])} - {info['large_clusters']}\n")
        f.write(f"  Small clusters (<50 cells, EXTENDED): {len(info['small_clusters'])} - {info['small_clusters']}\n\n")
    
    f.write(f"OVERALL SUMMARY:\n")
    f.write(f"  Total large clusters across all samples: {total_large_clusters}\n")
    f.write(f"  Total small clusters that will be extended: {total_small_clusters}\n")

print(f"\nExtension summary saved to: {summary_file}")
print("Script finished successfully.")
# %%