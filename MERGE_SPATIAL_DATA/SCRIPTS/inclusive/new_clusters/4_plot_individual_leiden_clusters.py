# %%
import os
import sys
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# %%
# Load the merged spatial data
merged_h5ad_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
print(f"Loading merged data from: {merged_h5ad_path}")
merged_adata = sc.read_h5ad(merged_h5ad_path)

print(f"Loaded data shape: {merged_adata.shape}")
print(f"Available observations: {list(merged_adata.obs.columns)}")

# %%
# Get leiden_new clusters
unique_clusters = sorted(merged_adata.obs['leiden_new'].unique())
print(f"Found {len(unique_clusters)} leiden_new clusters: {unique_clusters}")

# %%
# Create output directory
output_dir = 'MERGE_SPATIAL_DATA/plots_individual_leiden_new_clusters'
os.makedirs(output_dir, exist_ok=True)
print(f"Saving plots to {output_dir}")

# %%
# Calculate grid dimensions for subplots
n_clusters = len(unique_clusters)
cols = min(4, n_clusters)  # Maximum 4 columns
rows = (n_clusters + cols - 1) // cols  # Ceiling division

print(f"Creating {rows}x{cols} grid for {n_clusters} clusters")

# %%
# Create figure with subplots for individual cluster highlighting
fig = plt.figure(figsize=(5*cols, 4*rows))
gs = GridSpec(rows, cols, figure=fig, hspace=0.3, wspace=0.3)

# Plot each cluster
for i, cluster in enumerate(unique_clusters):
    row = i // cols
    col = i % cols
    
    # Create subplot
    ax = fig.add_subplot(gs[row, col])
    
    # Create highlighting: current cluster vs others
    highlight_col = f'highlight_{cluster}'
    merged_adata.obs[highlight_col] = 'Other'
    merged_adata.obs.loc[merged_adata.obs['leiden_new'] == cluster, highlight_col] = f'Cluster {cluster}'
    
    # Plot UMAP with current cluster highlighted
    sc.pl.umap(merged_adata, 
               color=highlight_col,
               ax=ax,
               show=False,
               frameon=False,
               title=f'Cluster {cluster}',
               legend_loc='none')
    
    # Customize the plot
    ax.set_title(f'Leiden_new Cluster {cluster}', fontsize=12, fontweight='bold')
    ax.set_xlabel('UMAP 1', fontsize=10)
    ax.set_ylabel('UMAP 2', fontsize=10)

# Remove empty subplots if any
for i in range(n_clusters, rows * cols):
    row = i // cols
    col = i % cols
    fig.delaxes(fig.add_subplot(gs[row, col]))

# Add main title
fig.suptitle('Individual Leiden_new Clusters - UMAP Visualization', fontsize=16, fontweight='bold', y=0.98)

# Save the figure
output_path = os.path.join(output_dir, 'individual_leiden_new_clusters_grid.png')
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Saved individual clusters grid plot: {output_path}")
plt.show()

# %%
# Create individual plots for each cluster (separate files)
print("Creating individual cluster plots...")

for cluster in unique_clusters:
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    
    # Create highlighting for individual plots
    highlight_col = f'highlight_individual_{cluster}'
    merged_adata.obs[highlight_col] = 'Other'
    merged_adata.obs.loc[merged_adata.obs['leiden_new'] == cluster, highlight_col] = f'Cluster {cluster}'
    
    # Plot UMAP
    sc.pl.umap(merged_adata, 
               color=highlight_col,
               ax=ax,
               show=False,
               frameon=True,
               title=f'Leiden_new Cluster {cluster}',
               legend_loc='right margin')
    
    # Customize
    ax.set_title(f'Leiden_new Cluster {cluster}', fontsize=14, fontweight='bold')
    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    
    # Save individual plot
    individual_path = os.path.join(output_dir, f'leiden_new_cluster_{cluster}_individual.png')
    plt.savefig(individual_path, dpi=300, bbox_inches='tight')
    print(f"Saved individual plot for cluster {cluster}: {individual_path}")
    plt.close()

# %%
# Create a summary plot showing all clusters with different colors
fig, ax = plt.subplots(1, 1, figsize=(10, 8))

sc.pl.umap(merged_adata, 
           color='leiden_new',
           ax=ax,
           show=False,
           frameon=True,
           title='All Leiden_new Clusters',
           legend_loc='right margin')

ax.set_title('All Leiden_new Clusters - Overview', fontsize=14, fontweight='bold')
ax.set_xlabel('UMAP 1', fontsize=12)
ax.set_ylabel('UMAP 2', fontsize=12)

# Save overview plot
overview_path = os.path.join(output_dir, 'all_leiden_new_clusters_overview.png')
plt.savefig(overview_path, dpi=300, bbox_inches='tight')
print(f"Saved overview plot: {overview_path}")
plt.show()

# %%
# Print summary
print("\n=== SUMMARY ===")
print(f"Total clusters analyzed: {len(unique_clusters)}")
print(f"Clusters: {unique_clusters}")
print(f"Output directory: {output_dir}")
print(f"Files created:")
print(f"  - Grid plot: individual_leiden_new_clusters_grid.png")
print(f"  - Overview plot: all_leiden_new_clusters_overview.png")
for cluster in unique_clusters:
    print(f"  - Individual plot: leiden_new_cluster_{cluster}_individual.png")

print("\nScript completed successfully!")

# %%