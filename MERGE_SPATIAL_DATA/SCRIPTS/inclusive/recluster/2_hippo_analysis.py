# %%
import scanpy as sc
import squidpy as sq
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Configure scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Define the path to the pre-computed reclustered data file
reclustered_h5ad_path = 'MERGE_SPATIAL_DATA/correlation_analysis_cluster0_recluster/merged_data_with_combined_cluster0_reclustering.h5ad'

# %%
# Load the pre-computed reclustered AnnData object
print(f"Loading pre-computed reclustered data from: {reclustered_h5ad_path}")
merged_adata = sc.read_h5ad(reclustered_h5ad_path)
print("Pre-computed reclustered data loaded successfully.")
print(f"Data shape: {merged_adata.shape}")
print(f"Available cluster columns: {[col for col in merged_adata.obs.columns if 'leiden_new' in col or 'cluster' in col]}")

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
output_dir = 'MERGE_SPATIAL_DATA/plots_hippo_recluster_highlight_inclusive'
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

print(f"Results will be saved to: {output_dir}")

# %%
# Using pre-computed reclustered data (assumed to always exist)
print("Using pre-computed reclustered data from correlation analysis.")
print(f"Available reclustered clusters: {sorted(merged_adata.obs['leiden_combined_cluster0_reclustered'].unique())}")

# %%
# Get unique samples
samples = merged_adata.obs['sample_id'].unique()
print(f"Found {len(samples)} samples: {samples}")

# %%
# Process samples using pre-computed reclustered data
print("\n" + "="*60)
print("USING PRE-COMPUTED RECLUSTERED DATA")
print("="*60)

recluster_summary = {}
all_reclustered_data = []

# Use pre-computed reclustered data
main_cluster_column = 'leiden_combined_cluster0_reclustered'

for sample_id in samples:
    print(f"\nProcessing sample: {sample_id}")
    print("-" * 40)
    
    # Filter data for current sample
    sample_adata = merged_adata[merged_adata.obs['sample_id'] == sample_id].copy()
    print(f"Sample {sample_id} shape: {sample_adata.shape}")
    
    # Remove unused categories from cluster columns
    if 'leiden_new' in sample_adata.obs.columns:
        sample_adata.obs['leiden_new'] = sample_adata.obs['leiden_new'].cat.remove_unused_categories()
    if main_cluster_column in sample_adata.obs.columns:
        sample_adata.obs[main_cluster_column] = sample_adata.obs[main_cluster_column].cat.remove_unused_categories()
    
    # Check cluster distribution
    if main_cluster_column in sample_adata.obs.columns:
        cluster_counts = sample_adata.obs[main_cluster_column].value_counts()
        print(f"Sample {sample_id} cluster counts ({main_cluster_column}):")
        for cluster, count in cluster_counts.sort_index().items():
            print(f"  Cluster {cluster}: {count} cells")
        
        # Analyze cluster types
        reclustered_clusters = [c for c in cluster_counts.index if str(c).startswith('0.')]
        original_clusters = [c for c in cluster_counts.index if not str(c).startswith('0.')]
        
        recluster_summary[sample_id] = {
            'had_cluster0': len(reclustered_clusters) > 0 or '0' in [str(c) for c in cluster_counts.index],
            'reclustering_successful': len(reclustered_clusters) > 0,
            'original_cluster0_cells': cluster_counts.get('0', 0) if '0' in [str(c) for c in cluster_counts.index] else sum(cluster_counts[c] for c in reclustered_clusters),
            'reclustered_subclusters': reclustered_clusters,
            'n_subclusters': len(reclustered_clusters),
            'cluster_column': main_cluster_column
        }
        
        if len(reclustered_clusters) > 0:
            print(f"  Found {len(reclustered_clusters)} reclustered subclusters: {reclustered_clusters}")
    else:
        print(f"Warning: Cluster column {main_cluster_column} not found in sample {sample_id}")
        recluster_summary[sample_id] = {
            'had_cluster0': False,
            'reclustering_successful': False,
            'original_cluster0_cells': 0,
            'reclustered_subclusters': [],
            'n_subclusters': 0,
            'cluster_column': 'leiden_new'
        }
    
    # Add metadata for tracking
    sample_adata.obs['recluster_info'] = f"{sample_id}_{main_cluster_column}"
    all_reclustered_data.append(sample_adata)

# %%
# Combine all reclustered sample data
print("\nCombining reclustered data from all samples...")
merged_reclustered_adata = all_reclustered_data[0].concatenate(
    all_reclustered_data[1:] if len(all_reclustered_data) > 1 else [],
    batch_key='batch_sample',
    batch_categories=[f"sample_{i}" for i in range(len(all_reclustered_data))]
)

# Create unified cluster column using pre-computed reclustered data
merged_reclustered_adata.obs['unified_clusters'] = merged_reclustered_adata.obs[main_cluster_column].astype(str).copy()
print(f"Using pre-computed cluster column: {main_cluster_column}")

merged_reclustered_adata.obs['unified_clusters'] = merged_reclustered_adata.obs['unified_clusters'].astype('category')

print(f"Combined data shape: {merged_reclustered_adata.shape}")
print(f"Unified clusters found: {sorted(merged_reclustered_adata.obs['unified_clusters'].unique())}")

# %%
# Analyze cluster sizes per sample for reporting (including reclustered subclusters)
print("\nAnalyzing reclustered cluster sizes per sample...")
cluster_summary = {}

for sample_id in samples:
    sample_data = merged_reclustered_adata[merged_reclustered_adata.obs['sample_id'] == sample_id]
    cluster_counts = sample_data.obs['unified_clusters'].value_counts().sort_index()
    
    # Separate original clusters from reclustered subclusters
    original_clusters = [c for c in cluster_counts.index if not c.startswith('0.')]
    reclustered_subclusters = [c for c in cluster_counts.index if c.startswith('0.')]
    
    small_clusters = cluster_counts[cluster_counts < 50].index.tolist()
    large_clusters = cluster_counts[cluster_counts >= 50].index.tolist()
    
    cluster_summary[sample_id] = {
        'original_clusters': original_clusters,
        'reclustered_subclusters': reclustered_subclusters,
        'small_clusters': small_clusters,
        'large_clusters': large_clusters,
        'total_clusters': len(cluster_counts)
    }
    
    print(f"\nSample {sample_id}:")
    print(f"  Total clusters: {len(cluster_counts)}")
    print(f"  Original clusters: {len(original_clusters)} - {original_clusters}")
    print(f"  Reclustered subclusters: {len(reclustered_subclusters)} - {reclustered_subclusters}")
    print(f"  Large clusters (≥50 cells): {len(large_clusters)} - {large_clusters}")
    print(f"  Small clusters (<50 cells): {len(small_clusters)} - {small_clusters}")

# %%
# Generate and save plots for each sample and cluster (including reclustered subclusters)
print("\n" + "="*60)
print("GENERATING INDIVIDUAL SUBCLUSTER PLOTS (RECLUSTERED ONLY)")
print("="*60)

for sample in samples:
    sample_adata = merged_reclustered_adata[merged_reclustered_adata.obs['sample_id'] == sample]
    
    # Get only reclustered subclusters for this sample (clusters starting with '0.')
    all_clusters = sample_adata.obs['unified_clusters'].unique()
    sample_clusters = sorted([c for c in all_clusters if c.startswith('0.')], key=lambda x: (len(x.split('.')), x))
    
    if len(sample_clusters) == 0:
        print(f"\nSample {sample}: No reclustered subclusters found, skipping plotting")
        continue
    
    print(f"\nProcessing sample {sample} with {len(sample_clusters)} reclustered subclusters: {sample_clusters}")
    
    for cluster in sample_clusters:
        cluster_mask = sample_adata.obs['unified_clusters'] == cluster
        cluster_cell_count = cluster_mask.sum()
        
        print(f"  Plotting sample {sample}, cluster {cluster} ({cluster_cell_count} cells)")

        # Create a figure with two subplots
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))
        
        # Get the color for this cluster
        if 'unified_clusters_colors' in merged_reclustered_adata.uns:
            cluster_colors = dict(zip(merged_reclustered_adata.obs['unified_clusters'].cat.categories, 
                                    merged_reclustered_adata.uns['unified_clusters_colors']))
            cluster_color = cluster_colors.get(cluster, 'blue')
        elif 'leiden_new_colors' in merged_reclustered_adata.uns and not cluster.startswith('0.'):
            # Use original leiden_new colors for non-reclustered clusters
            cluster_colors = dict(zip(merged_reclustered_adata.obs['leiden_new'].cat.categories, 
                                    merged_reclustered_adata.uns['leiden_new_colors']))
            cluster_color = cluster_colors.get(cluster, 'blue')
        else:
            # Generate a color for reclustered subclusters
            if cluster.startswith('0.'):
                # Use different shades of red for cluster 0 subclusters
                subcluster_idx = int(cluster.split('.')[1]) if '.' in cluster else 0
                colors = ['#FF0000', '#FF4444', '#FF8888', '#FFAAAA', '#FFCCCC', '#FF6666', '#FF2222']
                cluster_color = colors[subcluster_idx % len(colors)]
            else:
                cluster_color = 'blue'
        
        # UMAP plot
        umap_coords = merged_reclustered_adata.obsm['X_umap']
        sample_mask = merged_reclustered_adata.obs['sample_id'] == sample
        full_cluster_mask = (merged_reclustered_adata.obs['unified_clusters'] == cluster) & sample_mask
        
        # Plot all cells in gray on the first subplot
        axes[0].scatter(umap_coords[:, 0], umap_coords[:, 1], c='lightgray', s=1, alpha=0.3)
        
        # Plot highlighted sample-cluster combination in color on the first subplot
        if full_cluster_mask.sum() > 0:
            axes[0].scatter(umap_coords[full_cluster_mask, 0], umap_coords[full_cluster_mask, 1],
                          c=cluster_color, s=5, alpha=0.8)
        
        axes[0].set_xlabel('UMAP1')
        axes[0].set_ylabel('UMAP2')
        
        # Title formatting for reclustered subclusters only
        axes[0].set_title(f"UMAP: {sample} - Subcluster {cluster}\n(Reclustered from Cluster 0)", fontsize=10)
        
        # Spatial plot using matplotlib directly on the second subplot
        spatial_coords = sample_adata.obsm['spatial']
        
        # Plot all cells in gray on the second subplot
        axes[1].scatter(spatial_coords[:, 0], spatial_coords[:, 1], c='lightgray', s=1, alpha=0.5)
        
        # Plot highlighted cluster in the same color as UMAP on the second subplot
        axes[1].scatter(spatial_coords[cluster_mask, 0], spatial_coords[cluster_mask, 1],
                      c=cluster_color, s=3, alpha=1.0)
        
        axes[1].set_xlabel('X coordinate')
        axes[1].set_ylabel('Y coordinate')
        
        # Title formatting for reclustered subclusters only
        axes[1].set_title(f"Spatial: {sample} - Subcluster {cluster}\n(Reclustered from Cluster 0, {cluster_cell_count} cells)", fontsize=10)
        
        axes[1].set_aspect('equal')
        
        # Save the figure with subcluster naming
        plot_filename = f"sample_{sample}_subcluster_{cluster.replace('.', '_')}.png"
        
        fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight')
        plt.close(fig)

# %%
# Create a single UMAP plot highlighting only original cluster 0 cells
print("\n" + "="*60)
print("GENERATING SINGLE UMAP PLOT WITH CLUSTER 0 HIGHLIGHTED")
print("="*60)

# Create figure for the overview UMAP plot
fig, ax = plt.subplots(1, 1, figsize=(10, 8))

# Get UMAP coordinates for all cells
umap_coords = merged_reclustered_adata.obsm['X_umap']

# Plot all cells in light gray
ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c='lightgray', s=1, alpha=0.3, label='Other clusters')

# Highlight original cluster 0 cells (before reclustering) in red
original_cluster0_mask = merged_reclustered_adata.obs['leiden_new'] == '0'
if original_cluster0_mask.sum() > 0:
    ax.scatter(umap_coords[original_cluster0_mask, 0], umap_coords[original_cluster0_mask, 1], 
              c='red', s=3, alpha=0.8, label='Original Cluster 0')
    print(f"Highlighted {original_cluster0_mask.sum()} cells from original cluster 0")
else:
    print("No original cluster 0 cells found to highlight")

ax.set_xlabel('UMAP1', fontsize=12)
ax.set_ylabel('UMAP2', fontsize=12)
ax.set_title('UMAP Overview: Original Cluster 0 Highlighted\n(Before Reclustering)', fontsize=14)
ax.legend()

# Save the overview plot
overview_filename = 'umap_overview_cluster0_highlighted.png'
fig.savefig(os.path.join(output_dir, overview_filename), bbox_inches='tight', dpi=300)
plt.close(fig)

print(f"Overview UMAP plot saved as: {overview_filename}")

# %%
# Create a reclustered analysis summary report
print("\n" + "="*60)
print("PRE-COMPUTED RECLUSTERED ANALYSIS SUMMARY")
print("="*60)

total_samples_with_cluster0 = 0
total_samples_reclustered = 0
total_subclusters_created = 0
total_clusters_analyzed = 0

for sample_id, info in cluster_summary.items():
    recluster_info = recluster_summary[sample_id]
    
    print(f"\nSample: {sample_id}")
    print(f"  Total clusters found: {info['total_clusters']}")
    print(f"  Original clusters: {len(info['original_clusters'])} - {info['original_clusters']}")
    
    if recluster_info['had_cluster0']:
        total_samples_with_cluster0 += 1
        print(f"  Had cluster 0: Yes ({recluster_info['original_cluster0_cells']} cells)")
        
        if recluster_info['reclustering_successful']:
            total_samples_reclustered += 1
            total_subclusters_created += recluster_info['n_subclusters']
            print(f"  Reclustering successful: Yes")
            print(f"  Subclusters created: {recluster_info['n_subclusters']} - {recluster_info['reclustered_subclusters']}")
        else:
            print(f"  Reclustering successful: No")
    else:
        print(f"  Had cluster 0: No")
    
    print(f"  Large clusters (≥50 cells): {len(info['large_clusters'])} - {info['large_clusters']}")
    print(f"  Small clusters (<50 cells): {len(info['small_clusters'])} - {info['small_clusters']}")
    
    total_clusters_analyzed += info['total_clusters']

print(f"\nOVERALL SUMMARY:")
print(f"  Total samples processed: {len(samples)}")
print(f"  Samples with cluster 0: {total_samples_with_cluster0}")
print(f"  Samples successfully reclustered: {total_samples_reclustered}")
print(f"  Total subclusters created: {total_subclusters_created}")
print(f"  Total clusters analyzed across all samples: {total_clusters_analyzed}")

# Save reclustered analysis summary to file
summary_file = os.path.join(output_dir, 'reclustered_analysis_summary.txt')
with open(summary_file, 'w') as f:
    f.write("PRE-COMPUTED RECLUSTERED ANALYSIS SUMMARY\n")
    f.write("=" * 60 + "\n\n")
    f.write("This file shows results from analysis using pre-computed reclustered data.\n")
    f.write("The reclustered data was loaded from: merged_data_with_combined_cluster0_reclustering.h5ad\n")
    f.write("Cluster 0 cells were previously reclustered into subclusters (0.0, 0.1, 0.2, etc.)\n")
    f.write("and individual plots were generated ONLY for the new subclusters.\n")
    f.write("Original clusters (1, 2, 3, etc.) were not plotted.\n\n")
    
    for sample_id, info in cluster_summary.items():
        recluster_info = recluster_summary[sample_id]
        
        f.write(f"Sample: {sample_id}\n")
        f.write(f"  Total clusters found: {info['total_clusters']}\n")
        f.write(f"  Original clusters: {len(info['original_clusters'])} - {info['original_clusters']}\n")
        
        if recluster_info['had_cluster0']:
            f.write(f"  Had cluster 0: Yes ({recluster_info['original_cluster0_cells']} cells)\n")
            
            if recluster_info['reclustering_successful']:
                f.write(f"  Reclustering successful: Yes\n")
                f.write(f"  Subclusters created: {recluster_info['n_subclusters']} - {recluster_info['reclustered_subclusters']}\n")
            else:
                f.write(f"  Reclustering successful: No\n")
        else:
            f.write(f"  Had cluster 0: No\n")
        
        f.write(f"  Large clusters analyzed (≥50 cells): {len(info['large_clusters'])} - {info['large_clusters']}\n")
        f.write(f"  Small clusters analyzed (<50 cells): {len(info['small_clusters'])} - {info['small_clusters']}\n\n")
    
    f.write(f"OVERALL SUMMARY:\n")
    f.write(f"  Total samples processed: {len(samples)}\n")
    f.write(f"  Samples with cluster 0: {total_samples_with_cluster0}\n")
    f.write(f"  Samples successfully reclustered: {total_samples_reclustered}\n")
    f.write(f"  Total subclusters created: {total_subclusters_created}\n")
    f.write(f"  Total clusters analyzed: {total_clusters_analyzed}\n\n")
    
    f.write("PRE-COMPUTED DATA SOURCE:\n")
    f.write(f"  - Loaded from: {reclustered_h5ad_path}\n")
    f.write(f"  - Cluster column used: {main_cluster_column}\n")
    f.write("  - Reclustering was performed in the correlation analysis script\n")

print(f"\nAnalysis summary saved to: {summary_file}")
print("Script finished successfully.")

# %%