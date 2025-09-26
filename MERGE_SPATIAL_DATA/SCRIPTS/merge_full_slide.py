# %%
import os
import sys
import scanpy as sc
import pandas as pd
import anndata

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# %%
# Define the data structure and metadata
data_map = {
    'p0-p7': {
        'R1': {'time_point': 'p7', 'region': 'R1', 'type': 'mut'},
        'R2': {'time_point': 'p7', 'region': 'R2', 'type': 'ctrl'},
        'R3': {'time_point': 'p0', 'region': 'R3a', 'type': 'mut'},
        'R4': {'time_point': 'p30', 'region': 'R4p', 'type': 'ctrl'},
        'R5': {'time_point': 'p0', 'region': 'R5p', 'type': 'mut'},
        'R6': {'time_point': 'p0', 'region': 'R6', 'type': 'ctrl'}
    },
    'p30-E165': {
        'R1': {'time_point': 'p16', 'region': 'R1', 'type': 'mut'},
        'R2': {'time_point': 'p16', 'region': 'R2', 'type': 'ctrl'},
        'R3': {'time_point': 'p30', 'region': 'R3', 'type': 'mut'},
        'R4': {'time_point': 'p30', 'region': 'R4a', 'type': 'ctrl'}
    }
}

adatas = []
base_data_dir = 'DATA'

# %%
# Load each dataset and add metadata
for time_group, regions in data_map.items():
    print(f"region: {time_group}")
    for region, metadata in regions.items():
        print(f"region: {region}")
        h5ad_path = os.path.join(base_data_dir, time_group, region, 'data.h5ad')
        if os.path.exists(h5ad_path):
            print(f"Loading data from: {h5ad_path}")
            adata = sc.read_h5ad(h5ad_path)
            adata.obs['time_point'] = metadata['time_point']
            adata.obs['region'] = metadata['region']
            adata.obs['type'] = metadata['type']
            adata.obs['sample_id'] = f"{metadata['time_point']}_{metadata['region']}_{metadata['type']}"
            adatas.append(adata)
        else:
            print(f"File not found, skipping: {h5ad_path}")

if not adatas:
    print("No data loaded. Exiting.")
    sys.exit()

# %%
# Merge all AnnData objects
print("Merging AnnData objects...")
merged_adata = anndata.concat(adatas, join='outer', index_unique='-')

# %%
merged_adata

# %%
# calculate qc metrics
sc.pp.calculate_qc_metrics(merged_adata, percent_top=None, log1p=False, inplace=True)


# %%
print(f"n_genes: {merged_adata.obs.n_genes[:10]}\n")
print(f"n_counts: {merged_adata.obs.n_counts[:10]}\n")
print(f"n_genes_by_counts: {merged_adata.obs.n_genes_by_counts[:10]}\n")
print(f"total_counts: {merged_adata.obs.total_counts[:10]}\n")

# %%
print(f"n_cells_by_counts: {merged_adata.var.n_cells_by_counts[:10]}\n")
print(f"mean_counts: {merged_adata.var.mean_counts[:10]}\n")
print(f"total_counts: {merged_adata.var.total_counts[:10]}\n")

# %%
# QC metric visualizations
output_dir = 'MERGE_SPATIAL_DATA/plots_full_slide'
os.makedirs(output_dir, exist_ok=True)
print(f"Saving QC plots to {output_dir}")

# %%
# Violin plots for n_genes, n_counts, and total_counts across samples
sc.pl.violin(merged_adata, keys=['n_genes'], groupby='sample_id', rotation=90, save='_qc_violin_n_genes.png', show=True)
sc.pl.violin(merged_adata, keys=['n_counts'], groupby='sample_id', rotation=90, save='_qc_violin_n_counts.png', show=True)
sc.pl.violin(merged_adata, keys=['total_counts'], groupby='sample_id', rotation=90, save='_qc_violin_total_counts.png', show=True)

# %%
# Pre-processing and dimensionality reduction on the merged data
print("Performing PCA...")
sc.tl.pca(merged_adata, svd_solver='arpack')

# %%
print("Computing neighbors...")
sc.pp.neighbors(merged_adata, n_neighbors=10, n_pcs=40)

# %%
print("Running UMAP...")
sc.tl.umap(merged_adata)

# %%
# Plotting
print(f"Saving UMAP plots to {output_dir}")

# %%
# Plot UMAP colored by clusters (leiden)
sc.pl.umap(merged_adata, color='leiden', save='_leiden_clusters.png', show=False)

# %%
# Plot UMAP colored by time points
sc.pl.umap(merged_adata, color='time_point', save='_time_points.png', show=False)

# %%
# Plot UMAP colored by regions
sc.pl.umap(merged_adata, color='region', save='_regions.png', show=False)

# %%
# Move plots to the correct directory
plot_files_to_move = [
    'umap_leiden_clusters.png', 'umap_time_points.png', 'umap_regions.png',
    'violin_qc_violin.png', 'scatter_qc_scatter_total_vs_ngenes.png',
    'scatter_qc_scatter_total_vs_mt.png', 'scatter_qc_scatter_ngenes_vs_mt.png',
    'violin_qc_violin_n_genes.png', 'violin_qc_violin_n_counts.png', 'violin_qc_violin_total_counts.png'
]
for plot_file in plot_files_to_move:
    source_path = os.path.join('figures', plot_file)
    dest_path = os.path.join(output_dir, plot_file)
    if os.path.exists(source_path):
        os.rename(source_path, dest_path)
        print(f"Moved {plot_file} to {output_dir}")
    else:
        print(f"Plot file not found, skipping: {source_path}")

# %%
# Save the merged data
output_h5ad = 'MERGE_SPATIAL_DATA/merged_spatial_data_full_slide.h5ad'
print(f"Saving merged data to {output_h5ad}")
merged_adata.write(output_h5ad)

print("Script finished successfully.")
# %%
