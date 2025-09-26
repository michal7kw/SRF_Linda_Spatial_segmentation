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
# for ROIs defined in DATA/ROIs
data_map = {
    'p0-p7': {
        'R1': {'time_point': 'p7', 'region': 'R1', 'type': 'mut'}, #ok
        'R2': {'time_point': 'p7', 'region': 'R2', 'type': 'ctrl'}, #ok
        # 'R3': {'time_point': 'p0', 'region': 'R3a', 'type': 'mut'},
        'R4': {'time_point': 'p30', 'region': 'R4p', 'type': 'ctrl'}, #ok
        'R5': {'time_point': 'p0', 'region': 'R5p', 'type': 'mut'}, #ok
        'R6': {'time_point': 'p0', 'region': 'R6', 'type': 'ctrl'} #ok
    },
    'p30-E165': {
        'R1': {'time_point': 'e16', 'region': 'R1', 'type': 'mut'}, #ok
        'R2': {'time_point': 'e16', 'region': 'R2', 'type': 'ctrl'}, #todo
        'R3': {'time_point': 'p30', 'region': 'R3', 'type': 'mut'} #ok
        # 'R4': {'time_point': 'p30', 'region': 'R4a', 'type': 'ctrl'}
    }
}

adatas = []
base_data_dir = 'DATA'


# %%
# Load each dataset and add metadata
print("Starting data loading and metadata assignment...")
for time_group, regions in data_map.items():
    print(f"Processing time group: {time_group}")
    for region, metadata in regions.items():
        print(f"  Processing region: {region}")
        h5ad_path = os.path.join(base_data_dir, time_group, region, 'adata_by_region_new.h5ad')
        if os.path.exists(h5ad_path):
            print(f"    Loading data from: {h5ad_path}")
            adata = sc.read_h5ad(h5ad_path).copy() # Explicitly copy to avoid ImplicitModificationWarning
            
            print(f"    Original data shape: {adata.shape}")
            
            # Check roi_assignment distribution before filtering
            if 'roi_assignment' in adata.obs.columns:
                roi_counts_before = adata.obs['roi_assignment'].value_counts()
                print(f"    ROI assignment distribution: {dict(roi_counts_before)}")
            
            print(f"    KEEPING ALL CELLS from ROI selection")
            # Don't filter - keep all cells from the ROI selection
            
            adata.obs['time_point'] = metadata['time_point']
            adata.obs['region'] = metadata['region']
            adata.obs['type'] = metadata['type']
            adata.obs['sample_id'] = f"{metadata['time_point']}_{metadata['region']}_{metadata['type']}"
            adatas.append(adata)
            sample_id = f"{metadata['time_point']}_{metadata['region']}_{metadata['type']}"
            print(f"    Final shape for {sample_id}: {adata.shape}")
            
        else:
            print(f"    File not found, skipping: {h5ad_path}")

if not adatas:
    print("No data loaded. Exiting.")
    sys.exit()
print("Data loading and metadata assignment completed.")

# %%
# Print summary before merging
print("\n" + "="*60)
print("SUMMARY BEFORE MERGING:")
total_cells = sum(adata.n_obs for adata in adatas)
print(f"Total cells across all samples: {total_cells}")
for i, adata in enumerate(adatas):
    sample_id = adata.obs['sample_id'].iloc[0]
    print(f"  {sample_id}: {adata.n_obs} cells")
print("="*60)

# %%
# Merge all AnnData objects
print("Merging AnnData objects...")
merged_adata = anndata.concat(adatas, join='outer', index_unique='-')
print(f"Merged data shape: {merged_adata.shape}")
print("AnnData objects merged.")

# %%
merged_adata

# %%
# calculate qc metrics
print("Calculating QC metrics...")
sc.pp.calculate_qc_metrics(merged_adata, percent_top=None, log1p=False, inplace=True)
print("QC metrics calculated.")

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
output_dir = 'MERGE_SPATIAL_DATA/plots_hippo_extended'
os.makedirs(output_dir, exist_ok=True)
print(f"Saving QC plots to {output_dir}")

# Set scanpy's figure directory
sc.settings.figdir = output_dir

# %%
print("Generating violin plots...")
# Violin plots for n_genes, n_counts, and total_counts across samples
sc.pl.violin(merged_adata, keys=['n_genes'], groupby='sample_id', rotation=90, save='_qc_violin_n_genes.png', show=False)
sc.pl.violin(merged_adata, keys=['n_counts'], groupby='sample_id', rotation=90, save='_qc_violin_n_counts.png', show=False)
sc.pl.violin(merged_adata, keys=['total_counts'], groupby='sample_id', rotation=90, save='_qc_violin_total_counts.png', show=False)
print("Violin plots generated.")

# %%
# Pre-processing and dimensionality reduction on the merged data
print("Performing PCA...")
sc.tl.pca(merged_adata, svd_solver='arpack')
print("PCA completed.")

# %%
print("Computing neighbors...")
sc.pp.neighbors(merged_adata, n_neighbors=15, n_pcs=40)
print("Neighbors computed.")

# %%
print("Running Leiden clustering...")
sc.tl.leiden(merged_adata, key_added='leiden_new', resolution=0.5)
print("Leiden clustering completed.")

# %%
print("Running UMAP...")
# Key UMAP parameters for better separation
sc.tl.umap(merged_adata,
           min_dist=0.1,      # Increase for more separation (0.1-1.0)
           spread=1.0,        # Increase for more spread out embedding
           n_components=2,
           gamma=2.0,         # Repulsion strength - increase for separation
           random_state=42)
print("UMAP completed.")

# %%
# Plotting
print(f"Saving UMAP plots to {output_dir}")

print("Generating UMAP plots...")
# Plot UMAP colored by clusters (leiden)
sc.pl.umap(merged_adata, color='leiden', save='_leiden_clusters.png', show=False)

# Plot UMAP colored by clusters (leiden_new)
sc.pl.umap(merged_adata, color='leiden_new', save='_leiden_new_clusters.png', show=False)

# Plot UMAP colored by time points
sc.pl.umap(merged_adata, color='time_point', save='_time_points.png', show=False)

# Plot UMAP colored by regions
sc.pl.umap(merged_adata, color='region', save='_regions.png', show=False)
print("UMAP plots generated.")

# %%
# Save the merged data
output_h5ad = f'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_all_cells.h5ad'
print(f"Saving merged data to {output_h5ad}")
merged_adata.write(output_h5ad)
print("Merged data saved.")

# %%
# Final summary
print("\n" + "="*60)
print("FINAL SUMMARY:")
print(f"Final merged data shape: {merged_adata.shape}")
print(f"Total cells in final dataset: {merged_adata.n_obs}")

if 'roi_assignment' in merged_adata.obs.columns:
    final_roi_counts = merged_adata.obs['roi_assignment'].value_counts()
    print(f"Final ROI assignment distribution: {dict(final_roi_counts)}")

print(f"Data saved to: {output_h5ad}")
print("="*60)

print("Script finished successfully.")