# %%
import sys
import os
import scanpy as sc
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely import wkt
from matplotlib.patches import Patch

project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
os.chdir(project_root)

# --- 1. Load Data ---
sample_id = 'p0_R6_ctrl'
time_group = 'p0-p7'
region = 'R6'
base_data_dir = 'DATA'

h5ad_path = os.path.join(base_data_dir, time_group, region, 'data.h5ad')
adata = sc.read_h5ad(h5ad_path)
print(f"Shape of adata after loading: {adata.shape}")

roi_path = os.path.join(base_data_dir, "ROIs", f"{sample_id}_geometry.csv")
roi_polygon = pd.read_csv(roi_path)
roi_polygon['geometry'] = roi_polygon['geometry'].apply(wkt.loads)
roi_gdf = gpd.GeoDataFrame(roi_polygon, geometry='geometry')
print(f"Shape of roi_gdf after loading: {roi_gdf.shape}")

cell_boundaries_path = os.path.join(base_data_dir, time_group, region, 'cell_boundaries.parquet')
cell_boundaries_gdf = gpd.read_parquet(cell_boundaries_path)
print(f"Shape of cell_boundaries_gdf after loading: {cell_boundaries_gdf.shape}")
cell_boundaries_gdf['EntityID'] = cell_boundaries_gdf['EntityID'].astype(str).str.replace('-4$', '', regex=True)

# --- 2. Geometric Filtering ---
cells_in_rois_gdf = gpd.sjoin(cell_boundaries_gdf, roi_gdf, how="inner", predicate="intersects")
print(f"Shape of cells_in_rois_gdf after sjoin: {cells_in_rois_gdf.shape}")
cells_in_rois_gdf['EntityID_left'] = cells_in_rois_gdf['EntityID_left'].astype(str).str.replace('-4$', '', regex=True)
cells_to_keep = cells_in_rois_gdf['EntityID_left']
print(f"Number of cells to keep after ROI filtering: {len(cells_to_keep)}")
adata.obs.index = adata.obs.index.astype(str).str.replace('-4$', '', regex=True)
adata_filtered = adata[adata.obs.index.isin(cells_to_keep)].copy()
print(f"Shape of adata_filtered after ROI filtering: {adata_filtered.shape}")

# --- 3. Clustering ---
n_pcs = min(adata_filtered.n_obs - 1, adata_filtered.n_vars - 1, 50)
if n_pcs < 1:
    n_pcs = 1
sc.pp.pca(adata_filtered, n_comps=n_pcs)
sc.pp.neighbors(adata_filtered, n_neighbors=15, n_pcs=n_pcs)
sc.tl.leiden(adata_filtered, key_added='leiden', resolution=0.5)

# --- 4. Plotting ---
output_dir = "ROIs_SELECTION/plots_missing_cell_exploration"
os.makedirs(output_dir, exist_ok=True)

fig, ax = plt.subplots(1, 1, figsize=(12, 12), dpi=150)

# 1. Plot cells OUTSIDE the ROI as background
cells_outside_roi_gdf = cell_boundaries_gdf[~cell_boundaries_gdf['EntityID'].isin(cells_in_rois_gdf['EntityID_left'])]
cells_outside_roi_sample = cells_outside_roi_gdf.sample(frac=0.1) # Sample to avoid plotting too many cells
cells_outside_roi_sample.plot(ax=ax, color='lightgray', edgecolor='silver', alpha=0.3, label='Outside ROI')

# 2. Plot cells INSIDE the ROI that are NOT clustered (i.e., not in adata_filtered)
# These are cells in cells_in_rois_gdf but not in adata_filtered.obs.index
cells_in_roi_not_clustered_ids = cells_in_rois_gdf[~cells_in_rois_gdf['EntityID_left'].isin(adata_filtered.obs.index)]['EntityID_left']
cells_in_roi_not_clustered_gdf = cells_in_rois_gdf[cells_in_rois_gdf['EntityID_left'].isin(cells_in_roi_not_clustered_ids)]
cells_in_roi_not_clustered_gdf.plot(ax=ax, color='lightcoral', edgecolor='red', alpha=0.5, label='Unclustered in ROI')

# 3. Plot clustered cells
leiden_categories = adata_filtered.obs['leiden'].cat.categories
colors = [plt.cm.tab20(i/20) for i in range(len(leiden_categories))]
for i, cluster in enumerate(leiden_categories):
    mask = adata_filtered.obs['leiden'] == cluster
    print(f"Cluster {cluster}: Number of cells in cluster mask: {mask.sum()}")
    if mask.sum() > 0:
        cells_in_cluster = adata_filtered[mask].obs.index
        cells_to_plot = cells_in_rois_gdf[cells_in_rois_gdf['EntityID_left'].isin(cells_in_cluster)]
        print(f"Cluster {cluster}: Shape of cells_to_plot: {cells_to_plot.shape}")
        cells_to_plot.plot(ax=ax, color=colors[i], alpha=0.7, label=f'Cluster {cluster}')

# Plot ROI boundaries
roi_gdf.plot(ax=ax, facecolor='none', edgecolor='blue', linewidth=2)

ax.set_xlim(10500, 14500)
ax.set_ylim(0, 4000)
ax.set_aspect('equal', adjustable='box')

# Create custom legend handles
legend_handles = [Patch(color='lightcoral', label='Unclustered in ROI')]
for i, cluster in enumerate(leiden_categories):
    legend_handles.append(Patch(color=colors[i], label=f'Cluster {cluster}'))

ax.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2, fontsize=8)
plot_path = os.path.join(output_dir, f"ROI_selection_and_clustering_{sample_id}.png")
fig.savefig(plot_path, bbox_inches='tight')
plt.close(fig)