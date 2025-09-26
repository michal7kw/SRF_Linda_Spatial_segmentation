import sys
import os
import pandas as pd
import geopandas as gpd
from shapely import wkt
import scanpy as sc

project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
os.chdir(project_root)

# --- 1. Load Data ---
sample_id = 'p0_R6_ctrl'
time_group = 'p0-p7'
region = 'R6'
base_data_dir = 'DATA'

h5ad_path = os.path.join(base_data_dir, time_group, region, 'data.h5ad')
adata = sc.read_h5ad(h5ad_path)

roi_path = os.path.join(base_data_dir, "ROIs", f"{sample_id}_geometry.csv")
roi_polygon = pd.read_csv(roi_path)
roi_polygon['geometry'] = roi_polygon['geometry'].apply(wkt.loads)
roi_gdf = gpd.GeoDataFrame(roi_polygon, geometry='geometry')

cell_boundaries_path = os.path.join(base_data_dir, time_group, region, 'cell_boundaries.parquet')
cell_boundaries_gdf = gpd.read_parquet(cell_boundaries_path)

# --- 2. Geometric Filtering ---
print("ROI GDF head:")
print(roi_gdf.head())
print("\nCell Boundaries GDF head:")
print(cell_boundaries_gdf.head())

cells_in_rois_gdf = gpd.sjoin(cell_boundaries_gdf, roi_gdf, how="inner", predicate="intersects")

print(f"\nNumber of cells in ROI: {len(cells_in_rois_gdf)}")
print(cells_in_rois_gdf.columns)

adata.obs.index = adata.obs.index.astype(str).str.replace('-4$', '', regex=True)
cells_to_keep = cells_in_rois_gdf['EntityID_left'].astype(str)
adata_filtered = adata[adata.obs.index.isin(cells_to_keep)].copy()

print(f"\nNumber of cells in filtered adata: {len(adata_filtered)}")