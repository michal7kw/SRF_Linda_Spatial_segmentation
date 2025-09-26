# %% [markdown]
# # MERSCOPE DEGs Analysis

# %%
# Import necessary libraries
import sys
import os
import scanpy as sc
import anndata as ad
import pandas as pd
import geopandas as gpd # For .parquet file with geometries
import matplotlib.pyplot as plt
import matplotlib.image as mpimg # For displaying images
import seaborn as sns
import squidpy
import warnings
from shapely import wkt
import numpy as np

project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
# %%
os.chdir(project_root)
plot_enable = True

# %%
# Define the data structure and metadata
# # for ROIs defined in /DATA/ROIs/DemoSelections
# data_map = {
#     'p0-p7': {
#         'R1': {'time_point': 'p7', 'region': 'R1', 'type': 'mut'},
#         'R2': {'time_point': 'p7', 'region': 'R2', 'type': 'ctrl'},
#         'R3': {'time_point': 'p0', 'region': 'R3a', 'type': 'mut'},
#         'R4': {'time_point': 'p30', 'region': 'R4p', 'type': 'ctrl'},
#         'R5': {'time_point': 'p0', 'region': 'R5p', 'type': 'mut'},
#         'R6': {'time_point': 'p0', 'region': 'R6', 'type': 'ctrl'}
#     },
#     'p30-E165': {
#         'R1': {'time_point': 'p16', 'region': 'R1', 'type': 'mut'},
#         'R2': {'time_point': 'p16', 'region': 'R2', 'type': 'ctrl'},
#         'R3': {'time_point': 'p30', 'region': 'R3', 'type': 'mut'},
#         'R4': {'time_point': 'p30', 'region': 'R4a', 'type': 'ctrl'}
#     }
# }
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
        'R2': {'time_point': 'e16', 'region': 'R2', 'type': 'ctrl'}, #ok
        'R3': {'time_point': 'p30', 'region': 'R3', 'type': 'mut'} #ok
        # 'R4': {'time_point': 'p30', 'region': 'R4a', 'type': 'ctrl'}
    }
}

adatas_abd_geos = {}
base_data_dir = 'DATA'
sample_id_to_path_info = {}

# %%
# Load each dataset and add metadata
for time_group, regions in data_map.items():
    print(f"time_group: {time_group}")
    for region, metadata in regions.items():
        print(f"region: {region}")
        
        sample_id = f"{metadata['time_point']}_{metadata['region']}_{metadata['type']}"
        sample_id_to_path_info[sample_id] = {'time_group': time_group, 'region': region}
        
        roi_path = os.path.join(base_data_dir, "ROIs", f"{sample_id}_geometry.csv")
        if os.path.exists(roi_path):
            print(f"Loading ROI data from: {roi_path}")
            roi_polygon = pd.read_csv(roi_path)
        else:
            print(f"File not found: {roi_path}")
            continue

        h5ad_path = os.path.join(base_data_dir, time_group, region, 'data.h5ad')
        if os.path.exists(h5ad_path):
            print(f"Loading Transcriptomics data from: {h5ad_path}")
            adata = sc.read_h5ad(h5ad_path)
            adata.obs['time_point'] = metadata['time_point']
            adata.obs['region'] = metadata['region']
            adata.obs['type'] = metadata['type']
            adata.obs['sample_id'] = sample_id
        else:
            print(f"File not found: {h5ad_path}")
            continue

        cell_boundaries_path = os.path.join(base_data_dir, time_group, region, 'cell_boundaries.parquet')
        if os.path.exists(cell_boundaries_path):
            print(f"Loading Cell Boundaries data from: {cell_boundaries_path}")
            cell_boundaries_gdf = gpd.read_parquet(cell_boundaries_path)
            cell_boundaries_gdf = cell_boundaries_gdf.set_index('EntityID', drop=False)
        else:
            print(f"File not found: {cell_boundaries_path}")
            continue

        adatas_abd_geos[sample_id] = (adata, roi_polygon, cell_boundaries_gdf)

# %%
print(adatas_abd_geos.keys())

# %% [markdown]
# # Process ROI geometries and perform spatial selection

processed_adatas_abd_geos = {}
for sample_id, (adata, roi, cell_boundaries_gdf) in adatas_abd_geos.items():
    # Convert WKT string geometries to Shapely geometry objects
    roi['geometry'] = roi['geometry'].apply(wkt.loads)

    current_crs = None
    current_crs = cell_boundaries_gdf.crs

    roi_gdf = gpd.GeoDataFrame(roi, geometry='geometry', crs=current_crs)
    print(f"Successfully converted ROI geometries for {sample_id} to GeoDataFrame.")

    cell_boundaries_gdf_sjoin_ready = cell_boundaries_gdf.copy()
    cell_boundaries_gdf_sjoin_ready.index.name = 'original_cell_EntityID_idx' # Rename the index

    cells_in_rois_gdf = gpd.sjoin(
        cell_boundaries_gdf_sjoin_ready,
        roi_gdf,
        how="inner",
        predicate="intersects", # Changed from "within" to "intersects" to include cells that partially overlap with the ROI
        lsuffix='_cell',
        rsuffix='_roi'
    )
    
    processed_adatas_abd_geos[sample_id] = (adata, roi, cell_boundaries_gdf, cells_in_rois_gdf, roi_gdf)

    print(f"\nFound {cells_in_rois_gdf.shape[0]} cells within the defined ROIs for {sample_id}.")

adatas_abd_geos = processed_adatas_abd_geos

# %%
print(adatas_abd_geos.keys())

# %%
# Visualize selected cells with ROIs
if plot_enable:
    for sample_id, (adata, roi, cell_boundaries_gdf, cells_in_rois_gdf, roi_gdf) in adatas_abd_geos.items():
        print("\nCell counts per ROI group (from 'group' column of ROI CSV):")
        print(cells_in_rois_gdf.groupby('group').size())

        print("\nPlotting selected cells and ROIs...")
        fig, ax = plt.subplots(1, 1, figsize=(12, 12), dpi=150) # Use dpi for rasterization quality

        # Plot a random sample of original cell boundaries lightly for efficiency
        background_cells_sample = cell_boundaries_gdf.sample(frac=0.1)
        background_cells_sample.plot(ax=ax, color='lightgray', edgecolor='silver', alpha=0.3, label='All Cells (Original Sample)', rasterized=True)

        cells_in_rois_gdf.plot(ax=ax, column='group', legend=True, alpha=0.7, categorical=True,
                                            legend_kwds={'title': "ROI Group", 'loc': 'upper right', 'bbox_to_anchor': (1.25, 1)},
                                            rasterized=True)

        # Plot ROI boundaries
        roi_gdf.plot(ax=ax, facecolor='none', edgecolor='blue', linewidth=2, label='ROI Polygons')

        ax.set_title(f"Cells within Defined ROIs for {sample_id}")
        ax.set_xlabel("X-coordinate")
        ax.set_ylabel("Y-coordinate")

        # Adjust legend display
        handles, labels = ax.get_legend_handles_labels()

        plt.tight_layout(rect=[0, 0, 0.85, 1])

        # Save the plot
        # plot_dir = "ROIs_SELECTION/plots_old_rois"
        plot_dir = "ROIs_SELECTION/plots_new_rois"
        os.makedirs(plot_dir, exist_ok=True)
        plot_path = os.path.join(plot_dir, f"ROI_selection_{sample_id}.png")
        fig.savefig(plot_path, bbox_inches='tight')
        print(f"Plot saved: {plot_path}")
        
        plt.show()
        plt.close(fig) # Close the figure to free memory

# %%
for sample_id, (adata, roi, cell_boundaries_gdf, cells_in_rois_gdf, roi_gdf) in adatas_abd_geos.items():
    # Create a mapping of cell EntityID to its ROI group
    cell_to_roi_map = cells_in_rois_gdf["group"].to_dict()
    for key, value in list(cell_to_roi_map.items())[:5]:
        print(f"{key}: {value}")

    adata.obs['roi_assignment'] = 'unassigned'

    # Keys for cell_to_roi_map
    cell_to_roi_map = {str(key): value for key, value in cell_to_roi_map.items()}
    map_keys = list(cell_to_roi_map.keys())

    # Perform intersection ensuring
    adata_index_as_str = adata.obs.index
    map_keys_as_str_for_intersection = pd.Index(map_keys)

    valid_intersected_ids_str = adata_index_as_str.intersection(map_keys_as_str_for_intersection)

    # Assign ROI based on the intersection.
    for cell_id_str_from_intersection in valid_intersected_ids_str:
        if cell_id_str_from_intersection in adata.obs.index:
                adata.obs.loc[cell_id_str_from_intersection, 'roi_assignment'] = cell_to_roi_map[cell_id_str_from_intersection]

    # Check how many cells were assigned
    print(f"Cells assigned to ROIs in adata.obs: \n{adata.obs['roi_assignment'].value_counts()}")



# %%
# Save all modified adata objects
print("\n--- Saving all modified adata objects ---")
for sample_id, (adata, _, _, _, _) in adatas_abd_geos.items():
    path_info = sample_id_to_path_info.get(sample_id)
    if not path_info:
        print(f"Warning: Could not find path info for {sample_id}. Skipping save.")
        continue

    save_path_dir = os.path.join(base_data_dir, path_info['time_group'], path_info['region'])
    # save_path = os.path.join(save_path_dir, 'adata_by_region.h5ad')
    save_path = os.path.join(save_path_dir, 'adata_by_region_new.h5ad')
    
    try:
        os.makedirs(save_path_dir, exist_ok=True)
        print(f"Saving adata for {sample_id} to {save_path}")
        adata.write_h5ad(save_path)
    except Exception as e:
        print(f"Error saving file for {sample_id}: {e}")

print("\n--- All adata objects saved ---")

# %%
