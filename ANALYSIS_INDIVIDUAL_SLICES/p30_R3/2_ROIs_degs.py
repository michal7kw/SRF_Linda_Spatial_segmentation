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
sys.path.append(project_root)
from HELPERs.plotting import create_volcano_plot

# Suppress FutureWarning messages
warnings.filterwarnings('ignore', category=FutureWarning)

# %%
os.chdir(project_root)

base_path = 'DATA/p30-E165/R3'
h5ad_file = os.path.join(base_path, 'data.h5ad')

# %%
roi_csv_file_name = 'p30_R3_ROI_28-05-25_16-57_geometry.csv'
roi_geometry_file_path = os.path.join(base_path, roi_csv_file_name)
summary_image_file = os.path.join(base_path, 'summary.png')
plot_enable = True

# %% [markdown]
# ## 1. Data Loading

# %%
adata = sc.read_h5ad(h5ad_file)

# %%
keep_genes = [x for x in adata.var.index.tolist() if 'Blank' not in x]
print(len(keep_genes))
print(adata.shape[1])

# %%
min_expression = 25
ser_exp = adata.to_df().sum(axis=1)

keep_cells = ser_exp[ser_exp > min_expression].index.tolist()
print(len(keep_cells))
print(adata.shape[0])

# adata = adata[keep_cells]
# adata

# %%
# Cell boundaries
cell_boundaries_file = os.path.join(base_path, 'cell_boundaries.parquet')
cell_boundaries_gdf = None

# %%
cell_boundaries_gdf = gpd.read_parquet(cell_boundaries_file)
print(f"Loaded {cell_boundaries_file}. Shape: {cell_boundaries_gdf.shape}")

# %%
cell_boundaries_gdf.head()

# %%
cell_boundaries_gdf = cell_boundaries_gdf.set_index('EntityID', drop=False)
cell_boundaries_gdf.head()

# %% [markdown]
# # Subselect cells based on ROI polygons from CSV
# 
# This section loads ROI polygons from a CSV file, converts them to geometries,
# and then selects cells from `cell_boundaries_gdf` that fall within these ROIs.

# %%
# Load ROI geometry data from CSV
roi_polygons_df = pd.read_csv(roi_geometry_file_path)
print(f"Successfully loaded ROI geometry file: {roi_geometry_file_path}")
print("ROI CSV Head:")
print(roi_polygons_df.head())

# %% [markdown]
# # Process ROI geometries and perform spatial selection

# %%
# Convert WKT string geometries to Shapely geometry objects
roi_polygons_df['geometry'] = roi_polygons_df['geometry'].apply(wkt.loads)

# Create a GeoDataFrame from the ROI data
# Assume cell_boundaries_gdf is already defined and has a CRS.
# If cell_boundaries_gdf.crs is None, roi_gdf.crs will also be None,
# which is acceptable if coordinates are in the same arbitrary Cartesian system.
current_crs = None
if 'cell_boundaries_gdf' in locals() and cell_boundaries_gdf is not None and hasattr(cell_boundaries_gdf, 'crs'):
    current_crs = cell_boundaries_gdf.crs
    print(f"Using CRS from cell_boundaries_gdf: {current_crs}")
else:
    print("Warning: cell_boundaries_gdf not found or has no CRS. Assuming planar coordinates for ROIs.")

# %%
roi_gdf = gpd.GeoDataFrame(roi_polygons_df, geometry='geometry', crs=current_crs)
print("Successfully converted ROI geometries to GeoDataFrame.")
print("ROI GeoDataFrame Head:")
roi_gdf.head()

# %%
# Perform spatial selection of cells within ROIs
print(f"Shape of cell_boundaries_gdf before spatial join: {cell_boundaries_gdf.shape}")
print(f"Shape of roi_gdf before spatial join: {roi_gdf.shape}")

# %%
# Prepare the left GeoDataFrame for sjoin to avoid 'EntityID' column clash
# The original cell_boundaries_gdf has 'EntityID' as both index and column.
# Renaming the index ensures that when sjoin (or its internal functions)
# calls reset_index(), the new column from the index doesn't conflict.
cell_boundaries_gdf_sjoin_ready = cell_boundaries_gdf.copy()
cell_boundaries_gdf_sjoin_ready.index.name = 'original_cell_EntityID_idx' # Rename the index

# %%
# Spatial join: find cells whose geometries are 'within' the ROI polygons
# 'how="inner"' means only cells that are within an ROI are kept.
# 'predicate="within"' checks if cell geometry is entirely within ROI geometry.
# Added lsuffix and rsuffix to handle any potential column name overlaps clearly.
cells_in_rois_gdf = gpd.sjoin(
    cell_boundaries_gdf_sjoin_ready,
    roi_gdf,
    how="inner",
    predicate="within",
    lsuffix='_cell',
    rsuffix='_roi'
)

# %%
cells_in_rois_gdf.head()

# %%
print(f"\nFound {cells_in_rois_gdf.shape[0]} cells within the defined ROIs.")

# %%
print("Head of cells_in_rois_gdf (cells spatially selected by ROIs):")
print(cells_in_rois_gdf.head())

# %%
# Analyze the selected cells
print("\nCell counts per ROI:")
print(cells_in_rois_gdf.groupby('group').size())

# %%
# Visualize selected cells with ROIs

if plot_enable:

    print("\nCell counts per ROI group (from 'group' column of ROI CSV):")
    print(cells_in_rois_gdf.groupby('group').size())

    print("\nPlotting selected cells and ROIs...")
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))

    # Plot all original cell boundaries lightly
    cell_boundaries_gdf.plot(ax=ax, color='lightgray', edgecolor='silver', alpha=0.3, label='All Cells (Original)')

    unique_groups = cells_in_rois_gdf['group'].nunique()
    cells_in_rois_gdf.plot(ax=ax, column='group', legend=True, alpha=0.7, categorical=True,
                                           legend_kwds={'title': "ROI Group", 'loc': 'upper right', 'bbox_to_anchor': (1.25, 1)})

    # Plot ROI boundaries
    roi_gdf.plot(ax=ax, facecolor='none', edgecolor='blue', linewidth=2, label='ROI Polygons')

    ax.set_title("Cells within Defined ROIs")
    ax.set_xlabel("X-coordinate")
    ax.set_ylabel("Y-coordinate")

    # Adjust legend display
    handles, labels = ax.get_legend_handles_labels()

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()

# %% [markdown]
# # Differential Gene Expression (DGE) Analysis between ROIs
# 
# This section performs differential gene expression analysis between two selected ROIs.
# It assumes that `adata` (the AnnData object with gene expression) and
# `cells_in_rois_gdf` (GeoDataFrame with cells selected by ROIs and their group assignments)
# are available from the previous steps.

# %%
# DGE Analysis
print("\nStarting Differential Gene Expression Analysis...")
roi_group_column = 'group'

# Get unique ROI groups
available_roi_groups = cells_in_rois_gdf[roi_group_column].unique()
print(f"Available ROI groups for DGE: {available_roi_groups}")

# %%
roi_group_1_name = available_roi_groups[0]
roi_group_2_name = available_roi_groups[1]
print(f"Selected ROI groups for DGE: '{roi_group_1_name}' vs '{roi_group_2_name}'")

# %%
# Create a mapping of cell EntityID to its ROI group
cell_to_roi_map = cells_in_rois_gdf[roi_group_column].to_dict()
for key, value in list(cell_to_roi_map.items())[:5]:
    print(f"{key}: {value}")

# %%
# Add 'roi_assignment' to adata.obs
# Initialize with 'unassigned' or a suitable default
adata.obs['roi_assignment'] = 'unassigned'

# %%
print(f"Debug: adata.obs.index type: {adata.obs.index.dtype}, length: {len(adata.obs.index)}")
if len(adata.obs.index) > 0: print(f"Debug: First 3 adata.obs.index examples: {adata.obs.index[:3].tolist()}")

# %%
# Keys for cell_to_roi_map
cell_to_roi_map = {str(key): value for key, value in cell_to_roi_map.items()}
map_keys = list(cell_to_roi_map.keys())

# %%
# Perform intersection ensuring
adata_index_as_str = adata.obs.index
map_keys_as_str_for_intersection = pd.Index(map_keys)

valid_intersected_ids_str = adata_index_as_str.intersection(map_keys_as_str_for_intersection)
print(f"Debug: Number of common string IDs found by intersection: {len(valid_intersected_ids_str)}")
if len(valid_intersected_ids_str) > 0: print(f"Debug: Common string IDs example: {valid_intersected_ids_str[:3].tolist()}")

# %%
# Assign ROI based on the intersection.
for cell_id_str_from_intersection in valid_intersected_ids_str:
    if cell_id_str_from_intersection in adata.obs.index:
            adata.obs.loc[cell_id_str_from_intersection, 'roi_assignment'] = cell_to_roi_map[cell_id_str_from_intersection]

# %%
# Check how many cells were assigned
print(f"Cells assigned to ROIs in adata.obs: \n{adata.obs['roi_assignment'].value_counts()}")

# %%
# Filter adata for cells belonging to the two selected ROI groups
adata_dge = adata[adata.obs['roi_assignment'].isin([roi_group_1_name, roi_group_2_name])].copy()
print(f"Shape of AnnData object for DGE (cells from '{roi_group_1_name}' and '{roi_group_2_name}'): {adata_dge.shape}")

# %%
print(f"Performing DGE using Wilcoxon rank-sum test, comparing {roi_group_1_name} (reference) vs {roi_group_2_name}...")
# Ensure data is log-transformed for rank_genes_groups, as it expects log-normalized counts.
# This will also populate adata_dge.uns['log1p'] correctly.
print("Applying sc.pp.log1p(adata_dge) before rank_genes_groups...")
sc.pp.log1p(adata_dge)
sc.tl.rank_genes_groups(adata_dge, 
                        groupby='roi_assignment', 
                        groups=[roi_group_1_name], # Group to test
                        reference=roi_group_2_name, # Reference group
                        method='wilcoxon', # or 't-test', 'logreg'
                        corr_method='benjamini-hochberg', # multiple testing correction
                        key_added='rank_genes_dge') # Store results in adata_dge.uns['rank_genes_dge']

print("DGE calculation complete.")

# %%
# Visualize DGE results
print(f"Creating volcano plot for {roi_group_1_name}...")

# Create the volcano plot
fig, ax = create_volcano_plot(
    adata_dge, 
    group_name=roi_group_1_name,
    key='rank_genes_dge',
    pval_threshold=0.05,
    logfc_threshold=0.5,
    title=f"Volcano Plot: {roi_group_1_name} vs {roi_group_2_name}",
    show_gene_labels=True,
    top_genes=10
)

plt.show()

# %%
# Displaying the DGE results as a DataFrame
dge_results_df = pd.DataFrame(adata_dge.uns['rank_genes_dge']['names'])[roi_group_1_name]
dge_scores_df = pd.DataFrame(adata_dge.uns['rank_genes_dge']['scores'])[roi_group_1_name]
dge_pvals_df = pd.DataFrame(adata_dge.uns['rank_genes_dge']['pvals_adj'])[roi_group_1_name]
dge_logfc_df = pd.DataFrame(adata_dge.uns['rank_genes_dge']['logfoldchanges'])[roi_group_1_name]

summary_df = pd.DataFrame({
    'gene': dge_results_df,
    'score': dge_scores_df,
    'pval_adj': dge_pvals_df,
    'log2fc': dge_logfc_df
})
print(f"\nTop differentially expressed genes for {roi_group_1_name} compared to {roi_group_2_name}:")
print(summary_df.head(20))

# %%
print("\n--- End of Script ---")


