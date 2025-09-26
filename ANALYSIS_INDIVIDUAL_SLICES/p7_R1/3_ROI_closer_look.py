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
# import squidpy
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

base_path = 'DATA/p0-p7/R1'
h5ad_file = os.path.join(base_path, 'data.h5ad')

# %%
roi_csv_file_name = 'p0_R1_ROI_28-07-25_11-57_geometry.csv'
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
# Cell boundaries
cell_boundaries_file = os.path.join(base_path, 'cell_boundaries.parquet')
cell_boundaries_gdf = None

# %%
cell_boundaries_gdf = gpd.read_parquet(cell_boundaries_file)
print(f"Loaded {cell_boundaries_file}. Shape: {cell_boundaries_gdf.shape}")

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
# # 2. Analyze transcript expression directly within ROIs

# %%
# Load transcript data
# The transcript file is typically named 'detected_transcripts.csv' in Vizgen datasets.
transcripts_file = os.path.join(base_path, 'detected_transcripts.csv')
if os.path.exists(transcripts_file):
    print(f"Loading transcripts from {transcripts_file}...")
    transcripts_df = pd.read_csv(transcripts_file)
    print(f"Loaded {len(transcripts_df)} transcripts.")
    
    # Convert transcripts DataFrame to a GeoDataFrame
    # MERSCOPE data uses 'global_x' and 'global_y' for coordinates.
    transcripts_gdf = gpd.GeoDataFrame(
        transcripts_df,
        geometry=gpd.points_from_xy(transcripts_df['global_x'], transcripts_df['global_y']),
        crs=current_crs  # Use the same CRS as the cell boundaries
    )
    print("Successfully converted transcripts to GeoDataFrame.")

    # Spatially join transcripts with ROIs
    transcripts_in_rois_gdf = gpd.sjoin(transcripts_gdf, roi_gdf, how="inner", predicate="within")
    
    print(f"\nFound {len(transcripts_in_rois_gdf)} transcripts within all ROIs.")

    # Analyze and plot gene expression for each ROI group
    for group_name, group_df in transcripts_in_rois_gdf.groupby('group'):
        print(f"\n--- Analysis for ROI Group: {group_name} ---")
        
        # Get gene counts
        gene_counts = group_df['gene'].value_counts()
        
        print(f"Total transcripts in this ROI group: {len(group_df)}")
        print("Top 20 most expressed genes:")
        print(gene_counts.head(20))
        
        # Plotting
        if plot_enable:
            plt.figure(figsize=(12, 8))
            sns.barplot(x=gene_counts.head(20).values, y=gene_counts.head(20).index, palette="viridis")
            plt.title(f"Top 20 Gene Expression in ROI Group '{group_name}' (Transcript Count)")
            plt.xlabel("Transcript Count")
            plt.ylabel("Gene")
            plt.tight_layout()
            plt.show()
else:
    print(f"Transcript file not found at: {transcripts_file}")
    print("Skipping transcript-level analysis.")

# %%
