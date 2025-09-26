# %%
import os
import geopandas as gpd
import matplotlib.pyplot as plt
import rasterio
from rasterio.windows import Window
import json # Added for loading spec file
from shapely.affinity import affine_transform # Though GeoSeries.affine_transform is preferred

# Ensure the working directory is set if running this script directly and paths are relative
# For consistency with your original script, I'm keeping os.chdir, but be mindful if you move the script.
try:
    os.chdir("E:/Githubs/SPATIAL_data/Segmentation/Vpt-segmentation")
except FileNotFoundError:
    print("Warning: Could not change directory to E:/Githubs/SPATIAL_data/Segmentation/Vpt-segmentation. Make sure paths are correct.")

# --- Configuration ---
tile_index_to_check = 222
segmentation_spec_path = "E:/Githubs/SPATIAL_data/Segmentation/Vpt-segmentation/analysis_outputs_prepare_test/segmentation_specification.json"
parquet_file_path = f"E:/Githubs/SPATIAL_data/Segmentation/Vpt-segmentation/analysis_outputs_prepare_test/result_tiles/cell_{tile_index_to_check}.parquet"
# Path to the full mosaic DAPI image for the correct Z-plane (z3 in your case)
dapi_image_path = "E:/Githubs/SPATIAL_data/data_p30-E165/R1/images/mosaic_DAPI_z3.tif"


# --- Load Segmentation Specification for Transform Matrix and Tile Info ---
micron_to_mosaic_matrix_gp = None
tile_col_offset = 0
tile_row_offset = 0
tile_width = 2800 # Default, will be updated from spec if possible
tile_height = 2800 # Default

try:
    with open(segmentation_spec_path, 'r') as f:
        spec_data = json.load(f)
    
    # Extract micron to mosaic transform for geopandas
    T = spec_data['input_data']['micron_to_mosaic_tform']
    # Matrix for geopandas.GeoSeries.affine_transform is [a, b, d, e, xoff, yoff]
    # T = [[a, b, xoff], [d, e, yoff], [0, 0, 1]]
    micron_to_mosaic_matrix_gp = [T[0][0], T[0][1], T[1][0], T[1][1], T[0][2], T[1][2]]
    print(f"Micron to Mosaic transform matrix (for GeoPandas): {micron_to_mosaic_matrix_gp}")

    # Extract tile info for the specified tile_index_to_check
    # The 'windows' array in spec_data['window_grid']['windows'] is 0-indexed
    # Each window is [col_off, row_off, width, height]
    tile_info = spec_data['window_grid']['windows'][tile_index_to_check]
    tile_col_offset = tile_info[0]
    tile_row_offset = tile_info[1]
    tile_width = tile_info[2] # Use actual width from spec
    tile_height = tile_info[3] # Use actual height from spec
    print(f"Tile {tile_index_to_check} info: col_off={tile_col_offset}, row_off={tile_row_offset}, width={tile_width}, height={tile_height}")

except Exception as e:
    print(f"Error loading or parsing segmentation_specification.json: {e}")
    print("Proceeding with default/hardcoded tile offsets and no micron transform if matrix is None.")
    # Fallback to previously known values if spec load fails, though micron transform will be missing
    if micron_to_mosaic_matrix_gp is None: # Ensure we have some tile offsets if spec fails
        tile_col_offset = 0    # x_start from previous logs for tile 
        tile_row_offset = 20800 # y_start from previous logs for tile 
        # tile_width and tile_height remain default or from previous hardcoding

# Tile  coordinates and size (from vpt logs: Tile  [0, 20800, 2800, 2800])
# These are (col_off, row_off, width, height) for rasterio.windows.Window
# This is now dynamically loaded but kept here as a comment for reference
# tile_col_offset = 0
# tile_row_offset = 20800
# tile_width = 2800
# tile_height = 2800
# tile_col_offset = 0    # x_start # Overridden by spec file
# tile_row_offset = 20800 # y_start # Overridden by spec file
# tile_width = 2800 # Overridden by spec file
# tile_height = 2800 # Overridden by spec file

# --- Load Segmentation Polygons ---
try:
    gdf = gpd.read_parquet(parquet_file_path)
    print(f"Successfully loaded {len(gdf)} polygons from {parquet_file_path}")
    print("GeoDataFrame columns:", gdf.columns)
    print("GeoDataFrame head:\n", gdf.head())

    # --- Explicitly set geometry column if needed ---
    # Check if the active geometry column is correctly set.
    # If your geometry data is in a column NOT named 'geometry',
    # you might need to set it explicitly. For example, if it's in 'cell_boundary':
    # gdf = gdf.set_geometry('cell_boundary', inplace=False)
    # Or, if it's already the active geometry but has a different name,
    # gdf.geometry.name will give its name.
    
    # For now, we assume 'geometry' is the intended column or geopandas found it.
    # If KeyError persists, check gdf.columns output and adjust here.
    # The active geometry column is 'Geometry' (capital G)
    # No need to set it if geopandas found it, which it did.

except Exception as e:
    print(f"Error loading parquet file {parquet_file_path}: {e}")
    gdf = None

# --- Load Image Tile ---
image_tile = None
try:
    with rasterio.open(dapi_image_path) as src:
        # Define the window to read
        window = Window(col_off=tile_col_offset, row_off=tile_row_offset,
                        width=tile_width, height=tile_height)
        print(f"Reading window {window} from {dapi_image_path}")
        image_tile = src.read(1, window=window) # Read the first band
    print(f"Successfully loaded image tile. Shape: {image_tile.shape}, Min: {image_tile.min()}, Max: {image_tile.max()}")
except Exception as e:
    print(f"Error loading image tile from {dapi_image_path}: {e}")

# --- Plotting ---
if gdf is not None and image_tile is not None:
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    # Display the image tile
    # The extent should be from 0 to tile_width and 0 to tile_height for the cropped image
    ax.imshow(image_tile, cmap='gray', extent=(0, tile_width, tile_height, 0))
                                            # Note: extent y-axis is inverted for imshow to match typical image coordinates

    # Translate geometries to be relative to the tile's origin for plotting
    # The parquet files from vpt run-segmentation-on-tile usually store geometries
    # in global mosaic coordinates.
    # We need to shift them by the tile's offset.
    
    gdf_pixel_space = gdf.copy()

    if micron_to_mosaic_matrix_gp is not None:
        print("Applying micron to mosaic pixel space transformation...")
        # Ensure the active geometry column is used for transformation
        active_geom_col = gdf_pixel_space.geometry.name
        gdf_pixel_space[active_geom_col] = gdf_pixel_space[active_geom_col].affine_transform(micron_to_mosaic_matrix_gp)
        print(f"GeoDataFrame head after micron to pixel transform:\n{gdf_pixel_space.head()}")
    else:
        print("Skipping micron to mosaic transform as matrix was not loaded.")

    # Now translate these global pixel coordinates to be relative to the tile's origin
    print(f"Translating pixel coordinates by xoff={-tile_col_offset}, yoff={-tile_row_offset}")
    active_geom_col = gdf_pixel_space.geometry.name # Get active geom column name again (could have changed if set_geometry was used)
    gdf_translated_to_tile = gdf_pixel_space.copy() # Make a new copy for this step
    gdf_translated_to_tile[active_geom_col] = gdf_pixel_space[active_geom_col].translate(xoff=-tile_col_offset, yoff=-tile_row_offset)
    
    print(f"GeoDataFrame head after translating to tile origin:\n{gdf_translated_to_tile.head()}")
    print(f"Translated geometries total bounds: {gdf_translated_to_tile.total_bounds}")

    # Plot translated geometries.
    gdf_translated_to_tile.plot(ax=ax, facecolor='none', edgecolor='red', linewidth=0.5)
    
    ax.set_title(f"Tile {tile_index_to_check} Segmentation Overlay\n(Image: {os.path.basename(dapi_image_path)}, Tile Origin Pix: ({tile_col_offset},{tile_row_offset}))")
    # Set plot limits to match image dimensions and ensure y-axis is inverted like imshow
    ax.set_xlim(0, tile_width)
    ax.set_ylim(tile_height, 0) # Y-axis: 0 at top, tile_height at bottom

    ax.set_xlabel("X pixel (within tile)")
    ax.set_ylabel("Y pixel (within tile)")
    plt.tight_layout()
    plt.show()
elif gdf is None:
    print("Cannot plot because polygon data (gdf) could not be loaded.")
elif image_tile is None:
    print("Cannot plot because image tile data could not be loaded.")

# %%
