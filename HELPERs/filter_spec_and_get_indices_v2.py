import json
import argparse
import copy
import numpy as np
import pandas as pd

def load_transform_matrix(path: str) -> np.ndarray:
    """Loads the micron-to-mosaic transformation matrix."""
    # The matrix is space-delimited, not comma-delimited.
    return pd.read_csv(path, header=None, sep='\s+').values

def transform_coords(matrix: np.ndarray, x: float, y: float) -> tuple[float, float]:
    """Applies the transformation matrix to a coordinate."""
    point = np.array([x, y, 1])
    transformed_point = matrix.dot(point)
    return transformed_point[0], transformed_point[1]

def filter_spec_and_get_indices(spec_path: str, roi_path: str, transform_path: str, output_spec_path: str):
    """
    Filters a segmentation specification file to include only tiles that
    overlap with a given region of interest (ROI). It transforms the tile
    coordinates from pixel-space to micron-space for comparison.
    """
    with open(spec_path, 'r') as f:
        spec_data = json.load(f)

    with open(roi_path, 'r') as f:
        roi_coords = json.load(f)

    # Load the micron-to-pixel transformation matrix and calculate its inverse
    micron_to_pixel_matrix = load_transform_matrix(transform_path)
    pixel_to_micron_matrix = np.linalg.inv(micron_to_pixel_matrix)

    # ROI coordinates are already in micron space
    roi_x_min_micron = roi_coords['x_min']
    roi_y_min_micron = roi_coords['y_min']
    roi_x_max_micron = roi_coords['x_max']
    roi_y_max_micron = roi_coords['y_max']

    filtered_spec = copy.deepcopy(spec_data)
    filtered_windows = []
    original_indices = []

    for i, window in enumerate(spec_data['window_grid']['windows']):
        # Tile window coordinates are in pixel space
        tile_x_min_pixel = window[0]
        tile_y_min_pixel = window[1]
        tile_x_max_pixel = tile_x_min_pixel + window[2]
        tile_y_max_pixel = tile_y_min_pixel + window[3]

        # Transform tile corners to micron space for comparison
        tile_x_min_micron, tile_y_min_micron = transform_coords(pixel_to_micron_matrix, tile_x_min_pixel, tile_y_min_pixel)
        tile_x_max_micron, tile_y_max_micron = transform_coords(pixel_to_micron_matrix, tile_x_max_pixel, tile_y_max_pixel)

        # Check for overlap in micron space
        if (tile_x_min_micron < roi_x_max_micron and tile_x_max_micron > roi_x_min_micron and
                tile_y_min_micron < roi_y_max_micron and tile_y_max_micron > roi_y_min_micron):
            filtered_windows.append(window)
            original_indices.append(i)

    filtered_spec['window_grid']['windows'] = filtered_windows
    filtered_spec['window_grid']['num_tiles'] = len(filtered_windows)

    with open(output_spec_path, 'w') as f:
        json.dump(filtered_spec, f, indent=4)

    # Print original indices to stdout
    for index in original_indices:
        print(index)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter a spec file and get original tile indices for an ROI.")
    parser.add_argument('--input-spec', type=str, required=True,
                        help='Path to the full segmentation_specification.json file.')
    parser.add_argument('--input-roi', type=str, required=True,
                        help='Path to the roi_coords.json file.')
    parser.add_argument('--transform', type=str, required=True,
                        help='Path to the micron_to_mosaic_pixel_transform.csv file.')
    parser.add_argument('--output-spec', type=str, required=True,
                        help='Path to save the filtered segmentation_specification.json file.')

    args = parser.parse_args()

    filter_spec_and_get_indices(args.input_spec, args.input_roi, args.transform, args.output_spec)