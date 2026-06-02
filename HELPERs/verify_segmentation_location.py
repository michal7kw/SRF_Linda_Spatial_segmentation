import argparse
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import json

def plot_boundaries(boundaries_path: str, roi_path: str, output_plot_path: str):
    """
    Reads a parquet file with cell boundaries and plots them along with the ROI box.

    Args:
        boundaries_path: Path to the cellpose_micron_space.parquet file.
        roi_path: Path to the roi_coords.json file.
        output_plot_path: Path to save the output plot image.
    """
    try:
        # Read the cell boundaries
        gdf = gpd.read_parquet(boundaries_path)
    except Exception as e:
        print(f"Error reading parquet file: {e}")
        print("It's possible the compile step failed and this file was not created.")
        return

    # Read the ROI coordinates
    with open(roi_path, 'r') as f:
        roi_coords = json.load(f)

    # Create a figure and axes
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot the cell boundaries
    gdf.plot(ax=ax, edgecolor='blue', facecolor='none')

    # Create a rectangle for the ROI
    roi_rect = Rectangle(
        (roi_coords['x_min'], roi_coords['y_min']),
        roi_coords['x_max'] - roi_coords['x_min'],
        roi_coords['y_max'] - roi_coords['y_min'],
        edgecolor='red',
        facecolor='none',
        linewidth=2,
        label='Defined ROI'
    )
    ax.add_patch(roi_rect)

    # Set plot limits and labels
    ax.set_title('Verification of Segmented Cell Locations')
    ax.set_xlabel('X-coordinate (microns)')
    ax.set_ylabel('Y-coordinate (microns)')
    ax.legend()
    ax.set_aspect('equal', adjustable='box')
    plt.grid(True)

    # Save the plot
    plt.savefig(output_plot_path)
    print(f"Verification plot saved to: {output_plot_path}")
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Verify the location of segmented cell boundaries.")
    parser.add_argument('--input-boundaries', type=str, required=True,
                        help='Path to the cellpose_micron_space.parquet file from your ROI analysis.')
    parser.add_argument('--input-roi', type=str, required=True,
                        help='Path to the roi_coords.json file.')
    parser.add_argument('--output-plot', type=str, required=True,
                        help='Path to save the output verification plot.')

    args = parser.parse_args()

    plot_boundaries(args.input_boundaries, args.input_roi, args.output_plot)