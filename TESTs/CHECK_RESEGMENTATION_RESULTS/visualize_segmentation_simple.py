#!/usr/bin/env python3
"""
Simplified visualization script for cell segmentation results.
Focuses on tiles with detected cells to avoid memory issues.
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any

import numpy as np
import pandas as pd
import matplotlib
# Set backend for headless environments (cluster/server)
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from skimage import io, exposure
import geopandas as gpd
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

class SimpleSegmentationVisualizer:
    """Simplified class to visualize segmentation results on tiles with cells."""
    
    def __init__(self, 
                 output_dir: str = "2task_cellpose2_p30-E165_R1_roi_analysis_parallel",
                 data_dir: str = "DATA/p30-E165/R1"):
        """Initialize the visualizer."""
        self.output_dir = Path(output_dir)
        self.data_dir = Path(data_dir)
        
        # Load segmentation specification
        self.spec = self._load_segmentation_spec()
        
        # Set up visualization parameters
        self.boundary_color = 'red'
        self.boundary_width = 1.5
        self.figure_size = (15, 5)
        self.contrast_percentiles = (1, 99)
        
    def _load_segmentation_spec(self) -> Dict[str, Any]:
        """Load segmentation specification."""
        spec_file = self.output_dir / "segmentation_specification_roi.json"
        try:
            with open(spec_file, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            print(f"Warning: Segmentation specification {spec_file} not found.")
            return {}
    
    def load_image(self, channel: str, z_layer: int = 6) -> np.ndarray:
        """Load an image for a specific channel and z-layer."""
        image_path = self.data_dir / "images" / f"mosaic_{channel}_z{z_layer}.tif"
        
        if not image_path.exists():
            raise FileNotFoundError(f"Image not found: {image_path}")
            
        image = io.imread(str(image_path))
        return image
    
    def load_tile_segmentation(self, tile_index: int) -> Optional[gpd.GeoDataFrame]:
        """Load segmentation results for a specific tile."""
        tile_file = self.output_dir / "result_tiles" / f"cell_{tile_index}.parquet"
        
        if not tile_file.exists():
            return None
            
        tile_boundaries = gpd.read_parquet(str(tile_file))
        return tile_boundaries
    
    def enhance_image(self, image: np.ndarray) -> np.ndarray:
        """Enhance image contrast for better visualization."""
        # Apply percentile-based contrast stretching
        p_low, p_high = np.percentile(image, self.contrast_percentiles)
        image_enhanced = exposure.rescale_intensity(image, in_range=(p_low, p_high))
        
        # Apply adaptive histogram equalization
        image_enhanced = exposure.equalize_adapthist(image_enhanced, clip_limit=0.01)
        
        return image_enhanced
    
    def get_tile_window(self, tile_index: int) -> Tuple[int, int, int, int]:
        """Get the window coordinates for a specific tile."""
        if 'window_grid' not in self.spec or 'windows' not in self.spec['window_grid']:
            raise ValueError("Window grid information not found in specification")
            
        windows = self.spec['window_grid']['windows']
        if tile_index >= len(windows):
            raise ValueError(f"Tile index {tile_index} out of range (max: {len(windows)-1})")
            
        return tuple(windows[tile_index])
    
    def crop_image_to_tile(self, image: np.ndarray, tile_index: int) -> np.ndarray:
        """Crop image to a specific tile window."""
        x_min, y_min, width, height = self.get_tile_window(tile_index)
        x_max = x_min + width
        y_max = y_min + height
        
        # Ensure coordinates are within image bounds
        x_min = max(0, x_min)
        y_min = max(0, y_min)
        x_max = min(image.shape[1], x_max)
        y_max = min(image.shape[0], y_max)
        
        return image[y_min:y_max, x_min:x_max]
    
    def visualize_tile_segmentation(self, tile_index: int, channels: List[str] = ['DAPI', 'PolyT'], 
                                  z_layer: int = 6, save_path: Optional[str] = None):
        """Visualize segmentation results for a specific tile."""
        print(f"Visualizing tile {tile_index}...")
        
        # Load tile segmentation
        tile_boundaries = self.load_tile_segmentation(tile_index)
        if tile_boundaries is None or len(tile_boundaries) == 0:
            print(f"No segmentation data found for tile {tile_index}")
            return
        
        # Get tile window
        x_min, y_min, width, height = self.get_tile_window(tile_index)
        
        # Create figure
        n_channels = len(channels)
        fig, axes = plt.subplots(1, n_channels, figsize=self.figure_size)
        if n_channels == 1:
            axes = [axes]
        
        # Load and display each channel
        for i, channel in enumerate(channels):
            try:
                # Load full image
                print(f"Loading {channel} channel...")
                image = self.load_image(channel, z_layer)
                
                # Crop to tile
                tile_image = self.crop_image_to_tile(image, tile_index)
                
                # Enhance contrast
                tile_image_enhanced = self.enhance_image(tile_image)
                
                # Display image
                axes[i].imshow(tile_image_enhanced, cmap='gray', aspect='equal')
                axes[i].set_title(f'{channel} - Tile {tile_index} ({len(tile_boundaries)} cells)')
                axes[i].set_xlabel('X (pixels)')
                axes[i].set_ylabel('Y (pixels)')
                
                # Overlay segmentation boundaries
                print(f"Overlaying {len(tile_boundaries)} cell boundaries...")
                
                # Transform boundaries from micron space to pixel space, then to tile coordinates
                boundaries_pixel = self._transform_boundaries_to_tile_coords(tile_boundaries, tile_index)
                
                for idx, row in boundaries_pixel.iterrows():
                    if hasattr(row, 'geometry') and hasattr(row.geometry, 'exterior'):
                        coords = np.array(row.geometry.exterior.coords)
                        axes[i].plot(coords[:, 0], coords[:, 1],
                                   color=self.boundary_color,
                                   linewidth=self.boundary_width,
                                   alpha=0.8)
                
            except Exception as e:
                print(f"Error processing channel {channel}: {e}")
                axes[i].text(0.5, 0.5, f'Error loading\n{channel}', 
                           ha='center', va='center', transform=axes[i].transAxes)
                axes[i].set_title(f'{channel} (Error)')
        
        plt.tight_layout()
        
        # Add information text
        info_text = f"Tile {tile_index}: {len(tile_boundaries)} cells detected"
        fig.suptitle(info_text, fontsize=14, y=0.98)
        
        # Always save the image, either to specified path or default location
        if save_path:
            output_path = save_path
        else:
            # Create default output directory
            default_output_dir = Path("segmentation_visualizations")
            default_output_dir.mkdir(exist_ok=True)
            output_path = default_output_dir / f"tile_{tile_index}_segmentation.png"
        
        plt.savefig(output_path, dpi=100, bbox_inches='tight')
        print(f"Saved visualization to {output_path}")
        
        # Close the figure to free memory
        plt.close(fig)
    
    def _transform_boundaries_to_tile_coords(self, tile_boundaries: gpd.GeoDataFrame, tile_index: int) -> gpd.GeoDataFrame:
        """
        Transform cell boundaries from micron coordinates to tile pixel coordinates.
        
        Args:
            tile_boundaries: GeoDataFrame with boundaries in micron coordinates
            tile_index: Index of the tile
            
        Returns:
            GeoDataFrame with boundaries in tile pixel coordinates
        """
        if len(tile_boundaries) == 0:
            return tile_boundaries
        
        # Get transformation matrix from spec
        if 'input_data' not in self.spec or 'micron_to_mosaic_tform' not in self.spec['input_data']:
            print("Warning: No micron to mosaic transform found in spec")
            return tile_boundaries
        
        # Get tile window coordinates
        x_min, y_min, width, height = self.get_tile_window(tile_index)
        
        # Extract transformation matrix
        T = self.spec['input_data']['micron_to_mosaic_tform']
        # Matrix for geopandas.GeoSeries.affine_transform is [a, b, d, e, xoff, yoff]
        # T = [[a, b, xoff], [d, e, yoff], [0, 0, 1]]
        micron_to_mosaic_matrix = [T[0][0], T[0][1], T[1][0], T[1][1], T[0][2], T[1][2]]
        
        print(f"Applying micron to mosaic transform: {micron_to_mosaic_matrix}")
        
        # Create a copy for transformation
        boundaries_transformed = tile_boundaries.copy()
        
        # Transform from micron to mosaic pixel coordinates
        boundaries_transformed.geometry = boundaries_transformed.geometry.affine_transform(micron_to_mosaic_matrix)
        
        # Translate to tile coordinates (relative to tile origin)
        print(f"Translating to tile coordinates: offset=({-x_min}, {-y_min})")
        boundaries_transformed.geometry = boundaries_transformed.geometry.translate(xoff=-x_min, yoff=-y_min)
        
        print(f"Transformed boundaries bounds: {boundaries_transformed.total_bounds}")
        
        return boundaries_transformed
    
    def get_tiles_with_cells(self) -> List[Tuple[int, int]]:
        """Get list of tiles that contain cells."""
        tiles_with_cells = []
        
        tile_files = list((self.output_dir / "result_tiles").glob("cell_*.parquet"))
        
        for tile_file in sorted(tile_files):
            tile_index = int(tile_file.stem.split('_')[1])
            
            try:
                tile_boundaries = gpd.read_parquet(str(tile_file))
                n_cells = len(tile_boundaries)
                
                if n_cells > 0:
                    tiles_with_cells.append((tile_index, n_cells))
                    
            except Exception as e:
                print(f"Error reading tile {tile_index}: {e}")
        
        return sorted(tiles_with_cells, key=lambda x: x[1], reverse=True)  # Sort by cell count
    
    def print_summary_report(self):
        """Print a formatted summary report."""
        print("\n" + "="*60)
        print("SEGMENTATION RESULTS SUMMARY")
        print("="*60)
        
        tiles_with_cells = self.get_tiles_with_cells()
        
        print(f"Tiles with cells: {len(tiles_with_cells)}")
        
        total_cells = sum(count for _, count in tiles_with_cells)
        print(f"Total cells detected: {total_cells}")
        
        print(f"\nTile Details (sorted by cell count):")
        for tile_index, n_cells in tiles_with_cells:
            print(f"  Tile {tile_index:2d}: {n_cells:4d} cells")
        
        print("="*60)


def main():
    """Main function to run the simplified visualization script."""
    parser = argparse.ArgumentParser(
        description="Visualize cell segmentation results (simplified version)"
    )
    parser.add_argument(
        "--output-dir", 
        default="2task_cellpose2_p30-E165_R1_roi_analysis_parallel",
        help="Directory containing segmentation results"
    )
    parser.add_argument(
        "--data-dir",
        default="DATA/p30-E165/R1",
        help="Directory containing original image data"
    )
    parser.add_argument(
        "--tile-index",
        type=int,
        help="Specific tile index to visualize (optional)"
    )
    parser.add_argument(
        "--channels",
        nargs="+",
        default=["DAPI", "PolyT"],
        help="Image channels to visualize"
    )
    parser.add_argument(
        "--z-layer",
        type=int,
        default=6,
        help="Z-layer to visualize"
    )
    parser.add_argument(
        "--save-dir",
        help="Directory to save visualization images (optional)"
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Only print summary report without visualizations"
    )
    parser.add_argument(
        "--all-tiles",
        action="store_true",
        help="Visualize all tiles with cells"
    )
    
    args = parser.parse_args()
    
    # Initialize visualizer
    try:
        visualizer = SimpleSegmentationVisualizer(
            output_dir=args.output_dir,
            data_dir=args.data_dir
        )
    except Exception as e:
        print(f"Error initializing visualizer: {e}")
        sys.exit(1)
    
    # Print summary report
    visualizer.print_summary_report()
    
    if args.summary_only:
        return
    
    # Create save directory if specified
    save_dir = None
    if args.save_dir:
        save_dir = Path(args.save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
    
    # Get tiles with cells
    tiles_with_cells = visualizer.get_tiles_with_cells()
    
    if not tiles_with_cells:
        print("No tiles with cells found.")
        return
    
    # Visualize specific tile, all tiles, or ask user
    if args.tile_index is not None:
        # Check if the specified tile has cells
        tile_cell_counts = {tile_idx: count for tile_idx, count in tiles_with_cells}
        if args.tile_index in tile_cell_counts:
            save_path = None
            if save_dir:
                save_path = save_dir / f"tile_{args.tile_index}_segmentation.png"
            
            visualizer.visualize_tile_segmentation(
                tile_index=args.tile_index,
                channels=args.channels,
                z_layer=args.z_layer,
                save_path=str(save_path) if save_path else None
            )
        else:
            print(f"Tile {args.tile_index} has no cells. Available tiles with cells: {[t[0] for t in tiles_with_cells]}")
    
    elif args.all_tiles:
        # Visualize all tiles with cells
        for tile_index, n_cells in tiles_with_cells:
            save_path = None
            if save_dir:
                save_path = save_dir / f"tile_{tile_index}_segmentation.png"
            
            visualizer.visualize_tile_segmentation(
                tile_index=tile_index,
                channels=args.channels,
                z_layer=args.z_layer,
                save_path=str(save_path) if save_path else None
            )
    
    else:
        # Interactive mode - ask user which tiles to visualize
        print(f"\nAvailable tiles with cells: {[f'{t[0]} ({t[1]} cells)' for t in tiles_with_cells]}")
        
        try:
            response = input("Enter tile index to visualize (or 'all' for all tiles): ").strip()
            
            if response.lower() == 'all':
                tiles_to_show = [t[0] for t in tiles_with_cells]
            else:
                tile_idx = int(response)
                tiles_to_show = [tile_idx] if tile_idx in [t[0] for t in tiles_with_cells] else []
            
            for tile_idx in tiles_to_show:
                save_path = None
                if save_dir:
                    save_path = save_dir / f"tile_{tile_idx}_segmentation.png"
                
                visualizer.visualize_tile_segmentation(
                    tile_index=tile_idx,
                    channels=args.channels,
                    z_layer=args.z_layer,
                    save_path=str(save_path) if save_path else None
                )
        
        except (KeyboardInterrupt, ValueError) as e:
            print(f"\nVisualization cancelled: {e}")


if __name__ == "__main__":
    main()