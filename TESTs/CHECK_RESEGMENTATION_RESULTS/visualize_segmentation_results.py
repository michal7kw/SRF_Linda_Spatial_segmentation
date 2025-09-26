#!/usr/bin/env python3
"""
Visualization script for cell segmentation results from VPT pipeline.

This script loads segmentation results and overlays them on original images
for visual quality assessment of the segmentation pipeline.

Usage:
    python visualize_segmentation_results.py [--output-dir OUTPUT_DIR] [--tile-index TILE_INDEX]

Author: Generated for SRF Spatial Segmentation Pipeline
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
from matplotlib.colors import ListedColormap
import seaborn as sns
from skimage import io, exposure
from skimage.transform import resize
import geopandas as gpd
from shapely.geometry import Polygon
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

class SegmentationVisualizer:
    """Class to visualize segmentation results on original images."""
    
    def __init__(self, 
                 output_dir: str = "2task_cellpose2_p30-E165_R1_roi_analysis_parallel",
                 data_dir: str = "DATA/p30-E165/R1",
                 roi_coords_file: str = "roi_coords.json"):
        """
        Initialize the visualizer.
        
        Args:
            output_dir: Directory containing segmentation results
            data_dir: Directory containing original image data
            roi_coords_file: JSON file with ROI coordinates
        """
        self.output_dir = Path(output_dir)
        self.data_dir = Path(data_dir)
        self.roi_coords_file = Path(roi_coords_file)
        
        # Load ROI coordinates
        self.roi_coords = self._load_roi_coords()
        
        # Load segmentation specification
        self.spec = self._load_segmentation_spec()
        
        # Set up visualization parameters
        self.setup_visualization_params()
        
    def _load_roi_coords(self) -> Dict[str, float]:
        """Load ROI coordinates from JSON file."""
        try:
            with open(self.roi_coords_file, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            print(f"Warning: ROI coordinates file {self.roi_coords_file} not found.")
            return {}
    
    def _load_segmentation_spec(self) -> Dict[str, Any]:
        """Load segmentation specification."""
        spec_file = self.output_dir / "segmentation_specification_roi.json"
        try:
            with open(spec_file, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            print(f"Warning: Segmentation specification {spec_file} not found.")
            return {}
    
    def setup_visualization_params(self):
        """Set up visualization parameters."""
        # Color schemes
        from matplotlib import cm
        self.cell_colors = cm.get_cmap('Set3')(np.linspace(0, 1, 12))
        self.boundary_color = 'red'
        self.boundary_width = 1.5
        
        # Figure parameters
        self.figure_size = (15, 10)
        self.dpi = 100
        
        # Image enhancement parameters
        self.contrast_percentiles = (1, 99)
        
    def load_image(self, channel: str, z_layer: int = 6) -> np.ndarray:
        """
        Load an image for a specific channel and z-layer.
        
        Args:
            channel: Image channel (e.g., 'DAPI', 'PolyT')
            z_layer: Z-layer index
            
        Returns:
            Image array
        """
        image_path = self.data_dir / "images" / f"mosaic_{channel}_z{z_layer}.tif"
        
        if not image_path.exists():
            raise FileNotFoundError(f"Image not found: {image_path}")
            
        image = io.imread(str(image_path))
        return image
    
    def load_segmentation_boundaries(self, space: str = "mosaic") -> gpd.GeoDataFrame:
        """
        Load segmentation boundaries from parquet file.
        
        Args:
            space: Coordinate space ('mosaic' or 'micron')
            
        Returns:
            GeoDataFrame with cell boundaries
        """
        boundary_file = self.output_dir / f"cellpose2_{space}_space.parquet"
        
        if not boundary_file.exists():
            raise FileNotFoundError(f"Boundary file not found: {boundary_file}")
            
        boundaries = gpd.read_parquet(str(boundary_file))
        return boundaries
    
    def load_tile_segmentation(self, tile_index: int) -> Optional[gpd.GeoDataFrame]:
        """
        Load segmentation results for a specific tile.
        
        Args:
            tile_index: Index of the tile
            
        Returns:
            GeoDataFrame with cell boundaries for the tile, or None if not found
        """
        tile_file = self.output_dir / "result_tiles" / f"cell_{tile_index}.parquet"
        
        if not tile_file.exists():
            print(f"Warning: Tile file not found: {tile_file}")
            return None
            
        tile_boundaries = gpd.read_parquet(str(tile_file))
        return tile_boundaries
    
    def enhance_image(self, image: np.ndarray) -> np.ndarray:
        """
        Enhance image contrast for better visualization.
        
        Args:
            image: Input image array
            
        Returns:
            Enhanced image array
        """
        # Apply percentile-based contrast stretching
        p_low, p_high = np.percentile(image, self.contrast_percentiles)
        image_enhanced = exposure.rescale_intensity(image, in_range=(p_low, p_high))
        
        # Apply adaptive histogram equalization
        image_enhanced = exposure.equalize_adapthist(image_enhanced, clip_limit=0.01)
        
        return image_enhanced
    
    def get_tile_window(self, tile_index: int) -> Tuple[int, int, int, int]:
        """
        Get the window coordinates for a specific tile.
        
        Args:
            tile_index: Index of the tile
            
        Returns:
            Tuple of (x_min, y_min, width, height) in mosaic coordinates
        """
        if 'window_grid' not in self.spec or 'windows' not in self.spec['window_grid']:
            raise ValueError("Window grid information not found in specification")
            
        windows = self.spec['window_grid']['windows']
        if tile_index >= len(windows):
            raise ValueError(f"Tile index {tile_index} out of range (max: {len(windows)-1})")
            
        return tuple(windows[tile_index])
    
    def crop_image_to_tile(self, image: np.ndarray, tile_index: int) -> np.ndarray:
        """
        Crop image to a specific tile window.
        
        Args:
            image: Full mosaic image
            tile_index: Index of the tile
            
        Returns:
            Cropped image for the tile
        """
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
        """
        Visualize segmentation results for a specific tile.
        
        Args:
            tile_index: Index of the tile to visualize
            channels: List of image channels to display
            z_layer: Z-layer to visualize
            save_path: Optional path to save the figure
        """
        print(f"Visualizing tile {tile_index}...")
        
        # Load tile segmentation
        tile_boundaries = self.load_tile_segmentation(tile_index)
        if tile_boundaries is None:
            print(f"No segmentation data found for tile {tile_index}")
            return
        
        # Get tile window
        x_min, y_min, width, height = self.get_tile_window(tile_index)
        
        # Create figure
        n_channels = len(channels)
        fig, axes = plt.subplots(1, n_channels + 1, figsize=(5 * (n_channels + 1), 5))
        if n_channels == 0:
            axes = [axes]
        
        # Load and display each channel
        for i, channel in enumerate(channels):
            try:
                # Load full image
                image = self.load_image(channel, z_layer)
                
                # Crop to tile
                tile_image = self.crop_image_to_tile(image, tile_index)
                
                # Enhance contrast
                tile_image_enhanced = self.enhance_image(tile_image)
                
                # Display image
                axes[i].imshow(tile_image_enhanced, cmap='gray', aspect='equal')
                axes[i].set_title(f'{channel} (Tile {tile_index})')
                axes[i].set_xlabel('X (pixels)')
                axes[i].set_ylabel('Y (pixels)')
                
                # Overlay segmentation boundaries
                if tile_boundaries is not None and len(tile_boundaries) > 0:
                    for idx, row in tile_boundaries.iterrows():
                        if hasattr(row, 'geometry') and hasattr(row.geometry, 'exterior'):
                            # Convert coordinates relative to tile
                            coords = np.array(row.geometry.exterior.coords)
                            coords[:, 0] -= x_min  # Adjust x coordinates
                            coords[:, 1] -= y_min  # Adjust y coordinates
                            
                            axes[i].plot(coords[:, 0], coords[:, 1],
                                       color=self.boundary_color,
                                       linewidth=self.boundary_width,
                                       alpha=0.8)
                
            except Exception as e:
                print(f"Error processing channel {channel}: {e}")
                axes[i].text(0.5, 0.5, f'Error loading\n{channel}', 
                           ha='center', va='center', transform=axes[i].transAxes)
                axes[i].set_title(f'{channel} (Error)')
        
        # Create overlay visualization
        if len(channels) >= 2:
            try:
                # Load images for overlay
                dapi_image = self.load_image('DAPI', z_layer)
                polyt_image = self.load_image('PolyT', z_layer)
                
                # Crop to tile
                dapi_tile = self.crop_image_to_tile(dapi_image, tile_index)
                polyt_tile = self.crop_image_to_tile(polyt_image, tile_index)
                
                # Enhance contrast
                dapi_enhanced = self.enhance_image(dapi_tile)
                polyt_enhanced = self.enhance_image(polyt_tile)
                
                # Create RGB overlay
                overlay = np.zeros((*dapi_enhanced.shape, 3))
                overlay[:, :, 0] = polyt_enhanced  # Red channel
                overlay[:, :, 1] = dapi_enhanced   # Green channel
                overlay[:, :, 2] = dapi_enhanced   # Blue channel
                
                axes[-1].imshow(overlay, aspect='equal')
                axes[-1].set_title(f'Overlay with Segmentation (Tile {tile_index})')
                axes[-1].set_xlabel('X (pixels)')
                axes[-1].set_ylabel('Y (pixels)')
                
                # Overlay segmentation boundaries
                if tile_boundaries is not None and len(tile_boundaries) > 0:
                    for idx, row in tile_boundaries.iterrows():
                        if hasattr(row, 'geometry') and hasattr(row.geometry, 'exterior'):
                            coords = np.array(row.geometry.exterior.coords)
                            coords[:, 0] -= x_min
                            coords[:, 1] -= y_min
                            
                            axes[-1].plot(coords[:, 0], coords[:, 1],
                                        color='yellow',
                                        linewidth=self.boundary_width,
                                        alpha=0.9)
                
            except Exception as e:
                print(f"Error creating overlay: {e}")
                axes[-1].text(0.5, 0.5, 'Error creating\noverlay', 
                            ha='center', va='center', transform=axes[-1].transAxes)
                axes[-1].set_title('Overlay (Error)')
        
        plt.tight_layout()
        
        # Add information text
        info_text = f"Tile {tile_index}: {len(tile_boundaries) if tile_boundaries is not None else 0} cells detected"
        fig.suptitle(info_text, fontsize=14, y=0.98)
        
        # Always save the image, either to specified path or default location
        if save_path:
            output_path = save_path
        else:
            # Create default output directory
            default_output_dir = Path("segmentation_visualizations")
            default_output_dir.mkdir(exist_ok=True)
            output_path = default_output_dir / f"tile_{tile_index}_segmentation.png"
        
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        print(f"Saved visualization to {output_path}")
        
        # Close the figure to free memory
        plt.close(fig)
    
    def visualize_roi_overview(self, channels: List[str] = ['DAPI', 'PolyT'], 
                              z_layer: int = 6, save_path: Optional[str] = None):
        """
        Visualize the entire ROI with segmentation boundaries.
        
        Args:
            channels: List of image channels to display
            z_layer: Z-layer to visualize
            save_path: Optional path to save the figure
        """
        print("Creating ROI overview...")
        
        # Load full segmentation boundaries
        try:
            boundaries = self.load_segmentation_boundaries("mosaic")
        except FileNotFoundError:
            print("No compiled segmentation boundaries found.")
            return
        
        # Create figure
        n_channels = len(channels)
        fig, axes = plt.subplots(1, n_channels + 1, figsize=(6 * (n_channels + 1), 6))
        if n_channels == 0:
            axes = [axes]
        
        # Get ROI bounds for cropping
        if self.roi_coords:
            # Convert micron coordinates to mosaic coordinates using transform
            if 'input_data' in self.spec and 'micron_to_mosaic_tform' in self.spec['input_data']:
                transform = np.array(self.spec['input_data']['micron_to_mosaic_tform'])
                
                # Transform ROI coordinates
                roi_corners = np.array([
                    [self.roi_coords['x_min'], self.roi_coords['y_min'], 1],
                    [self.roi_coords['x_max'], self.roi_coords['y_max'], 1]
                ])
                
                mosaic_corners = roi_corners @ transform.T
                x_min_mosaic = int(mosaic_corners[0, 0])
                y_min_mosaic = int(mosaic_corners[0, 1])
                x_max_mosaic = int(mosaic_corners[1, 0])
                y_max_mosaic = int(mosaic_corners[1, 1])
            else:
                # Fallback: use tile windows to estimate ROI
                if 'window_grid' in self.spec and 'windows' in self.spec['window_grid']:
                    windows = np.array(self.spec['window_grid']['windows'])
                    x_min_mosaic = int(windows[:, 0].min())
                    y_min_mosaic = int(windows[:, 1].min())
                    x_max_mosaic = int((windows[:, 0] + windows[:, 2]).max())
                    y_max_mosaic = int((windows[:, 1] + windows[:, 3]).max())
                else:
                    x_min_mosaic = y_min_mosaic = 0
                    x_max_mosaic = y_max_mosaic = -1
        else:
            x_min_mosaic = y_min_mosaic = 0
            x_max_mosaic = y_max_mosaic = -1
        
        # Load and display each channel
        for i, channel in enumerate(channels):
            try:
                # Load full image
                image = self.load_image(channel, z_layer)
                
                # Crop to ROI if coordinates are available
                if x_max_mosaic > x_min_mosaic and y_max_mosaic > y_min_mosaic:
                    roi_image = image[y_min_mosaic:y_max_mosaic, x_min_mosaic:x_max_mosaic]
                else:
                    roi_image = image
                    x_min_mosaic = y_min_mosaic = 0
                
                # Enhance contrast
                roi_image_enhanced = self.enhance_image(roi_image)
                
                # Display image
                axes[i].imshow(roi_image_enhanced, cmap='gray', aspect='equal')
                axes[i].set_title(f'{channel} (ROI Overview)')
                axes[i].set_xlabel('X (pixels)')
                axes[i].set_ylabel('Y (pixels)')
                
                # Overlay segmentation boundaries
                if len(boundaries) > 0:
                    for idx, row in boundaries.iterrows():
                        if hasattr(row, 'geometry') and hasattr(row.geometry, 'exterior'):
                            coords = np.array(row.geometry.exterior.coords)
                            # Adjust coordinates relative to ROI crop
                            coords[:, 0] -= x_min_mosaic
                            coords[:, 1] -= y_min_mosaic
                            
                            # Only plot if within the cropped region
                            if (coords[:, 0].min() >= 0 and coords[:, 0].max() < roi_image.shape[1] and
                                coords[:, 1].min() >= 0 and coords[:, 1].max() < roi_image.shape[0]):
                                axes[i].plot(coords[:, 0], coords[:, 1],
                                           color=self.boundary_color,
                                           linewidth=0.8,
                                           alpha=0.7)
                
            except Exception as e:
                print(f"Error processing channel {channel}: {e}")
                axes[i].text(0.5, 0.5, f'Error loading\n{channel}', 
                           ha='center', va='center', transform=axes[i].transAxes)
                axes[i].set_title(f'{channel} (Error)')
        
        # Create overlay visualization
        if len(channels) >= 2:
            try:
                # Load images for overlay
                dapi_image = self.load_image('DAPI', z_layer)
                polyt_image = self.load_image('PolyT', z_layer)
                
                # Crop to ROI
                if x_max_mosaic > x_min_mosaic and y_max_mosaic > y_min_mosaic:
                    dapi_roi = dapi_image[y_min_mosaic:y_max_mosaic, x_min_mosaic:x_max_mosaic]
                    polyt_roi = polyt_image[y_min_mosaic:y_max_mosaic, x_min_mosaic:x_max_mosaic]
                else:
                    dapi_roi = dapi_image
                    polyt_roi = polyt_image
                
                # Enhance contrast
                dapi_enhanced = self.enhance_image(dapi_roi)
                polyt_enhanced = self.enhance_image(polyt_roi)
                
                # Create RGB overlay
                overlay = np.zeros((*dapi_enhanced.shape, 3))
                overlay[:, :, 0] = polyt_enhanced  # Red channel
                overlay[:, :, 1] = dapi_enhanced   # Green channel
                overlay[:, :, 2] = dapi_enhanced   # Blue channel
                
                axes[-1].imshow(overlay, aspect='equal')
                axes[-1].set_title(f'ROI Overlay with Segmentation')
                axes[-1].set_xlabel('X (pixels)')
                axes[-1].set_ylabel('Y (pixels)')
                
                # Overlay segmentation boundaries
                if len(boundaries) > 0:
                    for idx, row in boundaries.iterrows():
                        if hasattr(row, 'geometry') and hasattr(row.geometry, 'exterior'):
                            coords = np.array(row.geometry.exterior.coords)
                            coords[:, 0] -= x_min_mosaic
                            coords[:, 1] -= y_min_mosaic
                            
                            if (coords[:, 0].min() >= 0 and coords[:, 0].max() < dapi_roi.shape[1] and
                                coords[:, 1].min() >= 0 and coords[:, 1].max() < dapi_roi.shape[0]):
                                axes[-1].plot(coords[:, 0], coords[:, 1],
                                            color='yellow',
                                            linewidth=0.8,
                                            alpha=0.8)
                
            except Exception as e:
                print(f"Error creating overlay: {e}")
                axes[-1].text(0.5, 0.5, 'Error creating\noverlay', 
                            ha='center', va='center', transform=axes[-1].transAxes)
                axes[-1].set_title('Overlay (Error)')
        
        plt.tight_layout()
        
        # Add information text
        info_text = f"ROI Overview: {len(boundaries)} cells detected"
        fig.suptitle(info_text, fontsize=16, y=0.98)
        
        # Always save the image, either to specified path or default location
        if save_path:
            output_path = save_path
        else:
            # Create default output directory
            default_output_dir = Path("segmentation_visualizations")
            default_output_dir.mkdir(exist_ok=True)
            output_path = default_output_dir / "roi_overview_segmentation.png"
        
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        print(f"Saved ROI overview to {output_path}")
        
        # Close the figure to free memory
        plt.close(fig)
    
    def generate_summary_report(self) -> Dict[str, Any]:
        """
        Generate a summary report of the segmentation results.
        
        Returns:
            Dictionary containing summary statistics
        """
        print("Generating summary report...")
        
        summary = {
            'roi_coordinates': self.roi_coords,
            'total_tiles': 0,
            'tiles_with_cells': 0,
            'total_cells': 0,
            'cells_per_tile': [],
            'tile_details': []
        }
        
        # Check compiled results
        try:
            boundaries = self.load_segmentation_boundaries("mosaic")
            summary['total_cells'] = len(boundaries)
        except FileNotFoundError:
            print("No compiled boundaries found.")
        
        # Check individual tiles
        if 'window_grid' in self.spec:
            summary['total_tiles'] = self.spec['window_grid'].get('num_tiles', 0)
        
        tile_files = list((self.output_dir / "result_tiles").glob("cell_*.parquet"))
        
        for tile_file in sorted(tile_files):
            tile_index = int(tile_file.stem.split('_')[1])
            
            try:
                tile_boundaries = gpd.read_parquet(str(tile_file))
                n_cells = len(tile_boundaries)
                
                if n_cells > 0:
                    summary['tiles_with_cells'] += 1
                
                summary['cells_per_tile'].append(n_cells)
                summary['tile_details'].append({
                    'tile_index': tile_index,
                    'n_cells': n_cells,
                    'file_path': str(tile_file)
                })
                
            except Exception as e:
                print(f"Error reading tile {tile_index}: {e}")
        
        return summary
    
    def print_summary_report(self):
        """Print a formatted summary report."""
        summary = self.generate_summary_report()
        
        print("\n" + "="*60)
        print("SEGMENTATION RESULTS SUMMARY")
        print("="*60)
        
        if summary['roi_coordinates']:
            print(f"ROI Coordinates:")
            print(f"  X: {summary['roi_coordinates']['x_min']:.1f} - {summary['roi_coordinates']['x_max']:.1f}")
            print(f"  Y: {summary['roi_coordinates']['y_min']:.1f} - {summary['roi_coordinates']['y_max']:.1f}")
        
        print(f"\nTile Information:")
        print(f"  Total tiles: {summary['total_tiles']}")
        print(f"  Tiles with cells: {summary['tiles_with_cells']}")
        
        print(f"\nCell Counts:")
        print(f"  Total cells detected: {summary['total_cells']}")
        
        if summary['cells_per_tile']:
            cells_per_tile = np.array(summary['cells_per_tile'])
            print(f"  Cells per tile - Mean: {cells_per_tile.mean():.1f}, "
                  f"Std: {cells_per_tile.std():.1f}")
            print(f"  Cells per tile - Min: {cells_per_tile.min()}, "
                  f"Max: {cells_per_tile.max()}")
        
        print(f"\nTile Details:")
        for tile_info in summary['tile_details']:
            print(f"  Tile {tile_info['tile_index']:2d}: {tile_info['n_cells']:4d} cells")
        
        print("="*60)


def main():
    """Main function to run the visualization script."""
    parser = argparse.ArgumentParser(
        description="Visualize cell segmentation results from VPT pipeline"
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
        "--roi-coords",
        default="roi_coords.json",
        help="JSON file with ROI coordinates"
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
    
    args = parser.parse_args()
    
    # Initialize visualizer
    try:
        visualizer = SegmentationVisualizer(
            output_dir=args.output_dir,
            data_dir=args.data_dir,
            roi_coords_file=args.roi_coords
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
    
    # Visualize specific tile or all tiles
    if args.tile_index is not None:
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
        # Show ROI overview
        save_path = None
        if save_dir:
            save_path = save_dir / "roi_overview_segmentation.png"
        
        visualizer.visualize_roi_overview(
            channels=args.channels,
            z_layer=args.z_layer,
            save_path=str(save_path) if save_path else None
        )
        
        # Ask user if they want to see individual tiles
        try:
            response = input("\nWould you like to visualize individual tiles? (y/n): ").lower()
            if response.startswith('y'):
                # Get available tiles
                tile_files = list((Path(args.output_dir) / "result_tiles").glob("cell_*.parquet"))
                tile_indices = sorted([int(f.stem.split('_')[1]) for f in tile_files])
                
                print(f"Available tiles: {tile_indices}")
                tile_input = input("Enter tile indices to visualize (comma-separated, or 'all'): ").strip()
                
                if tile_input.lower() == 'all':
                    tiles_to_show = tile_indices
                else:
                    tiles_to_show = [int(x.strip()) for x in tile_input.split(',') if x.strip().isdigit()]
                
                for tile_idx in tiles_to_show:
                    if tile_idx in tile_indices:
                        save_path = None
                        if save_dir:
                            save_path = save_dir / f"tile_{tile_idx}_segmentation.png"
                        
                        visualizer.visualize_tile_segmentation(
                            tile_index=tile_idx,
                            channels=args.channels,
                            z_layer=args.z_layer,
                            save_path=str(save_path) if save_path else None
                        )
                    else:
                        print(f"Tile {tile_idx} not found in results.")
        
        except KeyboardInterrupt:
            print("\nVisualization interrupted by user.")
        except Exception as e:
            print(f"Error during interactive visualization: {e}")


if __name__ == "__main__":
    main()