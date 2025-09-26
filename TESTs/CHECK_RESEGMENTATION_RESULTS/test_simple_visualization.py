#!/usr/bin/env python3
"""
Test script for the simplified segmentation visualization.
"""

import sys
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

try:
    from visualize_segmentation_simple import SimpleSegmentationVisualizer
    
    print("Testing simplified segmentation visualizer...")
    
    # Initialize visualizer
    visualizer = SimpleSegmentationVisualizer()
    
    # Print summary report
    visualizer.print_summary_report()
    
    # Get tiles with cells
    tiles_with_cells = visualizer.get_tiles_with_cells()
    
    print(f"\nFound {len(tiles_with_cells)} tiles with cells:")
    for tile_idx, n_cells in tiles_with_cells:
        print(f"  Tile {tile_idx}: {n_cells} cells")
    
    print("\nSimplified visualization script is ready to use!")
    print("\nRecommended usage:")
    if tiles_with_cells:
        best_tile = tiles_with_cells[0][0]  # Tile with most cells
        print(f"python visualize_segmentation_simple.py --tile-index {best_tile}")
    print("python visualize_segmentation_simple.py --all-tiles")
    print("python visualize_segmentation_simple.py --save-dir ./output")
    
except ImportError as e:
    print(f"Import error: {e}")
    print("Please ensure all required packages are installed:")
    print("pip install numpy pandas matplotlib scikit-image geopandas")
    
except Exception as e:
    print(f"Error: {e}")
    print("Please check that the segmentation results exist in the expected directory.")