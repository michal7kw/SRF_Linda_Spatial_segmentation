#!/usr/bin/env python3
"""
Test script for the fixed visualization script.
"""

import sys
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path.cwd()))

from visualize_segmentation_fixed import FixedSegmentationVisualizer

def test_fixed_visualization():
    """Test the fixed visualization script."""
    print("Testing Fixed Segmentation Visualizer")
    print("=" * 50)
    
    try:
        # Initialize visualizer
        visualizer = FixedSegmentationVisualizer()
        
        # Print summary report
        visualizer.print_summary_report()
        
        # Get tiles with cells
        tiles_with_cells = visualizer.get_tiles_with_cells()
        
        if not tiles_with_cells:
            print("No tiles with cells found for testing.")
            return
        
        # Test visualization on the tile with most cells
        test_tile = tiles_with_cells[0][0]  # First tile (highest cell count)
        print(f"\nTesting visualization on tile {test_tile}...")
        
        # Create output directory
        output_dir = Path("test_fixed_visualizations")
        output_dir.mkdir(exist_ok=True)
        
        # Test visualization
        save_path = output_dir / f"test_tile_{test_tile}_fixed.png"
        
        visualizer.visualize_tile_segmentation(
            tile_index=test_tile,
            channels=['DAPI', 'PolyT'],
            z_layer=6,
            save_path=str(save_path)
        )
        
        print(f"\nTest completed successfully!")
        print(f"Check the output image: {save_path}")
        
    except Exception as e:
        print(f"Error during testing: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_fixed_visualization()