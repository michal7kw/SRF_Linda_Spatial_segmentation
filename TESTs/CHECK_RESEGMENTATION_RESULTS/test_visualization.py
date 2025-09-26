#!/usr/bin/env python3
"""
Simple test script to verify the segmentation visualization works.
"""

import sys
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

try:
    from visualize_segmentation_results import SegmentationVisualizer
    
    print("Testing segmentation visualizer...")
    
    # Initialize visualizer
    visualizer = SegmentationVisualizer()
    
    # Print summary report
    visualizer.print_summary_report()
    
    print("\nVisualization script is ready to use!")
    print("\nTo run the full visualization:")
    print("python visualize_segmentation_results.py")
    print("\nTo visualize a specific tile (e.g., tile 0):")
    print("python visualize_segmentation_results.py --tile-index 0")
    print("\nTo save visualizations:")
    print("python visualize_segmentation_results.py --save-dir ./visualization_output")
    
except ImportError as e:
    print(f"Import error: {e}")
    print("Please ensure all required packages are installed:")
    print("pip install numpy pandas matplotlib seaborn scikit-image geopandas shapely")
    
except Exception as e:
    print(f"Error: {e}")
    print("Please check that the segmentation results exist in the expected directory.")