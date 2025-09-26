# Enhanced Cell Segmentation Visualization Tool

This document provides an overview of the visualization tools for assessing cell segmentation results from the VPT pipeline.

## Overview and Recommendation

The segmentation pipeline has processed the ROI and generated cell boundaries. To visualize these results accurately, we strongly recommend using the **`visualize_segmentation_fixed.py`** script.

### Why Use the Fixed Script?
- **It Works**: Solves critical issues where cell boundaries were invisible or misaligned.
- **Clearer Output**: Creates a comprehensive 3-panel plot (DAPI, PolyT, RGB Overlay) for robust analysis.
- **High Visibility**: Renders cell boundaries with a thick, vibrant yellow line.

Older scripts (`visualize_segmentation_simple.py`, `visualize_segmentation_results.py`) are kept for legacy purposes but are not recommended for use.

## Quick Start

For a fast and reliable visualization of tile 4, run the following command:
```bash
python CHECK_RESEGMENTATION_RESULTS/visualize_segmentation_fixed.py --tile-index 4
```
The output will be saved as `segmentation_visualizations/tile_4_segmentation.png`.

## Requirements

Ensure you have the necessary Python packages installed:
```bash
pip install numpy pandas matplotlib scikit-image geopandas shapely
```

## Recommended Usage

## Recommended Usage with `visualize_segmentation_fixed.py`

This script is the new standard for all visualization tasks.

### 1. Generate a Summary Report
Get a quick overview of cell counts per tile without generating images.
```bash
python CHECK_RESEGMENTATION_RESULTS/visualize_segmentation_fixed.py --summary-only
```

### 2. Visualize a Specific Tile
The most common and recommended command.
```bash
# Visualize tile 4 and save it to the default directory
python CHECK_RESEGMENTATION_RESULTS/visualize_segmentation_fixed.py --tile-index 4 --save-dir segmentation_visualizations
```

### 3. Visualize All Tiles
Process all tiles that contain detected cells.
```bash
python CHECK_RESEGMENTATION_RESULTS/visualize_segmentation_fixed.py --all-tiles --save-dir segmentation_visualizations
```

### 4. Custom Save Location
Save the output to a specific directory.
```bash
python CHECK_RESEGMENTATION_RESULTS/visualize_segmentation_fixed.py --tile-index 4 --save-dir ./my_custom_results
```

## Command Line Options (`visualize_segmentation_fixed.py`)

| Option | Description | Default |
|---|---|---|
| `--output-dir`| Directory with segmentation results | `2task_cellpose2_p30-E165_R1_roi_analysis_parallel` |
| `--tile-index`| Specific tile to visualize | None (interactive mode) |
| `--save-dir`| Directory to save images | `segmentation_visualizations`|
| `--summary-only`| Only show summary report | False |
| `--all-tiles`| Visualize all tiles with cells | False |
| `--sample-n`| Randomly sample N cells for clarity | None |

## Key Visualization Features

The enhanced script provides:
- **3-Panel Plots**: DAPI, PolyT, and a merged RGB overlay for comprehensive analysis.
- **Accurate Overlays**: Cell boundaries are correctly transformed and displayed.
- **High Visibility**: Boundaries are rendered with a thick yellow line.
- **High Resolution**: Images are saved at 300 DPI for detailed inspection.

## Visualization Features

### ROI Overview
- Shows the entire region of interest
- Overlays all detected cell boundaries
- Displays multiple channels (DAPI, PolyT) and RGB overlay
- Provides summary statistics

### Tile-Specific Views
- Focuses on individual tiles for detailed inspection
- Shows original images with segmentation overlays
- Useful for quality assessment of segmentation results
- Highlights areas with high/low cell density

### Interactive Mode
When running without `--tile-index`, the script will:
1. Show ROI overview first
2. Ask if you want to see individual tiles
3. Allow selection of specific tiles or all tiles

## Segmentation Results Summary

Based on your pipeline run:

```
ROI Coordinates:
  X: 10691.9 - 11990.2 microns
  Y: 6190.1 - 6940.4 microns

Tile Information:
  Total tiles: 12
  Tiles with cells: 2
  
Cell Detection:
  Total cells: 56
  Average per tile: 4.7 ± 13.5
  Range: 0-49 cells per tile
```

## Quality Assessment Tips

1. **Check High-Density Areas**: Focus on tile 4 which has 49 cells
   ```bash
   python visualize_segmentation_results.py --tile-index 4
   ```

2. **Verify Boundary Accuracy**: Look for:
   - Proper cell boundary detection
   - No over-segmentation (cells split incorrectly)
   - No under-segmentation (multiple cells merged)

3. **Assess Empty Tiles**: Check why 10 tiles have 0 cells:
   ```bash
   python visualize_segmentation_results.py --tile-index 0
   ```

4. **Compare Channels**: Ensure segmentation uses both DAPI and PolyT appropriately

## Troubleshooting

### Common Issues

1. **Memory Issues / Script Killed**:
   - Use the simplified version: `python visualize_segmentation_simple.py`
   - Focus on individual tiles instead of full ROI
   - The simplified version only loads tiles with detected cells

2. **"'Series' object has no attribute 'geometry'" Error**:
   - This is fixed in the updated versions
   - Use the simplified version for better error handling

3. **Import Errors**:
   - Ensure all required packages are installed in the `vpt` environment
   - Run: `pip install numpy pandas matplotlib scikit-image geopandas shapely`

4. **File Not Found**:
   - Verify that segmentation results exist in the expected directory
   - Check that the output directory name matches your pipeline run

### Getting Help

Run the test scripts to verify everything is working:
```bash
# Test simplified version (recommended)
python test_simple_visualization.py

# Test full version
python test_visualization.py
```

### Performance Tips

1. **Start with the simplified version** - it's more memory efficient
2. **Focus on tiles with cells** - tiles 4 and 8 in your case
3. **Use `--summary-only` first** to check results without loading images
4. **Save visualizations** instead of displaying them to reduce memory usage

## Output Files

When using `--save-dir`, the following files are created:
- `roi_overview_segmentation.png` - Full ROI visualization
- `tile_X_segmentation.png` - Individual tile visualizations

## Technical Details

- **Coordinate System**: Uses mosaic pixel coordinates
- **Image Enhancement**: Applies CLAHE normalization and contrast stretching
- **Segmentation Format**: Reads parquet files with polygon geometries
- **Visualization**: Uses matplotlib with customizable color schemes

## Next Steps

After visual inspection, you may want to:
1. Adjust segmentation parameters if results are unsatisfactory
2. Re-run the pipeline with modified settings
3. Proceed with downstream analysis using the detected cells
4. Export results for further analysis in other tools

---

*Generated for the SRF Spatial Segmentation Pipeline*