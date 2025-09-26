# ROI-Based Temporal Markers Analysis

## Overview

This directory contains an alternative approach for temporal markers analysis that uses **ROI-based spatial binning** rather than individual segmented cells. This method is independent of cell segmentation quality and includes all transcripts detected within ROI boundaries.

## Key Differences from Cell-Based Analysis

| Aspect | Cell-Based Analysis | ROI-Based Analysis |
|--------|-------------------|-------------------|
| **Input Data** | Segmented cells from merged data | Original spatial transcripts |
| **Statistical Units** | Individual cells | Spatial bins (100μm squares) |
| **Transcript Coverage** | Only transcripts assigned to cells | All transcripts within ROI |
| **Segmentation Dependency** | ✗ Depends on segmentation quality | ✓ Independent of segmentation |
| **Spatial Context** | Cell-level | Bin-level within ROI |
| **Statistical Test** | Wilcoxon rank-sum (cell-level) | Mann-Whitney U (bin-level) |

## Analysis Workflow

### 1. ROI Geometry Loading
- Uses existing ROI definitions from `DATA/ROIs/{sample_id}_geometry.csv`
- Filters to "hippo" ROI group (matching existing pipeline)
- Leverages same geometric framework as `Subset_ROIs.py`

### 2. Spatial Binning
- Creates regular 100μm × 100μm spatial bins within ROI boundaries
- Only includes bins that intersect with ROI polygons
- Configurable bin size (50μm, 100μm, or 200μm options)

### 3. Transcript Aggregation
- Sums transcript counts from all points within each spatial bin
- Applies minimum transcript threshold (default: 10 transcripts per bin)
- Maintains sample metadata (time_point, type, region)

### 4. Statistical Analysis
- Compares ctrl vs mut samples at each time point using bin-level data
- Mann-Whitney U test for non-parametric comparisons
- FDR correction for multiple testing
- Same significance thresholds as cell-based analysis

### 5. Visualization
- Volcano plots showing differential expression
- Spatial maps of bin locations and transcript density
- Summary statistics across time points

## Files and Usage

### Main Scripts
- `roi_based_temporal_markers_analysis.py` - Main analysis script
- `run_roi_based_temporal_markers_analysis.sh` - SLURM submission script

### Configuration Parameters
```python
# Spatial binning parameters
BIN_SIZE = 100  # μm (adjustable: 50, 100, 200)
MIN_TRANSCRIPTS_PER_BIN = 10  # Minimum for statistical analysis
ROI_GROUP = "hippo"  # ROI group to analyze

# Statistical parameters
PVAL_THRESHOLD = 0.05
LOGFC_THRESHOLD = 0.5
MULTIPLE_TESTING_METHOD = "fdr_bh"
```

### Execution
```bash
# Submit to cluster
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation
sbatch MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/run_roi_based_temporal_markers_analysis.sh

# Monitor job
squeue -u $USER
```

## Output Files

### Expression Data
- `roi_bin_expression_{sample_id}.csv` - Expression matrix for spatial bins
- `roi_bin_metadata_{sample_id}.csv` - Bin coordinates and statistics

### Differential Expression Results
- `roi_de_results_{time_point}_mut_vs_ctrl.csv` - DE results per time point
- `roi_combined_results_all_timepoints.csv` - Combined results across all time points

### Visualizations
- `roi_volcano_plot_{time_point}_mut_vs_ctrl.png` - Volcano plots
- `roi_spatial_bins_{time_point}.png` - Spatial bin maps

### Summary
- `roi_summary_statistics.csv` - Summary statistics per time point
- `roi_analysis_summary.txt` - Complete analysis report

## Input Data Requirements

### Spatial Transcriptomics Data
- `DATA/{time_group}/{region}/data.h5ad` - Original Vizgen spatial data
- Must contain spatial coordinates in `obs['center_x']` and `obs['center_y']` columns
- Gene expression matrix in `adata.X`
- Data structure includes: `volume`, `center_x`, `center_y`, `n_counts`, `n_genes`, `leiden` columns

### ROI Definitions
- `DATA/ROIs/{sample_id}_geometry.csv` - ROI geometry files
- Must contain 'geometry' column with WKT format polygons
- Must contain 'group' column with ROI classifications

### Sample Structure
Uses same `data_map` as existing pipeline:
```python
data_map = {
    'p0-p7': {
        'R1': {'time_point': 'p7', 'region': 'R1', 'type': 'mut'},
        'R2': {'time_point': 'p7', 'region': 'R2', 'type': 'ctrl'},
        # ... additional samples
    },
    'p30-E165': {
        # ... additional time groups
    }
}
```

## Expected Results

### Time Point Comparisons
- **p0**: mut (R5) vs ctrl (R6)
- **p7**: mut (R1) vs ctrl (R2)
- **p30**: mut (R3) vs ctrl (R4)
- **e16**: mut (R1) vs ctrl (R2)

### Statistical Output
For each time point:
- Number of spatial bins analyzed
- Total genes tested
- Significantly upregulated genes
- Significantly downregulated genes
- Percentage of significant genes

### Spatial Context
- Bin-level expression maps showing spatial patterns
- ROI boundary overlays
- Transcript density visualization

## Advantages of ROI-Based Approach

1. **Segmentation Independence**: Results don't depend on cell segmentation algorithms
2. **Complete Coverage**: Includes all transcripts in tissue areas, not just cell-assigned ones
3. **Spatial Preservation**: Maintains spatial relationships within ROI boundaries
4. **Robust Statistics**: Aggregation reduces noise from individual transcript detection
5. **Consistency**: Uses same ROI definitions as existing pipeline
6. **Complementary**: Can be compared with cell-based results for validation

## Comparison with Cell-Based Analysis

To compare results:
1. Run both analyses: cell-based (`temporal_markers_analysis.py`) and ROI-based
2. Compare significant genes lists between approaches
3. Examine spatial patterns in ROI-based visualizations
4. Validate findings using both methods

## Troubleshooting

### Common Issues
1. **Missing ROI files**: Ensure ROI selection was completed using `Subset_ROIs.py`
2. **Spatial coordinates not found**: Error "AnnData object must contain spatial coordinates"
   - **Cause**: Script looking for coordinates in wrong location
   - **Solution**: Coordinates are stored in `obs['center_x']` and `obs['center_y']`, not `obsm['spatial']`
   - **Fixed in**: Current version of `roi_based_temporal_markers_analysis.py`
3. **Memory errors**: Reduce bin size or increase allocated memory
4. **No significant results**: Check if samples have sufficient spatial coverage
5. **Coordinate system mismatches**: Verify spatial coordinates match ROI geometries
6. **Path errors**: Ensure all paths use correct project root `/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation`

### Resource Requirements
- **Memory**: 128GB recommended for full dataset
- **Time**: ~2-6 hours depending on data size
- **Storage**: ~1-2GB for output files

## Integration with Existing Pipeline

This ROI-based analysis integrates seamlessly with the existing spatial analysis pipeline:

1. **Prerequisite**: ROI selection must be completed (`Subset_ROIs.py`)
2. **Parallel Analysis**: Can run alongside cell-based temporal markers analysis
3. **Complementary Results**: Provides alternative perspective on spatial gene expression
4. **Same Infrastructure**: Uses same compute environment and data structure

## Citation and Methods

### Methods Description for Publications
"ROI-based differential expression analysis was performed using spatial binning within predefined hippocampus regions. Original spatial transcriptomics data was divided into regular 100μm × 100μm spatial bins, and transcript counts were aggregated within each bin. Statistical comparisons between control and mutant samples were performed using Mann-Whitney U tests on bin-level expression data, with false discovery rate correction for multiple testing. This approach provides analysis independent of cell segmentation quality while preserving spatial context within regions of interest."

### Technical Parameters
- **Spatial bin size**: 100μm × 100μm
- **Minimum transcripts per bin**: 10
- **Statistical test**: Mann-Whitney U test
- **Multiple testing correction**: FDR (Benjamini-Hochberg)
- **Significance thresholds**: p-adj < 0.05, |log2FC| > 0.5

## Recent Updates

### Version 1.1 - Spatial Coordinates Fix
**Issue Resolved**: Initial version incorrectly looked for spatial coordinates in `adata.obsm['spatial']`

**Root Cause**: Vizgen spatial data stores coordinates in observation columns, not as multidimensional array

**Fix Applied**:
- Updated `aggregate_transcripts_to_bins()` function
- Now correctly reads coordinates from `adata.obs['center_x']` and `adata.obs['center_y']`
- Added proper error checking for coordinate column existence

**Impact**: Analysis now runs successfully and can process all spatial samples

**Files Updated**:
- `roi_based_temporal_markers_analysis.py` - Fixed coordinate access
- `README_ROI_BASED_ANALYSIS.md` - Updated documentation
- Troubleshooting section added for future reference


Summary of ROI-Based Analysis Implementation

Key Answer to Your Original Question

The original temporal markers analysis (temporal_markers_analysis.py) derives DEGs only from genes within segmented cells, making results dependent on cell segmentation quality.

The new ROI-based analysis (roi_based_temporal_markers_analysis.py) derives DEGs from all transcripts within the spatial ROI regions, independent of cell segmentation.

Files Created

1. roi_based_temporal_markers_analysis.py - Main analysis script that:
- Uses your existing ROI framework from Subset_ROIs.py
- Creates spatial bins within ROI boundaries
- Aggregates transcripts at the bin level
- Performs statistical analysis independent of segmentation
2. run_roi_based_temporal_markers_analysis.sh - SLURM submission script with:
- Proper resource allocation (128GB RAM, 6 hours)
- Integration with existing conda environment
- Comprehensive error checking and logging
3. README_ROI_BASED_ANALYSIS.md - Complete documentation including:
- Comparison with cell-based approach
- Usage instructions and parameters
- Expected outputs and interpretation
- Integration with existing pipeline

Key Features of the ROI-Based Approach

✅ Independent of Cell Segmentation - Uses original spatial transcript data
✅ Complete Transcript Coverage - Includes all transcripts within ROI boundaries✅ Spatial Context Preserved - Uses 100μm spatial bins within ROI regions
✅ Integrates with Existing Pipeline - Uses same ROI definitions and data structure
✅ Robust Statistics - Mann-Whitney U test on aggregated bin-level data
✅ Comprehensive Visualization - Volcano plots and spatial heatmaps

Usage

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation
sbatch MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/run_roi_based_temporal_markers_analysis.sh

The new approach provides an alternative analysis method that gives you insight into differential gene expression patterns within your ROI regions without being limited by cell segmentation boundaries, while still leveraging your
existing ROI selection infrastructure.