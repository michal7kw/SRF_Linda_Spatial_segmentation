#!/bin/bash
#SBATCH --job-name=roi_temporal_markers
#SBATCH --output=logs/roi_temporal_markers_analysis.out
#SBATCH --error=logs/roi_temporal_markers_analysis.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# ROI-Based Temporal Markers Analysis Pipeline
# This script runs the ROI-based temporal markers analysis for spatial transcriptomics data
# Compares ctrl vs mut samples at corresponding time points using spatial binning approach
# Independent of cell segmentation quality

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment (using same environment as existing pipeline)
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "Starting ROI-Based Temporal Markers Analysis Pipeline..."
echo "================================================================"

# Set the working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation

# Check if required directories and files exist
echo "Checking required files and directories..."

# Check for original spatial data files
required_files=(
    "DATA/p0-p7/R1/data.h5ad"
    "DATA/p0-p7/R2/data.h5ad"
    "DATA/p0-p7/R4/data.h5ad"
    "DATA/p0-p7/R5/data.h5ad"
    "DATA/p0-p7/R6/data.h5ad"
    "DATA/p30-E165/R1/data.h5ad"
    "DATA/p30-E165/R2/data.h5ad"
    "DATA/p30-E165/R3/data.h5ad"
)

missing_files=0
for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Warning: Data file not found: $file"
        ((missing_files++))
    fi
done

# Check for ROI geometry files
roi_files=(
    "DATA/ROIs/p7_R1_mut_geometry.csv"
    "DATA/ROIs/p7_R2_ctrl_geometry.csv"
    "DATA/ROIs/p30_R4p_ctrl_geometry.csv"
    "DATA/ROIs/p0_R5p_mut_geometry.csv"
    "DATA/ROIs/p0_R6_ctrl_geometry.csv"
    "DATA/ROIs/e16_R1_mut_geometry.csv"
    "DATA/ROIs/e16_R2_ctrl_geometry.csv"
    "DATA/ROIs/p30_R3_mut_geometry.csv"
)

missing_rois=0
for file in "${roi_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Warning: ROI file not found: $file"
        ((missing_rois++))
    fi
done

echo "Missing data files: $missing_files"
echo "Missing ROI files: $missing_rois"

if [ $missing_files -gt 0 ] || [ $missing_rois -gt 0 ]; then
    echo "Warning: Some required files are missing. Analysis will continue with available data."
fi

# Create output directories
echo "Creating output directories..."
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/roi_results
mkdir -p logs

# Display analysis information
echo "================================================================"
echo "ROI-BASED TEMPORAL MARKERS ANALYSIS"
echo "================================================================"
echo "This analysis uses spatial binning within ROI boundaries and compares:"
echo "- p0: mut (R5) vs ctrl (R6)"
echo "- p7: mut (R1) vs ctrl (R2)"
echo "- p30: mut (R3) vs ctrl (R4)"
echo "- e16: mut (R1) vs ctrl (R2)"
echo ""
echo "Analysis parameters:"
echo "- Spatial bin size: 100 μm"
echo "- Minimum transcripts per bin: 10"
echo "- ROI group: hippo"
echo "- Statistical test: Mann-Whitney U"
echo "- Multiple testing correction: FDR (Benjamini-Hochberg)"
echo "================================================================"

# Run the ROI-based analysis
echo "Running ROI-based temporal markers analysis..."
python MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/roi_based_temporal_markers_analysis.py

# Check if analysis completed successfully
if [ $? -eq 0 ]; then
    echo "================================================================"
    echo "ROI-BASED ANALYSIS COMPLETED SUCCESSFULLY!"
    echo "================================================================"
    echo "Results saved to:"
    echo "- Individual bin expression matrices: MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/roi_results/roi_bin_expression_*.csv"
    echo "- Bin metadata: MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/roi_results/roi_bin_metadata_*.csv"
    echo "- DE results per time point: MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/roi_results/roi_de_results_*_mut_vs_ctrl.csv"
    echo "- Volcano plots: MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/roi_results/roi_volcano_plot_*.png"
    echo "- Spatial bin plots: MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/roi_results/roi_spatial_bins_*.png"
    echo "- Combined results: MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/roi_results/roi_combined_results_all_timepoints.csv"
    echo "- Summary statistics: MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/roi_results/roi_summary_statistics.csv"
    echo "- Analysis report: MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/roi_results/roi_analysis_summary.txt"
    echo "================================================================"
    echo ""
    echo "Key differences from cell-based analysis:"
    echo "✓ Independent of cell segmentation quality"
    echo "✓ Includes all transcripts within ROI boundaries"
    echo "✓ Uses spatial binning for robust statistical analysis"
    echo "✓ Preserves spatial context within ROI regions"
    echo ""
    echo "Next steps:"
    echo "1. Compare results with cell-based analysis in MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/results/"
    echo "2. Review volcano plots and spatial visualizations"
    echo "3. Analyze combined results across all time points"
    echo "================================================================"
else
    echo "================================================================"
    echo "ROI-BASED ANALYSIS FAILED! Check the error logs for details."
    echo "================================================================"
    echo "Troubleshooting steps:"
    echo "1. Check that all required data files exist"
    echo "2. Verify ROI geometry files are properly formatted"
    echo "3. Ensure conda environment has all required packages"
    echo "4. Check available memory and disk space"
    echo "5. Review error log: logs/roi_temporal_markers_analysis.err"
    echo "================================================================"
    exit 1
fi