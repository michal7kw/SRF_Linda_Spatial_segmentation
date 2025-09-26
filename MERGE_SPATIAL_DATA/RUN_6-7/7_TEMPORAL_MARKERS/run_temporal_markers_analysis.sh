#!/bin/bash
#SBATCH --job-name=run_temporal_markers_analysis
#SBATCH --output=logs/run_temporal_markers_analysis.out
#SBATCH --error=logs/run_temporal_markers_analysis.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Temporal Markers Analysis Pipeline
# This script runs the temporal markers analysis for spatial transcriptomics data
# Compares ctrl vs mut samples at corresponding time points

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "Starting Temporal Markers Analysis Pipeline..."
echo "================================================================"

# Set the working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation

# Check if required files exist
echo "Checking required files..."

if [ ! -f "MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad" ]; then
    echo "Error: Merged spatial data file not found!"
    echo "Please run the merge script first: MERGE_SPATIAL_DATA/SCRIPTS/inclusive/1_merge_hippo.py"
    exit 1
fi

echo "All required files found."

# Create output directories
echo "Creating output directories..."
mkdir -p MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/results
mkdir -p logs

# Display analysis information
echo "================================================================"
echo "TEMPORAL MARKERS ANALYSIS"
echo "================================================================"
echo "This analysis will compare ctrl vs mut samples at each time point:"
echo "- p0: mut (R5) vs ctrl (R6)"
echo "- p7: mut (R1) vs ctrl (R2)"
echo "- p30: mut (R3) vs ctrl (R4)"
echo "- e16: mut (R1) vs ctrl (R2)"
echo "================================================================"

# Run the analysis
echo "Running temporal markers analysis..."
python MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/temporal_markers_analysis.py

# Check if analysis completed successfully
if [ $? -eq 0 ]; then
    echo "================================================================"
    echo "Pipeline completed successfully!"
    echo "================================================================"
    echo "Results saved to:"
    echo "- Individual results: MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/results/"
    echo "- Volcano plots: MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/results/volcano_plot_*.png"
    echo "- Heatmaps: MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/results/heatmap_*.png"
    echo "- Combined results: MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/results/combined_de_results_all_timepoints.csv"
    echo "- Summary statistics: MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/results/summary_statistics.csv"
    echo "================================================================"
else
    echo "================================================================"
    echo "Pipeline failed! Check the error logs for details."
    echo "================================================================"
    exit 1
fi