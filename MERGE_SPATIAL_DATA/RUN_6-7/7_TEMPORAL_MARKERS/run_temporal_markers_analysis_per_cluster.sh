#!/bin/bash
#SBATCH --job-name=run_temporal_markers_analysis_per_cluster
#SBATCH --output=logs/run_temporal_markers_analysis_per_cluster.out
#SBATCH --error=logs/run_temporal_markers_analysis_per_cluster.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Cluster-specific Temporal Markers Analysis Pipeline
# This script runs the cluster-specific temporal markers analysis for spatial transcriptomics data
# Compares ctrl vs mut samples within leiden_new clusters at corresponding time points
# Only analyzes clusters with at least 20 cells in both conditions

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "Starting Cluster-specific Temporal Markers Analysis Pipeline..."
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
mkdir -p MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/cluster_results
mkdir -p logs

# Display analysis information
echo "================================================================"
echo "CLUSTER-SPECIFIC TEMPORAL MARKERS ANALYSIS"
echo "================================================================"
echo "This analysis will compare ctrl vs mut samples within leiden_new clusters:"
echo "- Only clusters with ≥20 cells in both ctrl and mut will be analyzed"
echo "- Separate analysis for each valid cluster at each time point"
echo "- Time points: p0, p7, p30, e16"
echo "================================================================"

# Run the analysis
echo "Running cluster-specific temporal markers analysis..."
python MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/cluster_specific_markers_analysis.py

# Check if analysis completed successfully
if [ $? -eq 0 ]; then
    echo "================================================================"
    echo "Cluster-specific Pipeline completed successfully!"
    echo "================================================================"
    echo "Results saved to:"
    echo "- Individual cluster results: MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/cluster_results/"
    echo "- Volcano plots: MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/cluster_results/volcano_plot_*_cluster_*.png"
    echo "- Heatmaps: MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/cluster_results/heatmap_*_cluster_*.png"
    echo "- Combined results: MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/cluster_results/combined_cluster_de_results_all_timepoints.csv"
    echo "- Cluster summary: MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/cluster_results/cluster_summary_statistics.csv"
    echo "- Summary plots: MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/cluster_results/cluster_analysis_summary.png"
    echo "================================================================"
else
    echo "================================================================"
    echo "Pipeline failed! Check the error logs for details."
    echo "================================================================"
    exit 1
fi