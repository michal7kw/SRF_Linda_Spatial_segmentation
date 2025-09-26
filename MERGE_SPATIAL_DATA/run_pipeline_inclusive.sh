#!/bin/bash
#SBATCH --job-name=inclusive
#SBATCH --output=logs/pipeline_inclusive.out
#SBATCH --error=logs/pipeline_inclusive.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:30:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# --- Pipeline Header ---
echo "========================================================================"
echo "SPATIAL TRANSCRIPTOMICS ANALYSIS - INCLUSIVE PIPELINE"
echo "========================================================================"
echo "Description: This pipeline analyzes ALL clusters identified within the 'hippo' ROI."
echo "             It does not filter clusters based on cell count."
echo "========================================================================"

# --- Setup ---
# Create logs directory if it doesn't exist
mkdir -p logs

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

# Change to the working directory
cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation/MERGE_SPATIAL_DATA"

# --- Step 1: Merge Hippocampus Data (Inclusive) ---
MERGED_DATA_FILE="merged_spatial_data_hippo_inclusive.h5ad"
echo ""
echo "--- Step 1: Merging Hippocampus Data (Inclusive) ---"
if [ ! -f "${MERGED_DATA_FILE}" ]; then
    echo "Running merge_hippo_inclusive.py to create ${MERGED_DATA_FILE}"
    python ./SCRIPTS/merge_hippo_inclusive.py
    exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "ERROR: Merging failed with exit code ${exit_code}"
        exit ${exit_code}
    fi
    echo "Merging completed successfully."
else
    echo "Merged data file '${MERGED_DATA_FILE}' already exists. Skipping merge step."
fi

# --- Step 2: Spatial Analysis (Inclusive - all clusters) ---
echo ""
echo "--- Step 2: Spatial Analysis (Inclusive - all clusters) ---"
echo "Running hippo_analysis_inclusive.py"
echo "Output: plots_hippo_cluster_highlight_inclusive/"
python ./SCRIPTS/hippo_analysis_inclusive.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "ERROR: Spatial analysis failed with exit code ${exit_code}"
    exit ${exit_code}
fi
echo "Spatial analysis completed successfully."

# --- Step 3: Marker Gene Analysis (Inclusive - all clusters) ---
echo ""
echo "--- Step 3: Marker Gene Analysis (Inclusive - all clusters) ---"
echo "Running marker_genes_per_sample_inclusive.py"
echo "Output: marker_genes_per_sample_inclusive/"
python ./SCRIPTS/marker_genes_per_sample_inclusive.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "ERROR: Marker gene analysis failed with exit code ${exit_code}"
    exit ${exit_code}
fi
echo "Marker gene analysis completed successfully."

# --- Pipeline Completion Summary ---
echo ""
echo "========================================================================"
echo "INCLUSIVE PIPELINE COMPLETED SUCCESSFULLY"
echo "========================================================================"
echo "Results:"
echo "  - Merged data: ${MERGED_DATA_FILE}"
echo "  - Spatial plots: plots_hippo_cluster_highlight_inclusive/"
echo "  - Marker genes: marker_genes_per_sample_inclusive/"
echo ""
echo "Features of this analysis:"
echo "  ✓ Analyzes ALL clusters from the 'hippo' ROI"
echo "  ✓ No cluster size filtering"
echo "  ✓ Provides a complete picture of all cell populations in the ROI"
echo "  ⚠ Small clusters (<50 cells) are analyzed, but may have limited statistical power"
echo "========================================================================"