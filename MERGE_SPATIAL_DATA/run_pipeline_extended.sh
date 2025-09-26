#!/bin/bash
#SBATCH --job-name=extended
#SBATCH --output=logs/pipeline_extended.out
#SBATCH --error=logs/pipeline_extended.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# --- Pipeline Header ---
echo "========================================================================"
echo "SPATIAL TRANSCRIPTOMICS ANALYSIS - EXTENDED PIPELINE"
echo "========================================================================"
echo "Description: This pipeline analyzes ALL clusters with enhanced statistical power."
echo "             It uses ROI-selected cells and extends small clusters by adding"
echo "             cells from the same cluster across the entire dataset."
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

# --- Step 1: Initial Merge (Preserves all ROI cells) ---
INITIAL_MERGE_FILE="merged_spatial_data_hippo_all_cells.h5ad"
echo ""
echo "--- Step 1: Initial Merge (Preserves all ROI cells) ---"
if [ ! -f "${INITIAL_MERGE_FILE}" ]; then
    echo "Running merge_hippo_extended.py to create ${INITIAL_MERGE_FILE}"
    python ./SCRIPTS/merge_hippo_extended.py
    exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "ERROR: Initial merge failed with exit code ${exit_code}"
        exit ${exit_code}
    fi
    echo "Initial merge completed successfully."
else
    echo "Initial merged data file '${INITIAL_MERGE_FILE}' already exists. Skipping this merge step."
fi

# --- Step 2: Create Extended ROI Dataset ---
EXTENDED_DATA_FILE="merged_spatial_data_hippo_extended.h5ad"
echo ""
echo "--- Step 2: Create Extended ROI Dataset ---"
if [ ! -f "${EXTENDED_DATA_FILE}" ]; then
    echo "Running create_extended_roi_dataset.py to create ${EXTENDED_DATA_FILE}"
    echo "This will extend small clusters by adding cells from the original dataset."
    python ./SCRIPTS/create_extended_roi_dataset.py
    exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "ERROR: Dataset extension failed with exit code ${exit_code}"
        exit ${exit_code}
    fi
    echo "Dataset extension completed successfully."
else
    echo "Extended dataset file '${EXTENDED_DATA_FILE}' already exists. Skipping extension step."
fi

# --- Step 3: Spatial Analysis (Extended - all clusters) ---
echo ""
echo "--- Step 3: Spatial Analysis (Extended - all clusters) ---"
echo "Running hippo_analysis_extended.py"
echo "Output: plots_hippo_cluster_highlight_extended/"
python ./SCRIPTS/hippo_analysis_extended.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "ERROR: Spatial analysis failed with exit code ${exit_code}"
    exit ${exit_code}
fi
echo "Spatial analysis completed successfully."

# --- Step 4: Marker Gene Analysis (Extended - enhanced statistical power) ---
echo ""
echo "--- Step 4: Marker Gene Analysis (Extended - enhanced statistical power) ---"
echo "Running marker_genes_per_sample_extended.py"
echo "Output: marker_genes_per_sample_extended/"
python ./SCRIPTS/marker_genes_per_sample_extended.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "ERROR: Marker gene analysis failed with exit code ${exit_code}"
    exit ${exit_code}
fi
echo "Marker gene analysis completed successfully."

# --- Pipeline Completion Summary ---
echo ""
echo "========================================================================"
echo "EXTENDED PIPELINE COMPLETED SUCCESSFULLY"
echo "========================================================================"
echo "Results:"
echo "  - Initial merged data: ${INITIAL_MERGE_FILE}"
echo "  - Extended dataset: ${EXTENDED_DATA_FILE}"
echo "  - Spatial plots: plots_hippo_cluster_highlight_extended/"
echo "  - Marker genes: marker_genes_per_sample_extended/"
echo ""
echo "Features of this analysis:"
echo "  ✓ Analyzes ALL clusters (no size filtering)"
echo "  ✓ Enhanced statistical power for small clusters"
echo "  ✓ Comprehensive marker gene profiles"
echo "  ⚠ Uses additional cells beyond ROI boundaries for small clusters"
echo "========================================================================"