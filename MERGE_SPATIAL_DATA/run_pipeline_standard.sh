#!/bin/bash
#SBATCH --job-name=standard
#SBATCH --output=logs/pipeline_standard.out
#SBATCH --error=logs/pipeline_standard.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# --- Pipeline Header ---
echo "========================================================================"
echo "SPATIAL TRANSCRIPTOMICS ANALYSIS - STANDARD PIPELINE"
echo "========================================================================"
echo "Description: This pipeline analyzes hippocampus spatial transcriptomics data."
echo "             It filters out clusters with fewer than 50 cells per sample"
echo "             and uses only ROI-selected cells (no extensions)."
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

# --- Step 1: Merge Hippocampus Data ---
MERGED_DATA_FILE="merged_spatial_data_hippo_standard.h5ad"
echo ""
echo "--- Step 1: Merging Hippocampus Data ---"
if [ ! -f "${MERGED_DATA_FILE}" ]; then
    echo "Running merge_hippo_standard.py to create ${MERGED_DATA_FILE}"
    python ./SCRIPTS/merge_hippo_standard.py
    exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "ERROR: Merging failed with exit code ${exit_code}"
        exit ${exit_code}
    fi
    echo "Merging completed successfully."
else
    echo "Merged data file '${MERGED_DATA_FILE}' already exists. Skipping merge step."
fi

# --- Step 2: Spatial Analysis ---
echo ""
echo "--- Step 2: Spatial Analysis (Filters clusters <50 cells) ---"
echo "Running hippo_analysis_standard.py"
echo "Output: plots_hippo_cluster_highlight_standard/"
python ./SCRIPTS/hippo_analysis_standard.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "ERROR: Spatial analysis failed with exit code ${exit_code}"
    exit ${exit_code}
fi
echo "Spatial analysis completed."

# --- Step 3: Marker Gene Analysis ---
echo ""
echo "--- Step 3: Marker Gene Analysis (Filters clusters <50 cells) ---"
echo "Running marker_genes_per_sample_standard.py"
echo "Output: marker_genes_per_sample_standard/"
python ./SCRIPTS/marker_genes_per_sample_standard.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "ERROR: Marker gene analysis failed with exit code ${exit_code}"
    exit ${exit_code}
fi
echo "Marker gene analysis completed."

# --- Pipeline Completion Summary ---
echo ""
echo "========================================================================"
echo "STANDARD PIPELINE COMPLETED SUCCESSFULLY"
echo "========================================================================"
echo "Results:"
echo "  - Merged data: ${MERGED_DATA_FILE}"
echo "  - Spatial plots: plots_hippo_cluster_highlight_standard/"
echo "  - Marker genes: marker_genes_per_sample_standard/"
echo ""
echo "Features of this analysis:"
echo "  ✓ Uses only ROI-selected cells"
echo "  ✓ Filters out clusters with <50 cells per sample"
echo "  ✓ High statistical rigor for analyzed clusters"
echo "  ✗ May miss small but biologically relevant clusters"
echo "========================================================================"