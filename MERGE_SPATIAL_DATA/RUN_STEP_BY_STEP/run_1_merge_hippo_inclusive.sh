#!/bin/bash
#SBATCH --job-name=run_1_merge_hippo_inclusive
#SBATCH --output=logs/run_1_merge_hippo_inclusive.out
#SBATCH --error=logs/run_1_merge_hippo_inclusive.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:30:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# --- Pipeline Header ---
echo "========================================================================"
echo "SPATIAL TRANSCRIPTOMICS ANALYSIS - COMMON SETUP"
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

echo "========================================================================"
echo "SPATIAL TRANSCRIPTOMICS ANALYSIS - STEP 1: MERGE HIPPOCAMPUS DATA (INCLUSIVE)"
echo "========================================================================"

# --- Step 1: Merge Hippocampus Data (Inclusive) ---
MERGED_DATA_FILE="merged_spatial_data_hippo_inclusive.h5ad"
echo ""
echo "--- Step 1: Merging Hippocampus Data (Inclusive) ---"
if [ ! -f "${MERGED_DATA_FILE}" ]; then
    echo "Running merge_hippo_inclusive.py to create ${MERGED_DATA_FILE}"
    python /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation/MERGE_SPATIAL_DATA/SCRIPTS/inclusive/1_merge_hippo.py
    exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "ERROR: Merging failed with exit code ${exit_code}"
        exit ${exit_code}
    fi
    echo "Merging completed successfully."
else
    echo "Merged data file '${MERGED_DATA_FILE}' already exists. Skipping merge step."
fi

echo "========================================================================"
echo "STEP 1 COMPLETED"
echo "========================================================================"