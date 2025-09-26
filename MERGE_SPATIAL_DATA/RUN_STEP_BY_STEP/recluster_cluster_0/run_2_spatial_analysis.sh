#!/bin/bash
#SBATCH --job-name=run_2_spatial_analysis
#SBATCH --output=logs/run_2_spatial_analysis.out
#SBATCH --error=logs/run_2_spatial_analysis.err
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
echo "SPATIAL TRANSCRIPTOMICS ANALYSIS - STEP 2: SPATIAL ANALYSIS (INCLUSIVE)"
echo "========================================================================"

# --- Step 2: Spatial Analysis (Inclusive - all clusters) ---
echo ""
echo "--- Step 2: Spatial Analysis (Inclusive - all clusters) ---"
echo "Running hippo_analysis_inclusive.py"
echo "Output: plots_hippo_cluster_highlight_inclusive/"
python /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation/MERGE_SPATIAL_DATA/SCRIPTS/inclusive/recluster/2_hippo_analysis.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "ERROR: Spatial analysis failed with exit code ${exit_code}"
    exit ${exit_code}
fi
echo "Spatial analysis completed successfully."

echo "========================================================================"
echo "STEP 2 COMPLETED"
echo "========================================================================"