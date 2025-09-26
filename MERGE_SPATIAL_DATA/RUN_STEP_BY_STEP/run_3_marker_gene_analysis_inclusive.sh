#!/bin/bash
#SBATCH --job-name=run_3_marker_gene_analysis_inclusive
#SBATCH --output=logs/run_3_marker_gene_analysis_inclusive.out
#SBATCH --error=logs/run_3_marker_gene_analysis_inclusive.err
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
echo "SPATIAL TRANSCRIPTOMICS ANALYSIS - STEP 3: MARKER GENE ANALYSIS (INCLUSIVE)"
echo "========================================================================"

# --- Step 3: Marker Gene Analysis (Inclusive - all clusters) ---
echo ""
echo "--- Step 3: Marker Gene Analysis (Inclusive - all clusters) ---"
echo "Running marker_genes_per_sample_inclusive.py"
echo "Output: marker_genes_per_sample_inclusive/"
python /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation/MERGE_SPATIAL_DATA/SCRIPTS/inclusive/3_marker_genes_per_sample.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "ERROR: Marker gene analysis failed with exit code ${exit_code}"
    exit ${exit_code}
fi
echo "Marker gene analysis completed successfully."

echo "========================================================================"
echo "STEP 3 COMPLETED"
echo "========================================================================"