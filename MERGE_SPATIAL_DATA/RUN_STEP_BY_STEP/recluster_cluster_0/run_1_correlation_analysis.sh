#!/bin/bash
#SBATCH --job-name=run_1_correlation_analysis
#SBATCH --output=logs/run_1_correlation_analysis.out
#SBATCH --error=logs/run_1_correlation_analysis.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=03:00:00
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
echo "SPATIAL TRANSCRIPTOMICS ANALYSIS - COMBINED ANALYSIS WITH CLUSTER 0 RECLUSTERING"
echo "========================================================================"

# --- Combined Analysis: Cluster 0 Reclustering Across All Samples ---
echo ""
echo "--- Combined Analysis: Cluster 0 Reclustering Across All Samples ---"
echo "Running combined_analysis_cluster0_recluster.py"
echo "Output: combined_analysis_cluster0_recluster/"
python /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation/MERGE_SPATIAL_DATA/SCRIPTS/inclusive/recluster/1_correlation_analysis.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "ERROR: Combined analysis with cluster 0 reclustering failed with exit code ${exit_code}"
    exit ${exit_code}
fi
echo "Combined analysis with cluster 0 reclustering completed successfully."

echo "========================================================================"
echo "COMBINED ANALYSIS WITH CLUSTER 0 RECLUSTERING COMPLETED"
echo "========================================================================"