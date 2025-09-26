#!/bin/bash
#SBATCH --job-name=run_3_marker_gene_analysis
#SBATCH --output=logs/run_3_marker_gene_analysis.out
#SBATCH --error=logs/run_3_marker_gene_analysis.err
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
echo "SPATIAL TRANSCRIPTOMICS ANALYSIS - STEP 3: MARKER GENE ANALYSIS WITH CLUSTER 0 RECLUSTERING"
echo "========================================================================"

# --- Step 3: Marker Gene Analysis (Inclusive with cluster 0 reclustering) ---
echo ""
echo "--- Step 3: Marker Gene Analysis (Inclusive with cluster 0 reclustering) ---"
echo "Running marker_genes_per_sample_inclusive_cluster0_recluster.py"
echo "Output: marker_genes_per_sample_inclusive_cluster0_recluster/"
python /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation/MERGE_SPATIAL_DATA/SCRIPTS/inclusive/recluster/3_marker_genes_per_sample.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "ERROR: Marker gene analysis with cluster 0 reclustering failed with exit code ${exit_code}"
    exit ${exit_code}
fi
echo "Marker gene analysis with cluster 0 reclustering completed successfully."

echo "========================================================================"
echo "STEP 3 WITH CLUSTER 0 RECLUSTERING COMPLETED"
echo "========================================================================"