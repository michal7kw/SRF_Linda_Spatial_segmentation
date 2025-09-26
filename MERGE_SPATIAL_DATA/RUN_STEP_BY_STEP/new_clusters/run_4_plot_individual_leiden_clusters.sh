#!/bin/bash
#SBATCH --job-name=run_4_plot_individual_leiden_clusters
#SBATCH --output=logs/run_4_plot_individual_leiden_clusters.out
#SBATCH --error=logs/run_4_plot_individual_leiden_clusters.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:30:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

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

echo "Running 4_plot_individual_leiden_clusters.py"
python /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation/MERGE_SPATIAL_DATA/SCRIPTS/inclusive/new_clusters/4_plot_individual_leiden_clusters.py

echo "========================================================================"
echo "STEP 4 COMPLETED"
echo "========================================================================"