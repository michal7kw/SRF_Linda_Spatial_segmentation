#!/bin/bash
#SBATCH --job-name=4_improved_normalization
#SBATCH --output=logs/4_improved_normalization.out
#SBATCH --error=logs/4_improved_normalization.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

# Set the working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation/MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS

# Create logs directory if it doesn't exist
mkdir -p logs

# Run the focused analysis pipeline
echo "Running focused ctrl vs mut analysis pipeline..."
python 4_improved_normalization.py

echo "Analysis completed successfully!"