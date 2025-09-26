#!/bin/bash
#SBATCH --job-name=debug_spatial_coords
#SBATCH --output=logs/debug_spatial_coords.out
#SBATCH --error=logs/debug_spatial_coords.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Debug script to check spatial coordinates structure in data files
echo "Starting debug script to check spatial coordinates..."

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment (using same environment as existing pipeline)
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

# Set the working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation

echo "Running debug script..."
python MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/debug_spatial_coords.py

echo "Debug script completed."