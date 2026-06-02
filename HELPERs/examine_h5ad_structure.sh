#!/bin/bash
#SBATCH --job-name=examine_h5ad_structure
#SBATCH --output=logs/examine_h5ad_structure.out
#SBATCH --error=logs/examine_h5ad_structure.err
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
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation

# Create logs directory if it doesn't exist
mkdir -p logs

# Run the focused analysis pipeline
python examine_h5ad_structure.py

echo "Analysis completed successfully!"