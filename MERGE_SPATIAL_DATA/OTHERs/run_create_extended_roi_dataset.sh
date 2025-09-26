#!/bin/bash
#SBATCH --job-name=create_extended_roi
#SBATCH --output=logs/create_extended_roi.out
#SBATCH --error=logs/create_extended_roi.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Create logs directory
mkdir -p logs

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "Creating extended ROI dataset..."
cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation/MERGE_SPATIAL_DATA"

python ./create_extended_roi_dataset.py