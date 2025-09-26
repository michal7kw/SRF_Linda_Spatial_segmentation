#!/bin/bash
#SBATCH --job-name=roi_selection_missing_cell_exploration
#SBATCH --output=logs/roi_selection_missing_cell_exploration.out
#SBATCH --error=logs/roi_selection_missing_cell_exploration.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation/ROIs_SELECTION"

python ./roi_selection_missing_cell_exploration.py