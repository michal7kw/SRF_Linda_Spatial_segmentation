#!/bin/bash
#SBATCH --job-name=markers_extended
#SBATCH --output=logs/markers_extended.out
#SBATCH --error=logs/markers_extended.err
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

cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation/MERGE_SPATIAL_DATA"

python ./marker_genes_per_sample_extended.py