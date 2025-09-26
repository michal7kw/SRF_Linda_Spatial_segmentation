#!/bin/bash
#SBATCH --job-name=merge_hippo
#SBATCH --output=logs/merge_hippo.out
#SBATCH --error=logs/merge_hippo.err
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

python ./merge_hippo.py