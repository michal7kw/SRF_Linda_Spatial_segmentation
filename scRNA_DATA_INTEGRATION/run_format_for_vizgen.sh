#!/bin/bash
#SBATCH --job-name=format_for_vizgen
#SBATCH --output=logs/format_for_vizgen.out
#SBATCH --error=logs/format_for_vizgen.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=84G
#SBATCH --time=1:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# --- Step 1: Environment Setup ---
echo "--- Setting up environment ---"

echo "Initializing conda..."
eval "$(conda shell.bash hook)"
conda activate tangram

WORK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation/DATA_INTEGRATION"

python "${WORK_DIR}/format_for_vizgen.py"