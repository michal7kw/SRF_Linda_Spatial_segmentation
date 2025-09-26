#!/bin/bash
#SBATCH --job-name=tangram_integration
#SBATCH --output=logs/tangram_integration.out
#SBATCH --error=logs/tangram_integration.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=84G
#SBATCH --time=72:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=cuda
#SBATCH --gpus=v100:1

# --- Step 1: Environment Setup ---
echo "--- Setting up environment ---"

echo "Initializing conda..."
eval "$(conda shell.bash hook)"
conda activate tangram

# --- Check CUDA device detection ---
echo "--- Checking CUDA device detection ---"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES}"
if command -v nvidia-smi &> /dev/null; then
    echo "nvidia-smi output:"
    nvidia-smi
else
    echo "nvidia-smi command not found. CUDA might not be installed or accessible."
fi

WORK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation/DATA_INTEGRATION"

python "${WORK_DIR}/tangram_integration.py"