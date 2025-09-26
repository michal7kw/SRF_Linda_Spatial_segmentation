#!/bin/bash

# Activate conda environment
echo "Activating conda environment vpt-env..."
eval "$(conda shell.bash hook)"
conda activate vpt-env

# Set up environment variables
PROJ_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Spatial"
export PYTHONPATH="${PROJ_DIR}/python-packages_v2_original:${PYTHONPATH}"
export PATH="${PROJ_DIR}:${PATH}"

# Test if vpt command is found
if command -v vpt &> /dev/null; then
    echo "vpt command found!"
    echo "Path to vpt: $(which vpt)"
    echo "Testing vpt help output:"
    vpt --help
else
    echo "ERROR: vpt command not found."
    echo "Current PATH: $PATH"
    echo "Current PYTHONPATH: $PYTHONPATH"
    exit 1
fi