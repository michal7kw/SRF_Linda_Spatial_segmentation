#!/bin/bash

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate vpt-env

# Set up environment variables
export PYTHONPATH="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Spatial/python-packages_v2_original:${PYTHONPATH}"
export PATH="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Spatial:${PATH}"
export VPT_EXPERIMENTAL="true"
export VPT_DEBUG=1

# Define paths
PROJ_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Spatial"
FOLDER="p30-E165"
REGION="R1"
SAMPLE="${FOLDER}_${REGION}"
DATA_PATH="${PROJ_DIR}/DATA/${FOLDER}"
CONFIG_NAME="1task_cellpose2_${SAMPLE}"
CONFIG_FILE_PATH="${PROJ_DIR}/${CONFIG_NAME}.json"
ROI_OUTPUT_DIR="${CONFIG_NAME}_roi_analysis_parallel"

# Create output directory
mkdir -p ${ROI_OUTPUT_DIR}

# Run the prepare-segmentation command with verbose output
echo "Running prepare-segmentation with verbose output..."
vpt --verbose --processes 1 prepare-segmentation \
    --segmentation-algorithm "${CONFIG_FILE_PATH}" \
    --input-images "${DATA_PATH}/${REGION}/images/mosaic_(?P<stain>[\\w|-]+)_z(?P<z>[0-9]+).tif" \
    --input-micron-to-mosaic "${DATA_PATH}/${REGION}/images/micron_to_mosaic_pixel_transform.csv" \
    --output-path "${ROI_OUTPUT_DIR}" \
    --overwrite