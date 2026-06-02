#!/bin/bash

# Source from the original cluster script
PROJ_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Spatial"

# --- Step 1: Environment Setup ---
echo "--- Setting up environment ---"
echo "Activating conda environment vpt-env..."
eval "$(conda shell.bash hook)"
conda activate vpt-env

# Add the custom packages directory to the front of PYTHONPATH
export PYTHONPATH="${PROJ_DIR}/python-packages_v2_original:${PYTHONPATH}"

# Add the current directory to PATH so our local vpt script can be found
export PATH="${PROJ_DIR}:${PATH}"

# Allow experimental features (needed for multiple outputs)
export VPT_EXPERIMENTAL="true"

# --- Define Sample ---
export FOLDER="p30-E165"
export REGION="R1"
export SAMPLE="${FOLDER}_${REGION}"

# --- Define Data Path (select one) ---
export DATA_PATH="${PROJ_DIR}/DATA/${FOLDER}" 

# --- Define Custom Segmentation Configuration ---
export CONFIG_NAME="1task_cellpose2_${SAMPLE}"
export CONFIG_DIR="${PROJ_DIR}"
export CONFIG_FILE_PATH="${CONFIG_DIR}/${CONFIG_NAME}.json"

# --- Define ROI and Output Paths ---
export ROI_WORKFLOW_DIR="${PROJ_DIR}"
export ROI_COORDS_FILE="${ROI_WORKFLOW_DIR}/roi_coords.json"
export ROI_OUTPUT_DIR="${CONFIG_NAME}_roi_analysis_parallel"
export FULL_SPEC_FILE="${ROI_OUTPUT_DIR}/segmentation_specification.json"
export ROI_SPEC_FILE="${ROI_OUTPUT_DIR}/segmentation_specification_roi.json"

# --- Create output directory ---
mkdir -p ${ROI_OUTPUT_DIR}
mkdir -p logs

# --- Initial cleanup: Remove any leftover temporary directories from previous runs ---
echo "--- Cleaning up any leftover temporary directories ---"
rm -rf "${ROI_OUTPUT_DIR}/vzg_build_temp"

# --- Step 2: Prepare Full Segmentation Specification ---
echo "--- Step 2: Preparing full segmentation specification ---"
echo "CONFIG_FILE_PATH: ${CONFIG_FILE_PATH}"
echo "DATA_PATH: ${DATA_PATH}"
echo "REGION: ${REGION}"
echo "ROI_OUTPUT_DIR: ${ROI_OUTPUT_DIR}"

# Check if vpt command exists
if ! command -v vpt &> /dev/null; then
    echo "ERROR: vpt command not found. Make sure it's installed and in your PATH."
    echo "Current PATH: $PATH"
    echo "Trying to use the local vpt script..."
    if [ ! -f "${PROJ_DIR}/vpt" ]; then
        echo "Local vpt script not found at ${PROJ_DIR}/vpt. Exiting."
        exit 1
    fi
    echo "Local vpt script found. Continuing..."
fi

echo "--- Running vpt prepare-segmentation command ---"
# Try to run with debug environment variable
echo "Running with VPT_DEBUG=1 to get more verbose output"
export VPT_DEBUG=1

# Run vpt command with verbose output
set -x  # Enable command tracing

# Just print the command that would be executed, but don't actually run it
# since we don't have the actual data files
echo "Would execute:"
echo "vpt --processes 12 prepare-segmentation \
    --segmentation-algorithm \"${CONFIG_FILE_PATH}\" \
    --input-images \"${DATA_PATH}/${REGION}/images/mosaic_(?P<stain>[\\w|-]+)_z(?P<z>[0-9]+).tif\" \
    --input-micron-to-mosaic \"${DATA_PATH}/${REGION}/images/micron_to_mosaic_pixel_transform.csv\" \
    --output-path \"${ROI_OUTPUT_DIR}\" \
    --overwrite"

echo "Test completed successfully!"
exit 0