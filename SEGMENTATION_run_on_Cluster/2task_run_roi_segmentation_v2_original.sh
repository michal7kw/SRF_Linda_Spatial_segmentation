#!/bin/bash
#SBATCH --job-name=2task_v2_original
#SBATCH --output=logs/2task_v2_original.out
#SBATCH --error=logs/2task_v2_original.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=84G
#SBATCH --time=72:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=cuda
#SBATCH --gpus=v100:1

# This script performs segmentation on a defined Region of Interest (ROI)
# using a custom model specified in the CONFIG json file.

# --- Step 1: Environment Setup ---
echo "--- Setting up environment ---"
# export CUDA_VISIBLE_DEVICES=""

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

# --- Check CUDA device detection ---
echo "--- Checking CUDA device detection ---"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES}"
if command -v nvidia-smi &> /dev/null; then
    echo "nvidia-smi output:"
    nvidia-smi
else
    echo "nvidia-smi command not found. CUDA might not be installed or accessible."
fi

# --- Check conda environment ---
echo "--- Checking conda environment ---"
echo "CONDA_DEFAULT_ENV: ${CONDA_DEFAULT_ENV}"
if command -v conda &> /dev/null; then
    echo "Active conda environments:"
    conda info --envs
else
    echo "conda command not found or not in PATH."
fi

PROJ_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"

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
export CONFIG_NAME="2task_cellpose2_${SAMPLE}"
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

# # --- Initial cleanup: Remove any leftover temporary directories from previous runs ---
# echo "--- Cleaning up any leftover temporary directories ---"
# rm -rf "${ROI_OUTPUT_DIR}/vzg_build_temp"

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
# echo "Running with VPT_DEBUG=1 to get more verbose output"
# export VPT_DEBUG=1

# Run vpt command with verbose output
set -x  # Enable command tracing
echo "--- Step 2: Preparing full segmentation specification ---"
vpt --processes 12 prepare-segmentation \
    --segmentation-algorithm "${CONFIG_FILE_PATH}" \
    --input-images "${DATA_PATH}/${REGION}/images/mosaic_(?P<stain>[\\w|-]+)_z(?P<z>[0-9]+).tif" \
    --input-micron-to-mosaic "${DATA_PATH}/${REGION}/images/micron_to_mosaic_pixel_transform.csv" \
    --output-path "${ROI_OUTPUT_DIR}" \
    --overwrite

# --- Step 3: Filter Spec and Get Tile Indices for ROI ---
echo "--- Step 3: Filtering spec and getting tile indices for ROI ---"
echo "FULL_SPEC_FILE: ${FULL_SPEC_FILE}"
echo "ROI_COORDS_FILE: ${ROI_COORDS_FILE}"
echo "ROI_SPEC_FILE: ${ROI_SPEC_FILE}"

ROI_TILE_INDICES=$(python ${ROI_WORKFLOW_DIR}/filter_spec_and_get_indices_v2.py \
    --input-spec "${FULL_SPEC_FILE}" \
    --input-roi "${ROI_COORDS_FILE}" \
    --transform "${DATA_PATH}/${REGION}/images/micron_to_mosaic_pixel_transform.csv" \
    --output-spec "${ROI_SPEC_FILE}")

echo "Found tile indices for ROI: ${ROI_TILE_INDICES}"

# --- Step 4: Run Segmentation on ROI Tiles ---
echo "--- Step 4: Running segmentation on ROI tiles sequentially ---"
for TILE_INDEX in ${ROI_TILE_INDICES}
do
    echo "--- Running segmentation on tile ${TILE_INDEX} ---"
    vpt run-segmentation-on-tile \
        --input-segmentation-parameters "${FULL_SPEC_FILE}" \
        --tile-index ${TILE_INDEX} \
        --overwrite
done

echo "--- Contents of result_tiles directory after segmentation ---"
ls -l "${ROI_OUTPUT_DIR}/result_tiles"

# --- Step 5: Rename Output Files for Compile Step ---
echo "--- Step 5: Renaming output files for compile step ---"

# Find all generated cell parquet files and rename them sequentially
i=0
CELL_FILES=$(ls -v "${ROI_OUTPUT_DIR}/result_tiles/cell_"*.parquet)
for f in ${CELL_FILES}
do
    echo "Renaming ${f} to cell_${i}.parquet"
    mv "${f}" "${ROI_OUTPUT_DIR}/result_tiles/cell_${i}.parquet"
    i=$((i+1))
done

# Find all generated nucleus parquet files and rename them sequentially
i=0
# The 2>/dev/null suppresses the "No such file or directory" error if no nucleus files are found.
NUCLEUS_FILES=$(ls -v "${ROI_OUTPUT_DIR}/result_tiles/nucleus_"*.parquet 2>/dev/null)
for f in ${NUCLEUS_FILES}
do
    echo "Renaming ${f} to nucleus_${i}.parquet"
    mv "${f}" "${ROI_OUTPUT_DIR}/result_tiles/nucleus_${i}.parquet"
    i=$((i+1))
done

# --- Step 6: Compile Segmentation Results ---
echo "--- Step 6: Compiling segmentation results ---"
vpt --processes 12 compile-tile-segmentation \
    --input-segmentation-parameters "${ROI_SPEC_FILE}" \
    --overwrite

# --- Step 7: Downstream Analysis ---
echo "--- Step 7: Performing downstream analysis on ROI ---"
vpt --processes 12 partition-transcripts \
    --input-boundaries "${ROI_OUTPUT_DIR}/cellpose2_micron_space.parquet" \
    --input-transcripts "${DATA_PATH}/${REGION}/detected_transcripts.csv" \
    --output-entity-by-gene "${ROI_OUTPUT_DIR}/cell_by_gene.csv" \
    --overwrite

vpt --processes 12 derive-entity-metadata \
    --input-boundaries "${ROI_OUTPUT_DIR}/cellpose2_micron_space.parquet" \
    --input-entity-by-gene "${ROI_OUTPUT_DIR}/cell_by_gene.csv" \
    --output-metadata "${ROI_OUTPUT_DIR}/cell_metadata.csv" \
    --overwrite

# --- Step 8: Update VZG File with ROI Segmentation ---
echo "--- Step 8: Updating VZG file ---"

# Enhanced cleanup: Remove any existing temporary directories
# echo "--- Cleaning up all temporary directories before VZG update ---"
# rm -rf "${ROI_OUTPUT_DIR}/vzg_build_temp"

# # Wait a moment to ensure filesystem operations complete
# sleep 2

# # Verify cleanup was successful
# if [ -d "${ROI_OUTPUT_DIR}/vzg_build_temp" ]; then
#     echo "Warning: vzg_build_temp directory still exists, attempting forced removal..."
#     chmod -R 777 "${ROI_OUTPUT_DIR}/vzg_build_temp" 2>/dev/null || true
#     rm -rf "${ROI_OUTPUT_DIR}/vzg_build_temp"
#     sleep 1
# fi

echo "--- Starting VZG update ---"
vpt --processes 4 update-vzg \
    --input-vzg "${DATA_PATH}/${REGION}/data.vzg2" \
    --input-boundaries "${ROI_OUTPUT_DIR}/cellpose2_micron_space.parquet" \
    --input-entity-by-gene "${ROI_OUTPUT_DIR}/cell_by_gene.csv" \
    --input-metadata "${ROI_OUTPUT_DIR}/cell_metadata.csv" \
    --output-vzg "${ROI_OUTPUT_DIR}/${SAMPLE}_roi_resegmented.vzg2" \
    --temp-path "${ROI_OUTPUT_DIR}/vzg_build_temp" \
    --overwrite

# --- Final cleanup: Remove temporary directories after successful completion ---
echo "--- Final cleanup: Removing temporary directories ---"
rm -rf "${ROI_OUTPUT_DIR}/vzg_build_temp"

echo "--- ROI segmentation workflow finished successfully! ---"