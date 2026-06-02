#!/bin/bash

# This script performs segmentation on a defined Region of Interest (ROI)
# using a custom model specified in the CONFIG json file.

# --- Step 1: Environment Setup ---
echo "--- Setting up environment ---"
# export CUDA_VISIBLE_DEVICES=""
# mamba activate vizgen2

# Add the parent directory of the plugins to PYTHONPATH
# Add the custom packages directory to the front of PYTHONPATH
export PYTHONPATH="/home/michal/Github/SRF_Spatial/python-packages_v2_original:${PYTHONPATH}"

# Allow experimental features (needed for multiple outputs)
export VPT_EXPERIMENTAL="true"

# --- Define Sample ---
export FOLDER="p30-E165"
export REGION="R1"
export SAMPLE="${FOLDER}_${REGION}"

# --- Define Data Path (select one) ---
export DATA_PATH="/home/michal/Github/SRF_Spatial/DATA/${FOLDER}" 
# export DATA_PATH="/beegfs/scratch/ric.broccoli/kubacki.michal/SPATIAL_data/data_${FOLDER}"

# --- Define Custom Segmentation Configuration ---
# This should be the name of your segmentation JSON file with the custom model
export CONFIG_NAME="2task_cellpose2_${SAMPLE}"
export CONFIG_DIR="./" # Assumes JSON is in the Vpt-segmentation folder
export CONFIG_FILE_NAME="${CONFIG_DIR}/${CONFIG_NAME}.json"

# --- Define ROI and Output Paths ---
export ROI_WORKFLOW_DIR="."
export ROI_COORDS_FILE="${ROI_WORKFLOW_DIR}/roi_coords.json"
export ROI_OUTPUT_DIR="${CONFIG_NAME}_roi_analysis"
export FULL_SPEC_FILE="${ROI_OUTPUT_DIR}/segmentation_specification.json"
export ROI_SPEC_FILE="${ROI_OUTPUT_DIR}/segmentation_specification_roi.json"

# --- Create output directory ---
mkdir -p ${ROI_OUTPUT_DIR}
mkdir -p logs

# # --- Initial cleanup: Remove any leftover temporary directories from previous runs ---
# echo "--- Cleaning up any leftover temporary directories ---"
# rm -rf "${ROI_OUTPUT_DIR}/vzg_build_temp_v2"
# rm -rf "${ROI_OUTPUT_DIR}/vzg_build_temp"

# --- Step 2: Prepare Full Segmentation Specification ---
echo "--- Step 2: Preparing full segmentation specification ---"
vpt --processes 12 prepare-segmentation \
    --segmentation-algorithm "$(pwd)/${CONFIG_FILE_NAME}" \
    --input-images "${DATA_PATH}/${REGION}/images/mosaic_(?P<stain>[\\w|-]+)_z(?P<z>[0-9]+).tif" \
    --input-micron-to-mosaic "${DATA_PATH}/${REGION}/images/micron_to_mosaic_pixel_transform.csv" \
    --output-path "${ROI_OUTPUT_DIR}" \
    --overwrite

# --- Step 3: Filter Spec and Get Tile Indices for ROI ---
echo "--- Step 3: Filtering spec and getting tile indices for ROI ---"
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

--- Step 8: Update VZG File with ROI Segmentation ---
echo "--- Step 8: Updating VZG file ---"

# Enhanced cleanup: Remove any existing temporary directories
# echo "--- Cleaning up all temporary directories before VZG update ---"
# rm -rf "${ROI_OUTPUT_DIR}/vzg_build_temp_v2"
# rm -rf "${ROI_OUTPUT_DIR}/vzg_build_temp"

# # Wait a moment to ensure filesystem operations complete
# sleep 2

# # Verify cleanup was successful
# if [ -d "${ROI_OUTPUT_DIR}/vzg_build_temp_v2" ]; then
#     echo "Warning: vzg_build_temp_v2 directory still exists, attempting forced removal..."
#     chmod -R 777 "${ROI_OUTPUT_DIR}/vzg_build_temp_v2" 2>/dev/null || true
#     rm -rf "${ROI_OUTPUT_DIR}/vzg_build_temp_v2"
#     sleep 1
# fi

echo "--- Starting VZG update ---"
vpt --processes 4 update-vzg \
    --input-vzg "${DATA_PATH}/${REGION}/data.vzg2" \
    --input-boundaries "${ROI_OUTPUT_DIR}/cellpose2_micron_space.parquet" \
    --input-entity-by-gene "${ROI_OUTPUT_DIR}/cell_by_gene.csv" \
    --input-metadata "${ROI_OUTPUT_DIR}/cell_metadata.csv" \
    --output-vzg "${ROI_OUTPUT_DIR}/${SAMPLE}_roi_resegmented.vzg2" \
    --temp-path "${ROI_OUTPUT_DIR}/vzg_build_temp_v2" \
    --overwrite

# --- Final cleanup: Remove temporary directories after successful completion ---
# echo "--- Final cleanup: Removing temporary directories ---"
# rm -rf "${ROI_OUTPUT_DIR}/vzg_build_temp_v2"
# rm -rf "${ROI_OUTPUT_DIR}/vzg_build_temp"

echo "--- ROI segmentation workflow finished successfully! ---"