#!/bin/bash
#SBATCH --job-name=full_resegmentation
#SBATCH --output=logs/full_resegmentation_%A_%a.out
#SBATCH --error=logs/full_resegmentation_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=84G
#SBATCH --time=72:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=cuda
#SBATCH --gpus=v100:1
#SBATCH --array=0-9

# This script performs full image resegmentation for multiple samples using a SLURM job array.
# Each job in the array processes one sample.

# --- Step 1: Environment Setup ---
echo "--- Setting up environment ---"
eval "$(conda shell.bash hook)"
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

PROJ_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
export PYTHONPATH="${PROJ_DIR}/python-packages_v2_original:${PYTHONPATH}"
export PATH="${PROJ_DIR}:${PATH}"
export VPT_EXPERIMENTAL="true"
export PYTHONWARNINGS="ignore:The 'openvino.runtime' module is deprecated"

# --- Define Samples ---
SAMPLES=(
    "p0-p7/R1"
    "p0-p7/R2"
    "p0-p7/R3"
    "p0-p7/R4"
    "p0-p7/R5"
    "p0-p7/R6"
    "p30-E165/R1"
    "p30-E165/R2"
    "p30-E165/R3"
    "p30-E165/R4"
)

# --- Get current sample from SLURM array task ID ---
CURRENT_SAMPLE_PATH=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
FOLDER=$(dirname ${CURRENT_SAMPLE_PATH})
REGION=$(basename ${CURRENT_SAMPLE_PATH})
SAMPLE_NAME="${FOLDER//\//-}_${REGION}"

echo "--- Starting processing for sample: ${SAMPLE_NAME} ---"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Folder: ${FOLDER}, Region: ${REGION}"

# --- Define Data and Output Paths ---
DATA_PATH="${PROJ_DIR}/DATA/${FOLDER}"
OUTPUT_DIR="${PROJ_DIR}/Full_resegmentation_results/${SAMPLE_NAME}"
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# --- Define Custom Segmentation Configuration ---
CONFIG_NAME="2task_cellpose2"
CONFIG_FILE_PATH="${PROJ_DIR}/SEGMENTATION_configs/${CONFIG_NAME}.json"

echo "Data Path: ${DATA_PATH}"
echo "Output Dir: ${OUTPUT_DIR}"
echo "Config File: ${CONFIG_FILE_PATH}"

if [ ! -f "${CONFIG_FILE_PATH}" ]; then
    echo "ERROR: Config file not found at ${CONFIG_FILE_PATH}. Exiting."
    exit 1
fi

# --- Step 2: Run Full Segmentation ---
echo "--- Step 2: Running full segmentation ---"
vpt --processes 12 run-segmentation \
    --segmentation-algorithm "${CONFIG_FILE_PATH}" \
    --input-images "${DATA_PATH}/${REGION}/images/mosaic_(?P<stain>[\\w|-]+)_z(?P<z>[0-9]+).tif" \
    --input-micron-to-mosaic "${DATA_PATH}/${REGION}/images/micron_to_mosaic_pixel_transform.csv" \
    --output-path "${OUTPUT_DIR}" \
    --overwrite

# The output boundary file from run-segmentation is assumed to be named based on the task name in the config,
# which we assume results in 'cellpose2_micron_space.parquet'.
BOUNDARY_FILE="${OUTPUT_DIR}/cellpose2_micron_space.parquet"
if [ ! -f "${BOUNDARY_FILE}" ]; then
    # Fallback to default name if assumption is wrong
    BOUNDARY_FILE="${OUTPUT_DIR}/segmentation_boundaries.parquet"
fi

if [ ! -f "${BOUNDARY_FILE}" ]; then
    echo "ERROR: Segmentation boundary file not found after running segmentation. Exiting."
    exit 1
fi
echo "Boundary file found: ${BOUNDARY_FILE}"


# --- Step 3: Downstream Analysis ---
echo "--- Step 3: Performing downstream analysis ---"
vpt --processes 12 partition-transcripts \
    --input-boundaries "${BOUNDARY_FILE}" \
    --input-transcripts "${DATA_PATH}/${REGION}/detected_transcripts.csv" \
    --output-entity-by-gene "${OUTPUT_DIR}/cell_by_gene.csv" \
    --overwrite

vpt --processes 12 derive-entity-metadata \
    --input-boundaries "${BOUNDARY_FILE}" \
    --input-entity-by-gene "${OUTPUT_DIR}/cell_by_gene.csv" \
    --output-metadata "${OUTPUT_DIR}/cell_metadata.csv" \
    --overwrite

# --- Step 4: Update VZG File ---
echo "--- Step 4: Updating VZG file ---"
rm -rf "${OUTPUT_DIR}/vzg_build_temp"

vpt --processes 4 update-vzg \
    --input-vzg "${DATA_PATH}/${REGION}/data.vzg2" \
    --input-boundaries "${BOUNDARY_FILE}" \
    --input-entity-by-gene "${OUTPUT_DIR}/cell_by_gene.csv" \
    --input-metadata "${OUTPUT_DIR}/cell_metadata.csv" \
    --output-vzg "${OUTPUT_DIR}/${SAMPLE_NAME}_resegmented.vzg2" \
    --temp-path "${OUTPUT_DIR}/vzg_build_temp" \
    --overwrite

# --- Final cleanup ---
echo "--- Final cleanup: Removing temporary directories ---"
rm -rf "${OUTPUT_DIR}/vzg_build_temp"

echo "--- Full resegmentation workflow for ${SAMPLE_NAME} finished successfully! ---"