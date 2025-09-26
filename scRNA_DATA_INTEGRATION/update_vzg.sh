#!/bin/bash
#SBATCH --job-name=update_vzg
#SBATCH --output=logs/update_vzg.out
#SBATCH --error=logs/update_vzg.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=84G
#SBATCH --time=72:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# --- Step 1: Environment Setup ---
echo "--- Setting up environment ---"

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

PROJ_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"

vpt update-vzg \
    --input-vzg "${PROJ_DIR}/DATA/p30-E165/R4/data.vzg2" \
    --input-boundaries "${PROJ_DIR}/DATA/p30-E165/R4/cell_boundaries.parquet" \
    --input-entity-by-gene "${PROJ_DIR}/DATA_INTEGRATION/p30_R4_tangram_integrated.parquet" \
    --output-vzg "${PROJ_DIR}/DATA_INTEGRATION/p30_R4_tangram_updated.vzg" \
    --temp-path "${PROJ_DIR}/DATA_INTEGRATION/temp" \
    --overwrite

# echo "--- Starting VZG update ---"
# vpt --processes 4 update-vzg \
#     --input-vzg "${DATA_PATH}/${REGION}/data.vzg2" \
#     --input-boundaries "${ROI_OUTPUT_DIR}/cellpose2_micron_space.parquet" \
#     --input-entity-by-gene "${ROI_OUTPUT_DIR}/cell_by_gene.csv" \
#     --input-metadata "${ROI_OUTPUT_DIR}/cell_metadata.csv" \
#     --output-vzg "${ROI_OUTPUT_DIR}/${SAMPLE}_roi_resegmented.vzg2" \
#     --temp-path "${ROI_OUTPUT_DIR}/vzg_build_temp" \
#     --overwrite


