#!/bin/bash
#SBATCH --job-name=celltypist_annotation
#SBATCH --output=logs/celltypist_annotation.out
#SBATCH --error=logs/celltypist_annotation.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# --- Pipeline Header ---
echo "========================================================================"
echo "CELLTYPIST ANNOTATION PIPELINE FOR SPATIAL HIPPOCAMPUS DATA"
echo "========================================================================"
echo "Description: This pipeline annotates cell types in the merged spatial"
echo "             hippocampus data using the Mouse_Dentate_Gyrus model."
echo "========================================================================"

# --- Setup ---
# Create logs directory if it doesn't exist
mkdir -p logs

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/celltypist

# Change to the working directory
cd "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"

# --- Check Prerequisites ---
MERGED_DATA_FILE="MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad"
echo ""
echo "--- Checking Prerequisites ---"
if [ ! -f "${MERGED_DATA_FILE}" ]; then
    echo "ERROR: Required merged data file '${MERGED_DATA_FILE}' not found."
    echo "Please run the inclusive pipeline first (run_pipeline_inclusive.sh)"
    exit 1
else
    echo "✓ Found merged data file: ${MERGED_DATA_FILE}"
fi

# --- Install CellTypist if needed ---
echo ""
echo "--- Checking CellTypist Installation ---"
python -c "import celltypist; print('✓ CellTypist is available')" 2>/dev/null || {
    echo "Installing CellTypist..."
    pip install celltypist
    exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "ERROR: CellTypist installation failed with exit code ${exit_code}"
        exit ${exit_code}
    fi
    echo "✓ CellTypist installed successfully"
}

# --- Run CellTypist Annotation ---
echo ""
echo "--- Running CellTypist Annotation ---"
echo "Model: Mouse_Dentate_Gyrus"
echo "Input: ${MERGED_DATA_FILE}"
echo "Output: MERGE_SPATIAL_DATA/celltypist_annotation/"

python ./MERGE_SPATIAL_DATA/SCRIPTS/annotate_spatial_hippo_clusters.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "ERROR: CellTypist annotation failed with exit code ${exit_code}"
    exit ${exit_code}
fi

echo "✓ CellTypist annotation completed successfully"

# --- Verify Output Files ---
echo ""
echo "--- Verifying Output Files ---"
OUTPUT_DIR="MERGE_SPATIAL_DATA/celltypist_annotation"
ANNOTATED_FILE="${OUTPUT_DIR}/merged_spatial_hippo_annotated.h5ad"
PLOTS_DIR="${OUTPUT_DIR}/plots"

if [ -f "${ANNOTATED_FILE}" ]; then
    echo "✓ Annotated data file created: ${ANNOTATED_FILE}"
    file_size=$(ls -lh "${ANNOTATED_FILE}" | awk '{print $5}')
    echo "  File size: ${file_size}"
else
    echo "⚠ WARNING: Annotated data file not found: ${ANNOTATED_FILE}"
fi

if [ -d "${PLOTS_DIR}" ]; then
    plot_count=$(find "${PLOTS_DIR}" -name "*.png" | wc -l)
    echo "✓ Plots directory created: ${PLOTS_DIR}"
    echo "  Number of plot files: ${plot_count}"
else
    echo "⚠ WARNING: Plots directory not found: ${PLOTS_DIR}"
fi

# Check for summary files
summary_files=$(find "${OUTPUT_DIR}" -name "*.csv" | wc -l)
echo "✓ Summary files created: ${summary_files}"

# --- Pipeline Completion Summary ---
echo ""
echo "========================================================================"
echo "CELLTYPIST ANNOTATION PIPELINE COMPLETED SUCCESSFULLY"
echo "========================================================================"
echo "Results:"
echo "  - Annotated data: ${ANNOTATED_FILE}"
echo "  - Visualization plots: ${PLOTS_DIR}/"
echo "  - Summary statistics: ${OUTPUT_DIR}/*.csv"
echo ""
echo "Key outputs:"
echo "  ✓ Cell type annotations using Mouse_Dentate_Gyrus model"
echo "  ✓ UMAP and spatial visualizations"
echo "  ✓ Confidence scores and statistics"
echo "  ✓ Comparison with original Leiden clustering"
echo "  ✓ Sample-wise cell type distributions"
echo ""
echo "Next steps:"
echo "  - Review the generated plots in ${PLOTS_DIR}/"
echo "  - Check cell type summaries in the CSV files"
echo "  - Use the annotated data for downstream analysis"
echo "========================================================================"