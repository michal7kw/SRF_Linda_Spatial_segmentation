#!/bin/bash
#SBATCH --job-name=sol1_dir_global
#SBATCH --output=../logs/solution_1_directional.out
#SBATCH --error=../logs/solution_1_directional.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Solution 1: Global Expression Filtering with Directional Gene Sets
# Removes genes with low/no expression globally before scoring
# Calculates positive, negative, and aggregate scores

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "Starting Solution 1: Global Expression Filtering (Directional)..."
echo "================================================================"

# Set the working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation

# Check if required files exist
echo "Checking required files..."

# if [ ! -f "MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad" ]; then
#     echo "Error: Merged spatial data file not found!"
#     exit 1
# fi

# if [ ! -f "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/ANALYSIS/Enrichment/GSEA/GENE_LIST/gene_lists.csv" ]; then
#     echo "Error: Directional gene sets file not found!"
#     exit 1
# fi

echo "All required files found."

# Create output directory
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_1_global_filter/results

# Run the analysis
echo "Running Solution 1 with directional gene sets..."
echo "This will calculate:"
echo "  - Positive scores (upregulated genes)"
echo "  - Negative scores (downregulated genes)"
echo "  - Aggregate scores (positive - negative)"
echo ""

python MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_1_global_filter/gene_set_scoring_global_filter_directional.py

echo "================================================================"
echo "Solution 1 (Directional) completed! Check the results in:"
echo "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_1_global_filter/results/"
echo ""
echo "Key files:"
echo "  - filtering_statistics.csv"
echo "  - heatmap_up_global_filter.png"
echo "  - heatmap_down_global_filter.png"
echo "  - heatmap_aggregate_global_filter.png"
echo "  - *_up_global_filter_sample_*.png (spatial plots - positive scores)"
echo "  - *_down_global_filter_sample_*.png (spatial plots - negative scores)"
echo "  - *_aggregate_global_filter_sample_*.png (spatial plots - aggregate scores)"
echo "================================================================"
