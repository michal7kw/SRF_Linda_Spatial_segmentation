#!/bin/bash
#SBATCH --job-name=sol1_global_filter
#SBATCH --output=logs/solution_1_global_filter.out
#SBATCH --error=logs/solution_1_global_filter.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Solution 1: Global Expression Filtering
# This script removes genes with low/no expression across the entire dataset

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "Starting Solution 1: Global Expression Filtering..."
echo "================================================================"

# Set the working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation

# Check if required files exist
echo "Checking required files..."

if [ ! -f "MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad" ]; then
    echo "Error: Merged spatial data file not found!"
    exit 1
fi

if [ ! -f "DATA/GENE_LISTS/output_combined/import_ready.csv" ]; then
    echo "Error: Gene sets file not found!"
    exit 1
fi

echo "All required files found."

# Create logs directory
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_1_global_filter/logs

# Create output directory
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_1_global_filter/results

# Run the analysis
echo "Running Solution 1: Global Expression Filtering..."
python MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_1_global_filter/gene_set_scoring_global_filter.py

echo "================================================================"
echo "Solution 1 completed! Check the results in:"
echo "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_1_global_filter/results/"
echo ""
echo "Key files:"
echo "  - global_filtering_statistics.csv"
echo "  - cluster_scores_heatmap_global_filter.png"
echo "  - *_global_filter_sample_*.png (spatial plots)"
echo "================================================================"
