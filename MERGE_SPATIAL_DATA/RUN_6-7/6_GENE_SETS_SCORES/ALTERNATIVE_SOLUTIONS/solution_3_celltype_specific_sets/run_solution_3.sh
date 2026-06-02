#!/bin/bash
#SBATCH --job-name=sol3_gene_splitter
#SBATCH --output=logs/solution_3_celltype_specific_sets.out
#SBATCH --error=logs/solution_3_celltype_specific_sets.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Solution 3: Cell-Type Specific Gene Sets
# This script analyzes expression patterns and suggests gene set splits

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "Starting Solution 3: Cell-Type Specific Gene Sets..."
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
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_3_celltype_specific_sets/logs

# Create output directory
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_3_celltype_specific_sets/results

# Run the analysis
echo "Running Solution 3: Gene Set Splitter..."
python MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_3_celltype_specific_sets/gene_set_splitter.py

echo "================================================================"
echo "Solution 3 completed! Check the results in:"
echo "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_3_celltype_specific_sets/results/"
echo ""
echo "KEY FILES TO REVIEW:"
echo "  - suggested_gene_set_splits_import_ready.csv (READY TO USE!)"
echo "  - gene_set_split_summary.csv (human-readable)"
echo "  - expression_heatmap_*.png (expression patterns)"
echo "  - gene_clustering_*.png (dendrograms)"
echo ""
echo "NEXT STEPS:"
echo "  1. Review suggested splits in gene_set_split_summary.csv"
echo "  2. Examine heatmaps to understand expression patterns"
echo "  3. Manually refine splits based on biological knowledge"
echo "  4. Use suggested_gene_set_splits_import_ready.csv for scoring"
echo "================================================================"
