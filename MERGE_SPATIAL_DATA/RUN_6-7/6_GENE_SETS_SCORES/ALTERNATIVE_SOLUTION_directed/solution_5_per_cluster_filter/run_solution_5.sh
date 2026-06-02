#!/bin/bash
#SBATCH --job-name=sol5_dir_cluster
#SBATCH --output=../logs/solution_5_directional.out
#SBATCH --error=../logs/solution_5_directional.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Solution 5: Per-Cluster Gene Filtering with Directional Gene Sets
# THE RECOMMENDED SOLUTION for dilution problems!
# Each cluster scored on genes actually expressed in that cluster

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "Starting Solution 5: Per-Cluster Filtering (Directional) ⭐ RECOMMENDED"
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
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_5_per_cluster_filter/results

# Run the analysis
echo "Running Solution 5 with directional gene sets..."
echo ""
echo "This is THE SOLUTION for cell-type-specific pathway components!"
echo ""
echo "For each cluster, this will:"
echo "  1. Identify which upregulated genes are expressed"
echo "  2. Identify which downregulated genes are expressed"
echo "  3. Score ONLY the expressed genes per direction"
echo "  4. Calculate aggregate = (up score) - (down score)"
echo ""
echo "Different clusters will be scored on different gene subsets!"
echo ""

python MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_5_per_cluster_filter/gene_set_scoring_per_cluster_filter_directional.py

echo "================================================================"
echo "Solution 5 (Directional) completed! Check the results in:"
echo "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_5_per_cluster_filter/results/"
echo ""
echo "Key files:"
echo "  - per_cluster_filtering_stats.csv (CRITICAL - shows which genes used per cluster)"
echo "  - heatmap_up_per_cluster.png"
echo "  - heatmap_down_per_cluster.png"
echo "  - heatmap_aggregate_per_cluster.png"
echo "  - *_up_per_cluster_sample_*.png (spatial plots)"
echo "  - *_down_per_cluster_sample_*.png"
echo "  - *_aggregate_per_cluster_sample_*.png"
echo ""
echo "NEXT STEPS:"
echo "  1. Review per_cluster_filtering_stats.csv to see gene usage"
echo "  2. Compare aggregate score heatmap across clusters"
echo "  3. Validate biological interpretation with domain knowledge"
echo "================================================================"
