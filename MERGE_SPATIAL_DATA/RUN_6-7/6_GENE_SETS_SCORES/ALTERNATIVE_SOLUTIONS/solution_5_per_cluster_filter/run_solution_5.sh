#!/bin/bash
#SBATCH --job-name=sol5_per_cluster
#SBATCH --output=logs/solution_5_per_cluster_filter.out
#SBATCH --error=logs/solution_5_per_cluster_filter.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Solution 5: Per-Cluster Gene Filtering ⭐ RECOMMENDED FOR REELIN PATHWAY
# This directly solves the Reelin dilution problem!
# GCs scored on receptors, CR cells scored on ligands

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "Starting Solution 5: Per-Cluster Gene Filtering ⭐"
echo "================================================================"
echo "THIS SOLUTION DIRECTLY ADDRESSES YOUR REELIN DILUTION CONCERN!"
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
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_5_per_cluster_filter/logs

# Create output directory
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_5_per_cluster_filter/results

# Run the analysis
echo "Running Solution 5: Per-Cluster Gene Filtering..."
python MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_5_per_cluster_filter/gene_set_scoring_per_cluster_filter.py

echo "================================================================"
echo "Solution 5 completed! Check the results in:"
echo "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_5_per_cluster_filter/results/"
echo ""
echo "CRITICAL FILES TO CHECK:"
echo "  - per_cluster_filtering_statistics.csv (WHICH GENES USED WHERE)"
echo "  - gene_usage_*.png (visual gene usage matrix)"
echo "  - per_cluster_scores_heatmap.png (final scores)"
echo "  - per_cluster_n_genes_heatmap.png (N genes per cluster)"
echo "  - *_per_cluster_filter_sample_*.png (spatial plots)"
echo ""
echo "VALIDATION STEPS:"
echo "  1. Open per_cluster_filtering_statistics.csv"
echo "  2. Verify GCs use Reelin receptors (Vldlr, Dab1, ApoER2)"
echo "  3. Verify CR cells use Reelin ligand (Reln)"
echo "  4. Check gene_usage_*.png for visual confirmation"
echo ""
echo "This solution eliminates dilution by scoring each cluster on its"
echo "relevant pathway components only!"
echo "================================================================"
