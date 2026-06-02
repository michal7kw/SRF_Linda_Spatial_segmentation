#!/bin/bash
#SBATCH --job-name=sol4_weighted_scoring
#SBATCH --output=logs/solution_4_weighted_scoring.out
#SBATCH --error=logs/solution_4_weighted_scoring.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Solution 4: Weighted Scoring by Expression
# Scores only genes expressed above threshold (complete dilution fix)

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "Starting Solution 4: Weighted Scoring by Expression..."
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
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_4_weighted_scoring/logs

# Create output directory
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_4_weighted_scoring/results

# Run the analysis
echo "Running Solution 4: Weighted Scoring..."
python MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_4_weighted_scoring/gene_set_scoring_weighted.py

echo "================================================================"
echo "Solution 4 completed! Check the results in:"
echo "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_4_weighted_scoring/results/"
echo ""
echo "Key files:"
echo "  - weighted_scoring_statistics.csv (scores and N genes used)"
echo "  - weighted_scores_heatmap.png"
echo "  - n_genes_used_heatmap.png"
echo "  - diagnostic_weighted_*.png"
echo "  - *_weighted_sample_*.png (spatial plots)"
echo ""
echo "IMPORTANT: Check N_Genes_Used column - scores based on <3 genes may be unreliable"
echo "================================================================"
