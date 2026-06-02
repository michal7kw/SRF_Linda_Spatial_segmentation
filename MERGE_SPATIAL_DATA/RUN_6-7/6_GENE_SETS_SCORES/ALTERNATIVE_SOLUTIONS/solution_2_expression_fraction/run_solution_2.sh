#!/bin/bash
#SBATCH --job-name=sol2_expr_fraction
#SBATCH --output=logs/solution_2_expression_fraction.out
#SBATCH --error=logs/solution_2_expression_fraction.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Solution 2: Expression Fraction Metrics
# DIAGNOSTIC TOOL - Start here to identify if dilution is a problem!

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "Starting Solution 2: Expression Fraction Metrics (DIAGNOSTIC)..."
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
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_2_expression_fraction/logs

# Create output directory
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_2_expression_fraction/results

# Run the analysis
echo "Running Solution 2: Expression Fraction Metrics..."
python MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_2_expression_fraction/gene_set_scoring_with_fraction.py

echo "================================================================"
echo "Solution 2 completed! Check the results in:"
echo "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_2_expression_fraction/results/"
echo ""
echo "KEY DIAGNOSTIC FILES:"
echo "  - expression_fraction_metrics.csv (CHECK THIS FIRST!)"
echo "  - diagnostic_scatter_*.png (score vs fraction plots)"
echo "  - expression_fraction_heatmap.png"
echo "  - combined_score_fraction_heatmap.png"
echo ""
echo "INTERPRETATION:"
echo "  If Fraction_Expressed < 0.3 → SEVERE dilution → Use Solution 5"
echo "  If Fraction_Expressed > 0.5 → MINOR dilution → Solution 1 may suffice"
echo "================================================================"
