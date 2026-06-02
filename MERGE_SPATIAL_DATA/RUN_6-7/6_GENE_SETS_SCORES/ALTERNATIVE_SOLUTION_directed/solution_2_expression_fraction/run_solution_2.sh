#!/bin/bash
#SBATCH --job-name=sol2_dir_fraction
#SBATCH --output=../logs/solution_2_directional.out
#SBATCH --error=../logs/solution_2_directional.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Solution 2: Expression Fraction Metrics with Directional Gene Sets
# Calculates scores AND tracks expression fraction per cluster
# Diagnostic tool to identify dilution vs true biological signal

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "Starting Solution 2: Expression Fraction Metrics (Directional)..."
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
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_2_expression_fraction/results

# Run the analysis
echo "Running Solution 2 with directional gene sets..."
echo "This will calculate scores AND expression fractions for:"
echo "  - Upregulated genes"
echo "  - Downregulated genes"
echo "  - Aggregate scores"
echo ""
echo "Diagnostic plots will help identify dilution effects"
echo ""

python MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_2_expression_fraction/gene_set_scoring_with_fraction_directional.py

echo "================================================================"
echo "Solution 2 (Directional) completed! Check the results in:"
echo "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_2_expression_fraction/results/"
echo ""
echo "Key files:"
echo "  - expression_fraction_metrics.csv"
echo "  - diagnostic_*_up.png (score vs fraction for upregulated genes)"
echo "  - diagnostic_*_down.png (score vs fraction for downregulated genes)"
echo ""
echo "INTERPRETATION:"
echo "  Low score + Low fraction = Dilution problem"
echo "  Low score + High fraction = True biological downregulation"
echo "================================================================"
