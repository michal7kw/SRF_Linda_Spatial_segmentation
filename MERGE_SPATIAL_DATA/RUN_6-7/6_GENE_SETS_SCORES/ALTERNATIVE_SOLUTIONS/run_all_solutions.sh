#!/bin/bash
#SBATCH --job-name=all_solutions
#SBATCH --output=logs/run_all_solutions.out
#SBATCH --error=logs/run_all_solutions.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Master Script: Run All Alternative Solutions Sequentially
# This runs all 5 solutions one after another for comprehensive analysis

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "================================================================"
echo "RUNNING ALL ALTERNATIVE GENE SET SCORING SOLUTIONS"
echo "================================================================"
echo ""
echo "This will run all 5 solutions sequentially:"
echo "  1. Global Expression Filtering"
echo "  2. Expression Fraction Metrics (DIAGNOSTIC)"
echo "  3. Cell-Type Specific Gene Sets"
echo "  4. Weighted Scoring by Expression"
echo "  5. Per-Cluster Gene Filtering (RECOMMENDED)"
echo ""
echo "Estimated total runtime: 2-3 hours"
echo "================================================================"
echo ""

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
echo ""

# Create main logs directory
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/logs

# Function to run a solution and check for errors
run_solution() {
    local solution_num=$1
    local solution_name=$2
    local solution_dir=$3
    local python_script=$4

    echo "================================================================"
    echo "SOLUTION ${solution_num}: ${solution_name}"
    echo "================================================================"
    echo "Start time: $(date)"

    # Create solution-specific directories
    mkdir -p "${solution_dir}/logs"
    mkdir -p "${solution_dir}/results"

    # Run the solution
    python "${solution_dir}/${python_script}"

    if [ $? -eq 0 ]; then
        echo "✓ Solution ${solution_num} completed successfully"
        echo "  Results in: ${solution_dir}/results/"
    else
        echo "✗ Solution ${solution_num} failed with error code $?"
        echo "  Check logs for details"
    fi

    echo "End time: $(date)"
    echo ""
}

# Solution 1: Global Expression Filtering
run_solution 1 \
    "Global Expression Filtering" \
    "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_1_global_filter" \
    "gene_set_scoring_global_filter.py"

# Solution 2: Expression Fraction Metrics
run_solution 2 \
    "Expression Fraction Metrics (DIAGNOSTIC)" \
    "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_2_expression_fraction" \
    "gene_set_scoring_with_fraction.py"

# Solution 3: Cell-Type Specific Gene Sets
run_solution 3 \
    "Cell-Type Specific Gene Sets" \
    "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_3_celltype_specific_sets" \
    "gene_set_splitter.py"

# Solution 4: Weighted Scoring by Expression
run_solution 4 \
    "Weighted Scoring by Expression" \
    "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_4_weighted_scoring" \
    "gene_set_scoring_weighted.py"

# Solution 5: Per-Cluster Gene Filtering
run_solution 5 \
    "Per-Cluster Gene Filtering (RECOMMENDED)" \
    "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTIONS/solution_5_per_cluster_filter" \
    "gene_set_scoring_per_cluster_filter.py"

echo "================================================================"
echo "ALL SOLUTIONS COMPLETED!"
echo "================================================================"
echo "Completion time: $(date)"
echo ""
echo "Results are available in each solution's results/ directory:"
echo "  - solution_1_global_filter/results/"
echo "  - solution_2_expression_fraction/results/"
echo "  - solution_3_celltype_specific_sets/results/"
echo "  - solution_4_weighted_scoring/results/"
echo "  - solution_5_per_cluster_filter/results/"
echo ""
echo "RECOMMENDED NEXT STEPS:"
echo "  1. Check Solution 2 diagnostic: expression_fraction_metrics.csv"
echo "  2. Review Solution 5 results: per_cluster_filtering_statistics.csv"
echo "  3. Compare spatial plots across solutions"
echo "  4. Validate gene usage in Solution 5"
echo ""
echo "See README.md and SOLUTIONS_SUMMARY.md for interpretation guides"
echo "================================================================"
