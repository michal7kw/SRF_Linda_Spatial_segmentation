#!/bin/bash
#SBATCH --job-name=all_dir_solutions
#SBATCH --output=logs/run_all_directional_solutions.out
#SBATCH --error=logs/run_all_directional_solutions.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Master Script: Run All Directional Alternative Solutions Sequentially
# This runs all 3 directional solutions for comprehensive analysis
# Each solution calculates positive, negative, and aggregate scores

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "================================================================"
echo "RUNNING ALL DIRECTIONAL GENE SET SCORING SOLUTIONS"
echo "================================================================"
echo ""
echo "Using directional gene lists with up/down regulation information"
echo ""
echo "This will run 3 solutions sequentially:"
echo "  1. Global Expression Filtering (Directional)"
echo "  2. Expression Fraction Metrics (Directional) - DIAGNOSTIC"
echo "  5. Per-Cluster Gene Filtering (Directional) - ⭐ RECOMMENDED"
echo ""
echo "Each solution calculates three score types:"
echo "  - Positive scores (upregulated genes)"
echo "  - Negative scores (downregulated genes)"
echo "  - Aggregate scores (positive - negative)"
echo ""
echo "Estimated total runtime: 2-3 hours"
echo "================================================================"
echo ""

# Set the working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation

# Check if required files exist
# echo "Checking required files..."
# if [ ! -f "MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad" ]; then
#     echo "Error: Merged spatial data file not found!"
#     exit 1
# fi

# if [ ! -f "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/ANALYSIS/Enrichment/GSEA/GENE_LIST/gene_lists.csv" ]; then
#     echo "Error: Directional gene sets file not found!"
#     exit 1
# fi

echo "All required files found."
echo ""

# Create main logs directory
mkdir -p MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/logs

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

# Solution 1: Global Expression Filtering (Directional)
run_solution 1 \
    "Global Expression Filtering (Directional)" \
    "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_1_global_filter" \
    "gene_set_scoring_global_filter_directional.py"

# Solution 2: Expression Fraction Metrics (Directional)
run_solution 2 \
    "Expression Fraction Metrics (Directional) - DIAGNOSTIC" \
    "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_2_expression_fraction" \
    "gene_set_scoring_with_fraction_directional.py"

# Solution 5: Per-Cluster Gene Filtering (Directional) - RECOMMENDED
run_solution 5 \
    "Per-Cluster Gene Filtering (Directional) - ⭐ RECOMMENDED" \
    "MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/solution_5_per_cluster_filter" \
    "gene_set_scoring_per_cluster_filter_directional.py"

echo "================================================================"
echo "ALL DIRECTIONAL SOLUTIONS COMPLETED!"
echo "================================================================"
echo "Completion time: $(date)"
echo ""
echo "Results are available in each solution's results/ directory:"
echo "  - solution_1_global_filter/results/"
echo "  - solution_2_expression_fraction/results/"
echo "  - solution_5_per_cluster_filter/results/"
echo ""
echo "SCORE TYPES GENERATED:"
echo "  • Positive scores - Based on upregulated genes (Reds colormap)"
echo "  • Negative scores - Based on downregulated genes (Blues colormap)"
echo "  • Aggregate scores - Combined: positive - negative (coolwarm colormap)"
echo ""
echo "RECOMMENDED ANALYSIS WORKFLOW:"
echo "  1. START HERE: Solution 5 (Per-Cluster Filter)"
echo "     - Review: per_cluster_filtering_stats.csv"
echo "     - Check: Which genes used per cluster per direction"
echo "     - View: heatmap_aggregate_per_cluster.png"
echo ""
echo "  2. DIAGNOSTIC: Solution 2 (Expression Fraction)"
echo "     - Review: expression_fraction_metrics.csv"
echo "     - Check diagnostic plots for dilution patterns"
echo "     - Compare: score vs fraction for up/down genes"
echo ""
echo "  3. BASELINE: Solution 1 (Global Filter)"
echo "     - Compare with Solution 5 to see impact of per-cluster filtering"
echo "     - Review: filtering_statistics.csv"
echo ""
echo "INTERPRETATION GUIDE:"
echo "  • High aggregate score = Strong positive regulation"
echo "  • Low aggregate score = Strong negative regulation"
echo "  • Near-zero aggregate = Balanced/inactive pathway"
echo ""
echo "  • Positive score = Expression of promoter genes"
echo "  • Negative score = Expression of inhibitor genes"
echo ""
echo "For detailed interpretation, see README.md in:"
echo "  MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/ALTERNATIVE_SOLUTION_directed/"
echo "================================================================"
