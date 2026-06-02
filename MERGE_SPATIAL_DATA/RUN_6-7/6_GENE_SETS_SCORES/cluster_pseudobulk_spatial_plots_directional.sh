#!/bin/bash
#SBATCH --job-name=cluster_pseudobulk_spatial_plots_directional
#SBATCH --output=logs/cluster_pseudobulk_spatial_plots_directional.out
#SBATCH --error=logs/cluster_pseudobulk_spatial_plots_directional.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Cluster Pseudo-Bulk Gene Set Scoring and Spatial Visualization Pipeline
# This script calculates gene set scores at the cluster level (pseudo-bulk)
# and displays these scores on spatial plots

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/vpt

echo "Starting Cluster Pseudo-Bulk Gene Set Scoring Pipeline..."
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

# Create output directory and logs directory
mkdir -p MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/cluster_pseudobulk_spatial_plots_directional
mkdir -p logs

# Run the analysis
echo "Running cluster pseudo-bulk gene set scoring and spatial visualization..."
python MERGE_SPATIAL_DATA/RUN_6-7/6_GENE_SETS_SCORES/cluster_pseudobulk_spatial_plots_directional.py

echo "================================================================"
echo "Pipeline completed! Check the results in:"
echo "MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/cluster_pseudobulk_spatial_plots_directional/"
echo "================================================================"
