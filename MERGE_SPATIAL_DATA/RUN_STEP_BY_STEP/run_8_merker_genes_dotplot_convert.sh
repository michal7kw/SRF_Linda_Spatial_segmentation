#!/bin/bash
#SBATCH --job-name=run_8_merker_genes_dotplot_convert
#SBATCH --output=logs/run_8_merker_genes_dotplot_convert.out
#SBATCH --error=logs/run_8_merker_genes_dotplot_convert.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Cluster-specific Temporal Markers Analysis Pipeline
# This script runs the cluster-specific temporal markers analysis for spatial transcriptomics data
# Compares ctrl vs mut samples within leiden_new clusters at corresponding time points
# Only analyzes clusters with at least 20 cells in both conditions

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment vpt..."
conda activate seurat_full2

# Set the working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation

# Run the analysis
echo "Running cluster-specific temporal markers analysis..."
python MERGE_SPATIAL_DATA/SCRIPTS/convert_h5ad_to_rds.py

