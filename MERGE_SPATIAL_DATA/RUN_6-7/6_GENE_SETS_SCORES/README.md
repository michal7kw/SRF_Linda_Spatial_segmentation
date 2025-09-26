# Gene Set Scoring and Spatial Visualization

This directory contains code to perform gene set scoring at the single cell level and create spatial visualization plots showing gene set scores as red color intensity across spatial transcriptomics samples.

## Files Description

### Scripts
- `simple_gene_set_spatial_plots.py` - Main Python script for gene set scoring and spatial visualization
- `gene_set_scoring_spatial_visualization.py` - Alternative comprehensive implementation with additional features
- `run_gene_set_analysis.sh` - Bash script to run the analysis pipeline

### Input Requirements
- **Spatial Data**: `../merged_spatial_data_hippo_inclusive.h5ad` - Merged spatial transcriptomics data
- **Gene Sets**: `../../GENE_LISTS/output_combined/import_ready.csv` - Gene sets file with format:
  - Column 1: Gene set name
  - Column 2: Gene name  
  - Columns 3-4: Ignored (additional metadata)

## Usage

### Quick Start
```bash
# Make the script executable
chmod +x run_gene_set_analysis.sh

# Run the analysis
./run_gene_set_analysis.sh
```

### Manual Execution
```bash
# Navigate to project root
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation

# Run the Python script directly
python MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/simple_gene_set_spatial_plots.py
```

## How It Works

1. **Data Loading**: Loads the merged spatial transcriptomics data and gene sets
2. **Gene Set Scoring**: Uses scanpy's `score_genes` function to calculate gene set scores for each cell
3. **Spatial Plotting**: Creates spatial scatter plots with red color intensity representing gene set scores
4. **Output Generation**: Saves individual plots for each gene set and sample combination

## Output

The analysis generates:

### Spatial Plots
- `spatial_plots/` directory containing:
  - Individual spatial plots for each gene set and sample combination
  - Format: `{gene_set_name}_sample_{sample_id}.png`
  - Red color intensity represents gene set score magnitude

### Summary Statistics
- `gene_set_scores_summary.csv` - Comprehensive statistics including:
  - Mean, standard deviation, min/max scores per gene set and sample
  - Number of genes found vs. total genes in each set
  - Cell counts per sample

## Gene Set Scoring Method

The script uses scanpy's `score_genes` function which:
- Calculates the average expression of genes in the gene set
- Subtracts the average expression of a control gene set (randomly selected genes with similar expression)
- Results in a normalized score representing gene set activity

## Visualization Features

- **Color Scheme**: Uses 'Reds' colormap for intuitive interpretation
- **Dynamic Range**: Uses 1st and 99th percentiles for color scaling to avoid outlier effects
- **High Resolution**: Saves plots at 300 DPI for publication quality
- **Informative Titles**: Includes gene set name, sample ID, and score range

## Troubleshooting

### Common Issues
- **Missing Genes**: Warning messages will show genes not found in the dataset
- **Empty Gene Sets**: Gene sets with no matching genes will be skipped
- **Memory Issues**: For large datasets, consider processing subsets of gene sets

### Dependencies
Ensure these packages are installed:
- scanpy
- squidpy  
- pandas
- numpy
- matplotlib

## Example Gene Sets Included

The gene sets cover various biological processes:
- Neural stem cells
- Neuroblast markers
- Synaptic proteins
- Axonal guidance
- Cell adhesion
- Epilepsy-related genes
- GABA interneurons
- And more...

For questions or issues, check the log output or examine the summary statistics file for insights into scoring success rates.