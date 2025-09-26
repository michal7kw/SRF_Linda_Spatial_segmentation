# Summary: 2_hippo_analysis.py

This script performs **spatial visualization analysis of reclustered hippocampus data** using pre-computed cluster 0 reclustering results. 
It creates individual spatial and UMAP plots specifically for the newly identified subclusters (0.0, 0.1, 0.2, etc.) that were generated from the original cluster 0 cells.

## Outputs

### Individual Subcluster Plots
- **File Pattern**: `sample_{sample_id}_subcluster_{cluster_id}.png`
- **Content**: Side-by-side UMAP and spatial plots for each reclustered subcluster
- **Examples**: 
  - `sample_sample1_subcluster_0_0.png` (for subcluster 0.0)
  - `sample_sample1_subcluster_0_1.png` (for subcluster 0.1)
  - `sample_sample2_subcluster_0_0.png` (for subcluster 0.0 in sample 2)

### Overview Plots
- **umap_overview_cluster0_highlighted.png** - UMAP showing all original cluster 0 cells highlighted in red

### Analysis Reports
- **reclustered_analysis_summary.txt** - Comprehensive summary of reclustering results per sample
        
### Visualizations

#### Leiden cluster 0
![Combined Analysis](umap_overview_cluster0_highlighted.png)