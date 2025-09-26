Key Features Implemented:

1. 1.
   Data Loading & Reclustering : The script loads the merged spatial data and performs cluster 0 reclustering for each sample using the same methodology as the marker genes script
2. 2.
   Subcluster Generation : For samples containing cluster 0, the script reclusters those cells into subclusters labeled as "0.0", "0.1", "0.2", etc.
3. 3.
   Individual Plot Generation : Creates both UMAP and spatial plots specifically for new subclusters, including:
   
   - Reclustered subclusters (0.0, 0.1, 0.2, etc.) for samples that had cluster 0
   - Special color coding for subclusters using different shades of red
   - Enhanced titles indicating when clusters are reclustered subclusters
   - UMAP overview plot highlighting original cluster 0 cells in red against light gray background
4. 4.
   Comprehensive Summary Report : Generates a detailed analysis summary that includes:
   
   - Reclustering success/failure per sample
   - Number of subclusters created
   - Cluster size analysis for subclusters only
   - Overall statistics across all samples
   - Note that original clusters are not plotted individually
Output Location : The script saves all plots to the `plots_hippo_recluster_highlight_inclusive` directory and creates a summary report file. Key output files include:

- Individual subcluster plots (UMAP and spatial) for each new subcluster
- `umap_overview_cluster0_highlighted.png` - overview plot showing original cluster 0 location
- `hippo_recluster_analysis_summary.txt` - comprehensive analysis report

The script focuses exclusively on visualizing the newly defined cluster 0 subclusters, providing targeted analysis of the reclustering results while maintaining the same high-quality plotting style as the original hippo_analysis_inclusive.py.