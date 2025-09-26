# UMAP Parameter Optimization

## Parameters Tested
- **Min_dist** `[0.01, 0.1, 0.3, 0.5, 0.8]`
  - Minimum distance between embedded points
- **Spread** `[0.5, 1.0, 1.5, 2.0]`
  - Controls separation between points in final plot
- **Gamma** `[1.0, 2.0, 5.0, 10.0]`
  - Controls balance between attractive and repulsive forces during optimization

## Results Location
- Analysis summaries: `MERGE_SPATIAL_DATA/umap_optimization/summary/`
- Individual UMAP plots: `MERGE_SPATIAL_DATA/umap_optimization/plots/{param_combination}/`
- Combined UMAP plots: `MERGE_SPATIAL_DATA/umap_parameter_collages/`

## Key Metrics
- **Total_spread**: Measures embedding dispersion
  - Higher values indicate greater point separation
  - Helps with cluster distinction
- **Neighborhood_preservation**: Measures maintenance of local relationships
  - Scale: 0.0 to 1.0 (higher is better)
  - Indicates how well original high-dimensional relationships are preserved in 2D
