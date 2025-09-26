# %%
import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import geopandas as gpd
from shapely import wkt
from shapely.geometry import Point, Polygon
import warnings
warnings.filterwarnings('ignore')

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Configure scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# %%
# =============================================================================
# CONFIGURATION PARAMETERS
# =============================================================================

# Spatial binning parameters
BIN_SIZE = 100  # μm, adjustable (50, 100, 200)
MIN_TRANSCRIPTS_PER_BIN = 10  # Minimum transcripts for statistical analysis
ROI_GROUP = "hippo"  # ROI group to analyze (matching existing pipeline)

# Statistical parameters
PVAL_THRESHOLD = 0.05
LOGFC_THRESHOLD = 0.5
MULTIPLE_TESTING_METHOD = "fdr_bh"

# %%
# =============================================================================
# DATA MAPPING AND METADATA
# =============================================================================

# Define the data structure and metadata (same as existing pipeline)
data_map = {
    'p0-p7': {
        'R1': {'time_point': 'p7', 'region': 'R1', 'type': 'mut'},
        'R2': {'time_point': 'p7', 'region': 'R2', 'type': 'ctrl'},
        'R4': {'time_point': 'p30', 'region': 'R4p', 'type': 'ctrl'},
        'R5': {'time_point': 'p0', 'region': 'R5p', 'type': 'mut'},
        'R6': {'time_point': 'p0', 'region': 'R6', 'type': 'ctrl'}
    },
    'p30-E165': {
        'R1': {'time_point': 'e16', 'region': 'R1', 'type': 'mut'},
        'R2': {'time_point': 'e16', 'region': 'R2', 'type': 'ctrl'},
        'R3': {'time_point': 'p30', 'region': 'R3', 'type': 'mut'}
    }
}

# %%
# =============================================================================
# ROI GEOMETRY LOADING AND PROCESSING FUNCTIONS
# =============================================================================

def load_roi_geometry(sample_id, base_data_dir='DATA'):
    """
    Load ROI geometry for a specific sample

    Parameters:
    -----------
    sample_id : str
        Sample identifier in format {time_point}_{region}_{type}
    base_data_dir : str
        Base directory containing ROI files

    Returns:
    --------
    roi_gdf : GeoDataFrame
        ROI geometries as GeoDataFrame
    """
    roi_path = os.path.join(base_data_dir, "ROIs", f"{sample_id}_geometry.csv")

    if not os.path.exists(roi_path):
        raise FileNotFoundError(f"ROI file not found: {roi_path}")

    print(f"Loading ROI geometry from: {roi_path}")
    roi_df = pd.read_csv(roi_path)

    # Convert WKT string geometries to Shapely geometry objects
    roi_df['geometry'] = roi_df['geometry'].apply(wkt.loads)

    # Create GeoDataFrame (CRS will be set based on data)
    roi_gdf = gpd.GeoDataFrame(roi_df, geometry='geometry')

    return roi_gdf

def filter_roi_by_group(roi_gdf, target_group=ROI_GROUP):
    """
    Filter ROI geometries to specific group (e.g., "hippo")

    Parameters:
    -----------
    roi_gdf : GeoDataFrame
        ROI geometries
    target_group : str
        Target ROI group to filter for

    Returns:
    --------
    filtered_roi_gdf : GeoDataFrame
        Filtered ROI geometries
    """
    if 'group' not in roi_gdf.columns:
        raise ValueError("ROI GeoDataFrame must contain 'group' column")

    filtered_roi_gdf = roi_gdf[roi_gdf['group'] == target_group].copy()

    if len(filtered_roi_gdf) == 0:
        raise ValueError(f"No ROI polygons found for group '{target_group}'")

    print(f"Filtered to {len(filtered_roi_gdf)} ROI polygons for group '{target_group}'")
    return filtered_roi_gdf

# %%
# =============================================================================
# SPATIAL BINNING FUNCTIONS
# =============================================================================

def create_spatial_bins(roi_gdf, bin_size=BIN_SIZE):
    """
    Create regular spatial bins within ROI boundaries

    Parameters:
    -----------
    roi_gdf : GeoDataFrame
        ROI geometries to create bins within
    bin_size : float
        Size of each bin in micrometers

    Returns:
    --------
    bins_gdf : GeoDataFrame
        Spatial bins as square polygons
    """
    # Get the bounds of all ROI polygons
    total_bounds = roi_gdf.total_bounds  # [minx, miny, maxx, maxy]

    # Create grid of bin centers
    x_bins = np.arange(total_bounds[0], total_bounds[2] + bin_size, bin_size)
    y_bins = np.arange(total_bounds[1], total_bounds[3] + bin_size, bin_size)

    # Create bin polygons
    bin_polygons = []
    bin_ids = []
    bin_centers = []

    for i, x in enumerate(x_bins[:-1]):
        for j, y in enumerate(y_bins[:-1]):
            # Create square bin polygon
            bin_polygon = Polygon([
                (x, y),
                (x + bin_size, y),
                (x + bin_size, y + bin_size),
                (x, y + bin_size)
            ])

            # Check if bin intersects with any ROI polygon
            if roi_gdf.intersects(bin_polygon).any():
                bin_polygons.append(bin_polygon)
                bin_ids.append(f"bin_{i}_{j}")
                bin_centers.append((x + bin_size/2, y + bin_size/2))

    # Create GeoDataFrame for bins
    bins_gdf = gpd.GeoDataFrame({
        'bin_id': bin_ids,
        'center_x': [center[0] for center in bin_centers],
        'center_y': [center[1] for center in bin_centers],
        'geometry': bin_polygons
    }, crs=roi_gdf.crs)

    print(f"Created {len(bins_gdf)} spatial bins of size {bin_size}μm")
    return bins_gdf

def aggregate_transcripts_to_bins(adata, bins_gdf, min_transcripts=MIN_TRANSCRIPTS_PER_BIN):
    """
    Aggregate transcript counts to spatial bins

    Parameters:
    -----------
    adata : AnnData
        Original spatial transcriptomics data
    bins_gdf : GeoDataFrame
        Spatial bins for aggregation
    min_transcripts : int
        Minimum transcripts per bin for inclusion

    Returns:
    --------
    bin_expression_df : DataFrame
        Expression matrix with bins as rows and genes as columns
    bin_metadata_df : DataFrame
        Metadata for each bin
    """
    # Get spatial coordinates from obs columns
    if 'center_x' not in adata.obs.columns or 'center_y' not in adata.obs.columns:
        raise ValueError("AnnData object must contain spatial coordinates in obs['center_x'] and obs['center_y']")

    # Extract spatial coordinates
    spatial_coords = adata.obs[['center_x', 'center_y']].values

    # Create GeoDataFrame for cell/transcript locations
    cell_points = [Point(x, y) for x, y in spatial_coords]
    cells_gdf = gpd.GeoDataFrame({
        'cell_id': adata.obs.index,
        'geometry': cell_points
    }, crs=bins_gdf.crs)

    # Spatial join to assign cells to bins
    cells_with_bins = gpd.sjoin(cells_gdf, bins_gdf, how='inner', predicate='within')

    # Aggregate expression data by bin
    bin_expression_data = []
    bin_metadata = []

    for bin_id in bins_gdf['bin_id']:
        # Get cells in this bin
        bin_cells = cells_with_bins[cells_with_bins['bin_id'] == bin_id]

        if len(bin_cells) == 0:
            continue

        # Get cell indices for this bin
        cell_indices = [adata.obs.index.get_loc(cell_id) for cell_id in bin_cells['cell_id']]

        # Sum expression across cells in bin
        bin_expression = adata.X[cell_indices, :].sum(axis=0)

        # Convert to array if sparse
        if hasattr(bin_expression, 'A1'):
            bin_expression = bin_expression.A1

        # Check if bin meets minimum transcript threshold
        total_transcripts = bin_expression.sum()
        if total_transcripts >= min_transcripts:
            bin_expression_data.append(bin_expression)

            # Get bin center coordinates
            bin_info = bins_gdf[bins_gdf['bin_id'] == bin_id].iloc[0]
            bin_metadata.append({
                'bin_id': bin_id,
                'center_x': bin_info['center_x'],
                'center_y': bin_info['center_y'],
                'n_cells': len(bin_cells),
                'total_transcripts': total_transcripts
            })

    # Create expression DataFrame
    if len(bin_expression_data) == 0:
        raise ValueError("No bins met the minimum transcript threshold")

    bin_expression_df = pd.DataFrame(
        bin_expression_data,
        columns=adata.var.index,
        index=[meta['bin_id'] for meta in bin_metadata]
    )

    bin_metadata_df = pd.DataFrame(bin_metadata)
    bin_metadata_df.set_index('bin_id', inplace=True)

    print(f"Created expression matrix for {len(bin_expression_df)} bins")
    print(f"Bins contain {bin_metadata_df['n_cells'].sum()} total cells")
    print(f"Average transcripts per bin: {bin_metadata_df['total_transcripts'].mean():.1f}")

    return bin_expression_df, bin_metadata_df

# %%
# =============================================================================
# STATISTICAL ANALYSIS FUNCTIONS
# =============================================================================

def perform_roi_de_analysis(ctrl_expression, mut_expression, time_point, output_dir):
    """
    Perform differential expression analysis between ctrl and mut bins

    Parameters:
    -----------
    ctrl_expression : DataFrame
        Expression data for control bins
    mut_expression : DataFrame
        Expression data for mutant bins
    time_point : str
        Time point being analyzed
    output_dir : str
        Output directory for results

    Returns:
    --------
    de_results : DataFrame
        Differential expression results
    """
    print(f"\n=== Analyzing time point: {time_point} ===")
    print(f"Control bins: {len(ctrl_expression)}")
    print(f"Mutant bins: {len(mut_expression)}")

    if len(ctrl_expression) == 0 or len(mut_expression) == 0:
        print(f"Insufficient bins for analysis at {time_point}")
        return None

    # Initialize results list
    results = []

    # Perform statistical test for each gene
    for gene in ctrl_expression.columns:
        ctrl_values = ctrl_expression[gene].values
        mut_values = mut_expression[gene].values

        # Calculate basic statistics
        ctrl_mean = np.mean(ctrl_values)
        mut_mean = np.mean(mut_values)

        # Skip genes with zero expression in both groups
        if ctrl_mean == 0 and mut_mean == 0:
            continue

        # Calculate log2 fold change (add pseudocount to avoid log(0))
        pseudocount = 1
        log2fc = np.log2((mut_mean + pseudocount) / (ctrl_mean + pseudocount))

        # Perform Mann-Whitney U test
        try:
            statistic, p_value = stats.mannwhitneyu(
                mut_values, ctrl_values,
                alternative='two-sided'
            )
        except ValueError:
            # Handle cases where one group has constant values
            p_value = 1.0
            statistic = np.nan

        results.append({
            'gene': gene,
            'ctrl_mean': ctrl_mean,
            'mut_mean': mut_mean,
            'log2fc': log2fc,
            'pvalue': p_value,
            'statistic': statistic,
            'time_point': time_point
        })

    # Convert to DataFrame
    de_results = pd.DataFrame(results)

    if len(de_results) == 0:
        print("No genes passed filtering criteria")
        return None

    # Multiple testing correction
    from statsmodels.stats.multitest import multipletests

    _, p_adj, _, _ = multipletests(
        de_results['pvalue'],
        method=MULTIPLE_TESTING_METHOD
    )
    de_results['pvalue_adj'] = p_adj

    # Calculate additional statistics
    de_results['abs_log2fc'] = np.abs(de_results['log2fc'])
    de_results['-log10_pvalue'] = -np.log10(de_results['pvalue'] + 1e-300)
    de_results['-log10_pvalue_adj'] = -np.log10(de_results['pvalue_adj'] + 1e-300)

    # Define significance
    de_results['significant'] = (
        (de_results['pvalue_adj'] < PVAL_THRESHOLD) &
        (de_results['abs_log2fc'] > LOGFC_THRESHOLD)
    )

    # Classify regulation direction
    de_results['regulation'] = 'not_significant'
    de_results.loc[
        (de_results['significant']) & (de_results['log2fc'] > 0),
        'regulation'
    ] = 'upregulated'
    de_results.loc[
        (de_results['significant']) & (de_results['log2fc'] < 0),
        'regulation'
    ] = 'downregulated'

    # Sort by significance
    de_results = de_results.sort_values('pvalue_adj')

    print(f"Total genes analyzed: {len(de_results)}")
    print(f"Significantly upregulated: {sum(de_results['regulation'] == 'upregulated')}")
    print(f"Significantly downregulated: {sum(de_results['regulation'] == 'downregulated')}")

    # Save results
    results_file = os.path.join(output_dir, f'roi_de_results_{time_point}_mut_vs_ctrl.csv')
    de_results.to_csv(results_file, index=False)
    print(f"Results saved to: {results_file}")

    return de_results

# %%
# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def create_roi_volcano_plot(de_results, time_point, output_dir):
    """
    Create volcano plot for ROI-based differential expression results
    """
    if de_results is None or len(de_results) == 0:
        return

    plt.figure(figsize=(10, 8))

    # Create scatter plot
    colors = {'upregulated': 'red', 'downregulated': 'blue', 'not_significant': 'gray'}

    for reg_type in ['not_significant', 'downregulated', 'upregulated']:
        subset = de_results[de_results['regulation'] == reg_type]
        if len(subset) > 0:
            plt.scatter(
                subset['log2fc'],
                subset['-log10_pvalue_adj'],
                c=colors[reg_type],
                alpha=0.6,
                s=20,
                label=f"{reg_type.replace('_', ' ').title()} ({len(subset)})"
            )

    # Add threshold lines
    plt.axhline(y=-np.log10(PVAL_THRESHOLD), color='black', linestyle='--', alpha=0.5,
                label=f'p-adj = {PVAL_THRESHOLD}')
    plt.axvline(x=LOGFC_THRESHOLD, color='black', linestyle='--', alpha=0.5)
    plt.axvline(x=-LOGFC_THRESHOLD, color='black', linestyle='--', alpha=0.5,
                label=f'|log2FC| = {LOGFC_THRESHOLD}')

    # Customize plot
    plt.xlabel('Log2 Fold Change (mut vs ctrl)')
    plt.ylabel('-Log10 Adjusted P-value')
    plt.title(f'ROI-Based Volcano Plot - {time_point.upper()} (mut vs ctrl)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)

    # Save plot
    plot_file = os.path.join(output_dir, f'roi_volcano_plot_{time_point}_mut_vs_ctrl.png')
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Volcano plot saved to: {plot_file}")

def create_roi_spatial_heatmap(bin_metadata_df, de_results, time_point, output_dir, top_n=10):
    """
    Create spatial heatmap showing top DE genes across ROI bins
    """
    if de_results is None or len(de_results) == 0:
        return

    # Get top upregulated and downregulated genes
    top_up = de_results[
        de_results['regulation'] == 'upregulated'
    ].head(top_n//2)

    top_down = de_results[
        de_results['regulation'] == 'downregulated'
    ].head(top_n//2)

    top_genes = list(top_up['gene']) + list(top_down['gene'])

    if len(top_genes) == 0:
        print(f"No significant genes found for spatial heatmap at {time_point}")
        return

    # Create summary spatial plot
    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot bin centers colored by significance
    significant_bins = bin_metadata_df.index

    scatter = ax.scatter(
        bin_metadata_df.loc[significant_bins, 'center_x'],
        bin_metadata_df.loc[significant_bins, 'center_y'],
        c=bin_metadata_df.loc[significant_bins, 'total_transcripts'],
        s=50,
        alpha=0.7,
        cmap='viridis'
    )

    plt.colorbar(scatter, label='Total Transcripts per Bin')
    ax.set_xlabel('X Coordinate (μm)')
    ax.set_ylabel('Y Coordinate (μm)')
    ax.set_title(f'ROI Spatial Bins - {time_point.upper()} (n={len(significant_bins)} bins)')

    # Save plot
    plot_file = os.path.join(output_dir, f'roi_spatial_bins_{time_point}.png')
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Spatial bins plot saved to: {plot_file}")

# %%
# =============================================================================
# MAIN ANALYSIS PIPELINE
# =============================================================================

def main():
    """
    Main function to run ROI-based temporal markers analysis
    """
    print("Starting ROI-Based Temporal Markers Analysis")
    print("=" * 60)

    # Create output directory
    output_dir = 'MERGE_SPATIAL_DATA/RUN_6-7/7_TEMPORAL_MARKERS/roi_results'
    os.makedirs(output_dir, exist_ok=True)
    print(f"Results will be saved to: {output_dir}")

    # Process each sample and collect bin expression data
    sample_bin_data = {}
    base_data_dir = 'DATA'

    for time_group, regions in data_map.items():
        for region, metadata in regions.items():
            sample_id = f"{metadata['time_point']}_{metadata['region']}_{metadata['type']}"
            print(f"\nProcessing sample: {sample_id}")

            try:
                # Load original spatial data
                h5ad_path = os.path.join(base_data_dir, time_group, region, 'data.h5ad')
                if not os.path.exists(h5ad_path):
                    print(f"Data file not found: {h5ad_path}")
                    continue

                print(f"Loading spatial data from: {h5ad_path}")
                adata = sc.read_h5ad(h5ad_path)

                # Add metadata
                adata.obs['time_point'] = metadata['time_point']
                adata.obs['region'] = metadata['region']
                adata.obs['type'] = metadata['type']
                adata.obs['sample_id'] = sample_id

                # Load ROI geometry
                roi_gdf = load_roi_geometry(sample_id, base_data_dir)

                # Filter to target ROI group
                roi_gdf_filtered = filter_roi_by_group(roi_gdf, ROI_GROUP)

                # Create spatial bins
                bins_gdf = create_spatial_bins(roi_gdf_filtered, BIN_SIZE)

                # Aggregate transcripts to bins
                bin_expression_df, bin_metadata_df = aggregate_transcripts_to_bins(
                    adata, bins_gdf, MIN_TRANSCRIPTS_PER_BIN
                )

                # Store data for this sample
                sample_bin_data[sample_id] = {
                    'expression': bin_expression_df,
                    'metadata': bin_metadata_df,
                    'sample_metadata': metadata
                }

                # Save bin expression matrix
                expression_file = os.path.join(output_dir, f'roi_bin_expression_{sample_id}.csv')
                bin_expression_df.to_csv(expression_file)

                metadata_file = os.path.join(output_dir, f'roi_bin_metadata_{sample_id}.csv')
                bin_metadata_df.to_csv(metadata_file)

                print(f"Saved bin data for {sample_id}")

            except Exception as e:
                print(f"Error processing {sample_id}: {str(e)}")
                continue

    # Perform differential expression analysis for each time point
    print("\n" + "=" * 60)
    print("PERFORMING DIFFERENTIAL EXPRESSION ANALYSIS")
    print("=" * 60)

    # Identify unique time points and their corresponding ctrl/mut pairs
    time_points = set()
    for sample_data in sample_bin_data.values():
        time_points.add(sample_data['sample_metadata']['time_point'])

    all_de_results = []

    for time_point in time_points:
        print(f"\nAnalyzing time point: {time_point}")

        # Get ctrl and mut samples for this time point
        ctrl_samples = []
        mut_samples = []

        for sample_id, data in sample_bin_data.items():
            if data['sample_metadata']['time_point'] == time_point:
                if data['sample_metadata']['type'] == 'ctrl':
                    ctrl_samples.append(data['expression'])
                elif data['sample_metadata']['type'] == 'mut':
                    mut_samples.append(data['expression'])

        if len(ctrl_samples) == 0 or len(mut_samples) == 0:
            print(f"Skipping {time_point}: missing ctrl or mut samples")
            continue

        # Combine expression data for ctrl and mut
        ctrl_expression = pd.concat(ctrl_samples, ignore_index=True)
        mut_expression = pd.concat(mut_samples, ignore_index=True)

        # Perform DE analysis
        de_results = perform_roi_de_analysis(
            ctrl_expression, mut_expression, time_point, output_dir
        )

        if de_results is not None:
            all_de_results.append(de_results)

            # Create visualizations
            create_roi_volcano_plot(de_results, time_point, output_dir)

            # Create spatial heatmap using metadata from first sample
            sample_metadata = None
            for sample_id, data in sample_bin_data.items():
                if data['sample_metadata']['time_point'] == time_point:
                    sample_metadata = data['metadata']
                    break

            if sample_metadata is not None:
                create_roi_spatial_heatmap(
                    sample_metadata, de_results, time_point, output_dir
                )

    # Combine all results and create summary
    if all_de_results:
        print("\n" + "=" * 60)
        print("CREATING SUMMARY")
        print("=" * 60)

        # Combine all DE results
        combined_results = pd.concat(all_de_results, ignore_index=True)

        # Save combined results
        combined_file = os.path.join(output_dir, 'roi_combined_results_all_timepoints.csv')
        combined_results.to_csv(combined_file, index=False)
        print(f"Combined results saved to: {combined_file}")

        # Create summary statistics
        summary_stats = []
        for time_point in time_points:
            tp_results = combined_results[combined_results['time_point'] == time_point]
            if len(tp_results) > 0:
                stats_dict = {
                    'time_point': time_point,
                    'total_genes': len(tp_results),
                    'significant_genes': sum(tp_results['significant']),
                    'upregulated': sum(tp_results['regulation'] == 'upregulated'),
                    'downregulated': sum(tp_results['regulation'] == 'downregulated'),
                    'pct_significant': (sum(tp_results['significant']) / len(tp_results)) * 100
                }
                summary_stats.append(stats_dict)

        summary_df = pd.DataFrame(summary_stats)
        summary_file = os.path.join(output_dir, 'roi_summary_statistics.csv')
        summary_df.to_csv(summary_file, index=False)
        print(f"Summary statistics saved to: {summary_file}")

        # Print summary
        print("\n=== ROI-BASED ANALYSIS SUMMARY ===")
        print(summary_df.to_string(index=False))

        # Create summary analysis report
        report_file = os.path.join(output_dir, 'roi_analysis_summary.txt')
        with open(report_file, 'w') as f:
            f.write("ROI-Based Temporal Markers Analysis Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Analysis Parameters:\n")
            f.write(f"- Bin size: {BIN_SIZE} μm\n")
            f.write(f"- Minimum transcripts per bin: {MIN_TRANSCRIPTS_PER_BIN}\n")
            f.write(f"- ROI group analyzed: {ROI_GROUP}\n")
            f.write(f"- P-value threshold: {PVAL_THRESHOLD}\n")
            f.write(f"- Log2FC threshold: {LOGFC_THRESHOLD}\n\n")
            f.write("Results Summary:\n")
            f.write(summary_df.to_string(index=False))
            f.write(f"\n\nTotal samples processed: {len(sample_bin_data)}\n")
            f.write(f"Output directory: {output_dir}\n")

        print(f"Analysis report saved to: {report_file}")

    else:
        print("No differential expression results generated.")

    print("\n=== ROI-BASED ANALYSIS COMPLETE ===")
    print(f"All results saved to: {output_dir}")

# %%
# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

if __name__ == "__main__":
    main()

# %%