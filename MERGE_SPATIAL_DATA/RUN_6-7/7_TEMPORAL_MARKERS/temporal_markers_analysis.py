# %%
import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Configure scanpy settings
sc.settings.verbosity = 3  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')

# %%
# Define the data structure and metadata (same as in merge script)
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
# Load the merged spatial data
print("Loading merged spatial data...")
merged_data_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
if not os.path.exists(merged_data_path):
    print(f"Error: Merged data file not found at {merged_data_path}")
    print("Please run the merge script first.")
    sys.exit(1)

adata = sc.read_h5ad(merged_data_path)
print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} genes")

# %%
# Create output directory
output_dir = 'MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/results'
os.makedirs(output_dir, exist_ok=True)
print(f"Results will be saved to: {output_dir}")

# Set scanpy figure directory
sc.settings.figdir = output_dir

# %%
# Identify unique time points and their corresponding ctrl/mut pairs
time_points = adata.obs['time_point'].unique()
print(f"Available time points: {time_points}")

# Define comparison pairs for each time point
comparisons = []
for tp in time_points:
    tp_data = adata[adata.obs['time_point'] == tp]
    types_available = tp_data.obs['type'].unique()
    if 'ctrl' in types_available and 'mut' in types_available:
        comparisons.append(tp)
        print(f"Time point {tp}: {len(tp_data[tp_data.obs['type'] == 'ctrl'])} ctrl cells, {len(tp_data[tp_data.obs['type'] == 'mut'])} mut cells")
    else:
        print(f"Skipping time point {tp}: missing ctrl or mut samples")

print(f"\nWill perform differential expression analysis for: {comparisons}")

# %%
# Function to perform differential expression analysis
def perform_de_analysis(adata_subset, time_point, output_dir):
    """
    Perform differential expression analysis between ctrl and mut for a specific time point
    """
    print(f"\n=== Analyzing time point: {time_point} ===")
    
    # Ensure we have both ctrl and mut
    types = adata_subset.obs['type'].unique()
    if len(types) < 2:
        print(f"Skipping {time_point}: insufficient sample types")
        return None
    
    # Perform differential expression using Wilcoxon rank-sum test
    sc.tl.rank_genes_groups(
        adata_subset, 
        groupby='type', 
        groups=['mut'], 
        reference='ctrl',
        method='wilcoxon',
        key_added=f'de_{time_point}'
    )
    
    # Extract results
    de_results = sc.get.rank_genes_groups_df(
        adata_subset, 
        group='mut', 
        key=f'de_{time_point}'
    )
    
    # Add time point information
    de_results['time_point'] = time_point
    
    # Calculate additional statistics
    de_results['abs_logfoldchanges'] = np.abs(de_results['logfoldchanges'])
    de_results['-log10_pvals_adj'] = -np.log10(de_results['pvals_adj'] + 1e-300)  # Add small value to avoid log(0)
    
    # Define significance thresholds
    pval_threshold = 0.05
    logfc_threshold = 0.5
    
    de_results['significant'] = (
        (de_results['pvals_adj'] < pval_threshold) & 
        (de_results['abs_logfoldchanges'] > logfc_threshold)
    )
    
    # Classify genes as upregulated or downregulated
    de_results['regulation'] = 'not_significant'
    de_results.loc[
        (de_results['significant']) & (de_results['logfoldchanges'] > 0), 
        'regulation'
    ] = 'upregulated'
    de_results.loc[
        (de_results['significant']) & (de_results['logfoldchanges'] < 0), 
        'regulation'
    ] = 'downregulated'
    
    print(f"Total genes analyzed: {len(de_results)}")
    print(f"Significantly upregulated: {sum(de_results['regulation'] == 'upregulated')}")
    print(f"Significantly downregulated: {sum(de_results['regulation'] == 'downregulated')}")
    
    # Save results
    results_file = os.path.join(output_dir, f'de_results_{time_point}_mut_vs_ctrl.csv')
    de_results.to_csv(results_file, index=False)
    print(f"Results saved to: {results_file}")
    
    return de_results

# %%
# Function to create volcano plot
def create_volcano_plot(de_results, time_point, output_dir):
    """
    Create volcano plot for differential expression results
    """
    plt.figure(figsize=(10, 8))
    
    # Create scatter plot
    colors = {'upregulated': 'red', 'downregulated': 'blue', 'not_significant': 'gray'}
    
    for reg_type in ['not_significant', 'downregulated', 'upregulated']:
        subset = de_results[de_results['regulation'] == reg_type]
        plt.scatter(
            subset['logfoldchanges'], 
            subset['-log10_pvals_adj'],
            c=colors[reg_type],
            alpha=0.6,
            s=20,
            label=f"{reg_type.replace('_', ' ').title()} ({len(subset)})"
        )
    
    # Add threshold lines
    plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5, label='p-adj = 0.05')
    plt.axvline(x=0.5, color='black', linestyle='--', alpha=0.5)
    plt.axvline(x=-0.5, color='black', linestyle='--', alpha=0.5, label='|log2FC| = 0.5')
    
    # Customize plot
    plt.xlabel('Log2 Fold Change (mut vs ctrl)')
    plt.ylabel('-Log10 Adjusted P-value')
    plt.title(f'Volcano Plot - {time_point.upper()} (mut vs ctrl)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    
    # Save plot
    plot_file = os.path.join(output_dir, f'volcano_plot_{time_point}_mut_vs_ctrl.png')
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Volcano plot saved to: {plot_file}")

# %%
# Function to create heatmap of top genes
def create_heatmap(adata_subset, de_results, time_point, output_dir, n_genes=20):
    """
    Create heatmap of top differentially expressed genes
    """
    # Get top upregulated and downregulated genes
    top_up = de_results[
        de_results['regulation'] == 'upregulated'
    ].nsmallest(n_genes//2, 'pvals_adj')['names'].tolist()
    
    top_down = de_results[
        de_results['regulation'] == 'downregulated'
    ].nsmallest(n_genes//2, 'pvals_adj')['names'].tolist()
    
    top_genes = top_up + top_down
    
    if len(top_genes) == 0:
        print(f"No significant genes found for {time_point}")
        return
    
    # Filter genes that exist in the data
    top_genes = [g for g in top_genes if g in adata_subset.var_names]
    
    if len(top_genes) == 0:
        print(f"No top genes found in data for {time_point}")
        return
    
    print(f"Creating heatmap with {len(top_genes)} genes for {time_point}")
    
    # Create heatmap
    sc.pl.heatmap(
        adata_subset,
        var_names=top_genes,
        groupby='type',
        cmap='RdBu_r',
        dendrogram=False,
        save=f'_heatmap_{time_point}_top_genes.png',
        show=False
    )

# %%
# Perform differential expression analysis for each time point
all_de_results = []

for time_point in comparisons:
    print(f"\n{'='*50}")
    print(f"Processing time point: {time_point}")
    print(f"{'='*50}")
    
    # Subset data for this time point
    adata_subset = adata[adata.obs['time_point'] == time_point].copy()
    
    # Perform DE analysis
    de_results = perform_de_analysis(adata_subset, time_point, output_dir)
    
    if de_results is not None:
        all_de_results.append(de_results)
        
        # Create volcano plot
        create_volcano_plot(de_results, time_point, output_dir)
        
        # Create heatmap
        create_heatmap(adata_subset, de_results, time_point, output_dir)

# %%
# Combine all results and create summary
if all_de_results:
    print("\n=== Creating Summary ===")
    
    # Combine all DE results
    combined_results = pd.concat(all_de_results, ignore_index=True)
    
    # Save combined results
    combined_file = os.path.join(output_dir, 'combined_de_results_all_timepoints.csv')
    combined_results.to_csv(combined_file, index=False)
    print(f"Combined results saved to: {combined_file}")
    
    # Create summary statistics
    summary_stats = []
    for time_point in comparisons:
        tp_results = combined_results[combined_results['time_point'] == time_point]
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
    summary_file = os.path.join(output_dir, 'summary_statistics.csv')
    summary_df.to_csv(summary_file, index=False)
    print(f"Summary statistics saved to: {summary_file}")
    
    # Print summary
    print("\n=== SUMMARY ===")
    print(summary_df.to_string(index=False))
    
    # Create summary bar plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot 1: Number of significant genes per time point
    ax1.bar(summary_df['time_point'], summary_df['significant_genes'], color='skyblue')
    ax1.set_xlabel('Time Point')
    ax1.set_ylabel('Number of Significant Genes')
    ax1.set_title('Significant Genes per Time Point')
    ax1.tick_params(axis='x', rotation=45)
    
    # Plot 2: Up vs Down regulation
    x = np.arange(len(summary_df))
    width = 0.35
    ax2.bar(x - width/2, summary_df['upregulated'], width, label='Upregulated', color='red', alpha=0.7)
    ax2.bar(x + width/2, summary_df['downregulated'], width, label='Downregulated', color='blue', alpha=0.7)
    ax2.set_xlabel('Time Point')
    ax2.set_ylabel('Number of Genes')
    ax2.set_title('Up vs Down Regulated Genes')
    ax2.set_xticks(x)
    ax2.set_xticklabels(summary_df['time_point'])
    ax2.legend()
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    summary_plot_file = os.path.join(output_dir, 'summary_barplot.png')
    plt.savefig(summary_plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Summary plot saved to: {summary_plot_file}")
    
else:
    print("No differential expression results generated.")

print("\n=== Analysis Complete ===")
print(f"All results saved to: {output_dir}")

# %%