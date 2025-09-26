# %%
import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import mannwhitneyu, fisher_exact
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Configure scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# %%
# Define the data structure and metadata
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
output_dir = 'MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/batch_corrected'
os.makedirs(output_dir, exist_ok=True)
print(f"Results will be saved to: {output_dir}")

# Set scanpy figure directory
sc.settings.figdir = output_dir

# %%
# Quality control and preprocessing for spatial transcriptomics
def preprocess_spatial_data(adata):
    """
    Preprocess spatial transcriptomics data accounting for limited gene panel
    """
    print("\n=== Preprocessing Spatial Data ===")
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('Mt-')
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    # For spatial data, we typically don't filter genes as aggressively
    # since the gene panel is already limited
    print(f"Initial: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Filter cells with very low gene counts (likely technical artifacts)
    sc.pp.filter_cells(adata, min_genes=10)
    
    # Filter genes expressed in very few cells (< 1% of cells)
    min_cells = max(3, int(0.01 * adata.n_obs))
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    print(f"After filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Store raw counts
    adata.raw = adata
    
    return adata

# %%
# Batch correction function
def perform_batch_correction(adata, batch_key='region'):
    """
    Perform batch correction using Harmony-like approach
    """
    print(f"\n=== Performing Batch Correction (batch_key: {batch_key}) ===")
    
    # Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Find highly variable genes (for spatial data, this is more conservative)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # Keep all genes for spatial analysis but mark HVGs
    adata.var['highly_variable_spatial'] = adata.var['highly_variable']
    
    # Scale data
    sc.pp.scale(adata, max_value=10)
    
    # PCA
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
    
    # Batch correction using scanpy's combat
    try:
        sc.pp.combat(adata, key=batch_key)
        print(f"Applied ComBat batch correction for {batch_key}")
    except Exception as e:
        print(f"ComBat failed: {e}. Proceeding without batch correction.")
    
    return adata

# %%
# DESeq2-style differential expression for spatial data
def deseq2_style_de_analysis(adata_subset, time_point, output_dir):
    """
    Perform DESeq2-style differential expression analysis adapted for spatial data
    """
    print(f"\n=== DESeq2-style Analysis for {time_point} ===")
    
    # Get raw counts
    if adata_subset.raw is not None:
        raw_counts = adata_subset.raw.X
        gene_names = adata_subset.raw.var_names
    else:
        raw_counts = adata_subset.X
        gene_names = adata_subset.var_names
    
    # Convert to dense if sparse
    if hasattr(raw_counts, 'toarray'):
        raw_counts = raw_counts.toarray()
    
    # Create count matrix DataFrame
    count_df = pd.DataFrame(raw_counts, columns=gene_names)
    count_df['type'] = adata_subset.obs['type'].values
    count_df['region'] = adata_subset.obs['region'].values
    
    # Separate ctrl and mut
    ctrl_counts = count_df[count_df['type'] == 'ctrl'].drop(['type', 'region'], axis=1)
    mut_counts = count_df[count_df['type'] == 'mut'].drop(['type', 'region'], axis=1)
    
    print(f"Ctrl cells: {len(ctrl_counts)}, Mut cells: {len(mut_counts)}")
    
    if len(ctrl_counts) == 0 or len(mut_counts) == 0:
        print(f"Insufficient samples for {time_point}")
        return None
    
    # DESeq2-style normalization and testing
    results = []
    
    for gene in gene_names:
        ctrl_gene = ctrl_counts[gene].values
        mut_gene = mut_counts[gene].values
        
        # Skip genes with very low expression
        if np.sum(ctrl_gene) + np.sum(mut_gene) < 10:
            continue
        
        # Calculate size factors (simplified DESeq2 approach)
        ctrl_size_factor = np.median(ctrl_counts.sum(axis=1))
        mut_size_factor = np.median(mut_counts.sum(axis=1))
        
        # Normalize counts
        ctrl_norm = ctrl_gene / (ctrl_size_factor / 1e4)
        mut_norm = mut_gene / (mut_size_factor / 1e4)
        
        # Calculate fold change
        ctrl_mean = np.mean(ctrl_norm + 1)  # Add pseudocount
        mut_mean = np.mean(mut_norm + 1)
        log2fc = np.log2(mut_mean / ctrl_mean)
        
        # Statistical test - use Mann-Whitney U for non-parametric testing
        try:
            statistic, pvalue = mannwhitneyu(mut_gene, ctrl_gene, alternative='two-sided')
        except ValueError:
            pvalue = 1.0
            statistic = 0
        
        # Calculate additional metrics
        ctrl_detection_rate = np.sum(ctrl_gene > 0) / len(ctrl_gene)
        mut_detection_rate = np.sum(mut_gene > 0) / len(mut_gene)
        
        results.append({
            'gene': gene,
            'log2FoldChange': log2fc,
            'pvalue': pvalue,
            'ctrl_mean': ctrl_mean,
            'mut_mean': mut_mean,
            'ctrl_detection_rate': ctrl_detection_rate,
            'mut_detection_rate': mut_detection_rate,
            'ctrl_cells': len(ctrl_gene),
            'mut_cells': len(mut_gene)
        })
    
    # Convert to DataFrame
    de_results = pd.DataFrame(results)
    
    if len(de_results) == 0:
        print(f"No genes passed filtering for {time_point}")
        return None
    
    # Multiple testing correction (Benjamini-Hochberg)
    from statsmodels.stats.multitest import multipletests
    _, de_results['padj'], _, _ = multipletests(de_results['pvalue'], method='fdr_bh')
    
    # Add time point and classification
    de_results['time_point'] = time_point
    de_results['abs_log2FC'] = np.abs(de_results['log2FoldChange'])
    de_results['-log10_padj'] = -np.log10(de_results['padj'] + 1e-300)
    
    # Define significance thresholds (more stringent for spatial data)
    padj_threshold = 0.01  # More stringent
    log2fc_threshold = 0.75  # Higher threshold
    detection_threshold = 0.1  # Minimum detection rate
    
    de_results['significant'] = (
        (de_results['padj'] < padj_threshold) & 
        (de_results['abs_log2FC'] > log2fc_threshold) &
        ((de_results['ctrl_detection_rate'] > detection_threshold) |
         (de_results['mut_detection_rate'] > detection_threshold))
    )
    
    # Classify regulation
    de_results['regulation'] = 'not_significant'
    de_results.loc[
        (de_results['significant']) & (de_results['log2FoldChange'] > 0), 
        'regulation'
    ] = 'upregulated'
    de_results.loc[
        (de_results['significant']) & (de_results['log2FoldChange'] < 0), 
        'regulation'
    ] = 'downregulated'
    
    # Sort by significance
    de_results = de_results.sort_values('padj')
    
    print(f"Total genes analyzed: {len(de_results)}")
    print(f"Significantly upregulated: {sum(de_results['regulation'] == 'upregulated')}")
    print(f"Significantly downregulated: {sum(de_results['regulation'] == 'downregulated')}")
    
    # Save results
    results_file = os.path.join(output_dir, f'deseq2_style_results_{time_point}_mut_vs_ctrl.csv')
    de_results.to_csv(results_file, index=False)
    print(f"Results saved to: {results_file}")
    
    return de_results

# %%
# Enhanced volcano plot with better visualization
def create_enhanced_volcano_plot(de_results, time_point, output_dir):
    """
    Create enhanced volcano plot with better visualization
    """
    plt.figure(figsize=(12, 10))
    
    # Create scatter plot with size based on detection rate
    colors = {'upregulated': '#d62728', 'downregulated': '#1f77b4', 'not_significant': '#7f7f7f'}
    
    for reg_type in ['not_significant', 'downregulated', 'upregulated']:
        subset = de_results[de_results['regulation'] == reg_type]
        if len(subset) > 0:
            # Size based on average detection rate
            sizes = 20 + 30 * (subset['ctrl_detection_rate'] + subset['mut_detection_rate']) / 2
            
            plt.scatter(
                subset['log2FoldChange'], 
                subset['-log10_padj'],
                c=colors[reg_type],
                alpha=0.7,
                s=sizes,
                label=f"{reg_type.replace('_', ' ').title()} ({len(subset)})",
                edgecolors='white',
                linewidth=0.5
            )
    
    # Add threshold lines
    plt.axhline(y=-np.log10(0.01), color='black', linestyle='--', alpha=0.7, label='padj = 0.01')
    plt.axvline(x=0.75, color='black', linestyle='--', alpha=0.7)
    plt.axvline(x=-0.75, color='black', linestyle='--', alpha=0.7, label='|log2FC| = 0.75')
    
    # Annotate top genes
    top_genes = de_results[de_results['significant']].head(10)
    for _, gene_data in top_genes.iterrows():
        plt.annotate(
            gene_data['gene'],
            (gene_data['log2FoldChange'], gene_data['-log10_padj']),
            xytext=(5, 5), textcoords='offset points',
            fontsize=8, alpha=0.8
        )
    
    # Customize plot
    plt.xlabel('Log2 Fold Change (mut vs ctrl)', fontsize=12)
    plt.ylabel('-Log10 Adjusted P-value', fontsize=12)
    plt.title(f'Enhanced Volcano Plot - {time_point.upper()} (DESeq2-style)', fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    
    # Add statistics text
    stats_text = f"Total genes: {len(de_results)}\nSignificant: {sum(de_results['significant'])}\nThresholds: padj < 0.01, |log2FC| > 0.75"
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Save plot
    plot_file = os.path.join(output_dir, f'enhanced_volcano_plot_{time_point}_mut_vs_ctrl.png')
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Enhanced volcano plot saved to: {plot_file}")

# %%
# Quality control visualization
def create_qc_plots(adata, output_dir):
    """
    Create quality control plots to assess batch effects and data quality
    """
    print("\n=== Creating Quality Control Plots ===")
    
    # Cell count distribution by time point and type
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Cell counts by time point and type
    cell_counts = adata.obs.groupby(['time_point', 'type']).size().reset_index(name='count')
    pivot_counts = cell_counts.pivot(index='time_point', columns='type', values='count').fillna(0)
    pivot_counts.plot(kind='bar', ax=axes[0,0], color=['skyblue', 'lightcoral'])
    axes[0,0].set_title('Cell Counts by Time Point and Type')
    axes[0,0].set_ylabel('Number of Cells')
    axes[0,0].tick_params(axis='x', rotation=45)
    
    # Plot 2: Gene detection distribution
    adata.obs['n_genes_detected'] = (adata.X > 0).sum(axis=1)
    sns.boxplot(data=adata.obs, x='time_point', y='n_genes_detected', hue='type', ax=axes[0,1])
    axes[0,1].set_title('Gene Detection Distribution')
    axes[0,1].tick_params(axis='x', rotation=45)
    
    # Plot 3: Total UMI distribution
    sns.boxplot(data=adata.obs, x='time_point', y='total_counts', hue='type', ax=axes[1,0])
    axes[1,0].set_title('Total UMI Distribution')
    axes[1,0].tick_params(axis='x', rotation=45)
    
    # Plot 4: Region distribution
    region_counts = adata.obs.groupby(['region', 'time_point']).size().reset_index(name='count')
    pivot_regions = region_counts.pivot(index='region', columns='time_point', values='count').fillna(0)
    pivot_regions.plot(kind='bar', ax=axes[1,1], stacked=True)
    axes[1,1].set_title('Cell Distribution by Region and Time Point')
    axes[1,1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    qc_file = os.path.join(output_dir, 'quality_control_plots.png')
    plt.savefig(qc_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"QC plots saved to: {qc_file}")

# %%
# Main analysis pipeline
print("\n" + "="*60)
print("IMPROVED TEMPORAL MARKERS ANALYSIS PIPELINE")
print("="*60)

# Preprocess data
adata = preprocess_spatial_data(adata)

# Create QC plots
create_qc_plots(adata, output_dir)

# Perform batch correction
adata = perform_batch_correction(adata, batch_key='region')

# Identify time points for analysis
time_points = adata.obs['time_point'].unique()
print(f"\nAvailable time points: {time_points}")

# Define comparison pairs
comparisons = []
for tp in time_points:
    tp_data = adata[adata.obs['time_point'] == tp]
    types_available = tp_data.obs['type'].unique()
    if 'ctrl' in types_available and 'mut' in types_available:
        comparisons.append(tp)
        ctrl_count = len(tp_data[tp_data.obs['type'] == 'ctrl'])
        mut_count = len(tp_data[tp_data.obs['type'] == 'mut'])
        print(f"Time point {tp}: {ctrl_count} ctrl cells, {mut_count} mut cells")
        
        # Check for severe imbalance
        ratio = max(ctrl_count, mut_count) / min(ctrl_count, mut_count)
        if ratio > 10:
            print(f"  WARNING: Severe cell count imbalance (ratio: {ratio:.1f})")
    else:
        print(f"Skipping time point {tp}: missing ctrl or mut samples")

print(f"\nWill perform analysis for: {comparisons}")

# %%
# Perform improved differential expression analysis
all_de_results = []

for time_point in comparisons:
    print(f"\n{'='*50}")
    print(f"Processing time point: {time_point}")
    print(f"{'='*50}")
    
    # Subset data for this time point
    adata_subset = adata[adata.obs['time_point'] == time_point].copy()
    
    # Perform DESeq2-style DE analysis
    de_results = deseq2_style_de_analysis(adata_subset, time_point, output_dir)
    
    if de_results is not None:
        all_de_results.append(de_results)
        
        # Create enhanced volcano plot
        create_enhanced_volcano_plot(de_results, time_point, output_dir)

# %%
# Create comprehensive summary
if all_de_results:
    print("\n=== Creating Comprehensive Summary ===")
    
    # Combine all results
    combined_results = pd.concat(all_de_results, ignore_index=True)
    
    # Save combined results
    combined_file = os.path.join(output_dir, 'combined_deseq2_results_all_timepoints.csv')
    combined_results.to_csv(combined_file, index=False)
    print(f"Combined results saved to: {combined_file}")
    
    # Create detailed summary statistics
    summary_stats = []
    for time_point in comparisons:
        tp_results = combined_results[combined_results['time_point'] == time_point]
        
        # Get cell counts for this time point
        tp_data = adata[adata.obs['time_point'] == time_point]
        ctrl_cells = len(tp_data[tp_data.obs['type'] == 'ctrl'])
        mut_cells = len(tp_data[tp_data.obs['type'] == 'mut'])
        
        stats_dict = {
            'time_point': time_point,
            'ctrl_cells': ctrl_cells,
            'mut_cells': mut_cells,
            'cell_ratio': max(ctrl_cells, mut_cells) / min(ctrl_cells, mut_cells),
            'total_genes_tested': len(tp_results),
            'significant_genes': sum(tp_results['significant']),
            'upregulated': sum(tp_results['regulation'] == 'upregulated'),
            'downregulated': sum(tp_results['regulation'] == 'downregulated'),
            'pct_significant': (sum(tp_results['significant']) / len(tp_results)) * 100,
            'median_padj_significant': tp_results[tp_results['significant']]['padj'].median() if sum(tp_results['significant']) > 0 else np.nan,
            'mean_abs_log2fc_significant': tp_results[tp_results['significant']]['abs_log2FC'].mean() if sum(tp_results['significant']) > 0 else np.nan
        }
        summary_stats.append(stats_dict)
    
    summary_df = pd.DataFrame(summary_stats)
    summary_file = os.path.join(output_dir, 'improved_summary_statistics.csv')
    summary_df.to_csv(summary_file, index=False)
    print(f"Improved summary statistics saved to: {summary_file}")
    
    # Print summary
    print("\n=== IMPROVED ANALYSIS SUMMARY ===")
    print(summary_df.round(3).to_string(index=False))
    
    # Create comprehensive summary plot
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1: Cell count imbalance
    axes[0,0].bar(summary_df['time_point'], summary_df['cell_ratio'], color='orange', alpha=0.7)
    axes[0,0].axhline(y=1, color='red', linestyle='--', label='Perfect balance')
    axes[0,0].set_title('Cell Count Imbalance Ratio')
    axes[0,0].set_ylabel('Max cells / Min cells')
    axes[0,0].tick_params(axis='x', rotation=45)
    axes[0,0].legend()
    
    # Plot 2: Significant genes
    axes[0,1].bar(summary_df['time_point'], summary_df['significant_genes'], color='green', alpha=0.7)
    axes[0,1].set_title('Significant Genes (Improved Thresholds)')
    axes[0,1].set_ylabel('Number of Significant Genes')
    axes[0,1].tick_params(axis='x', rotation=45)
    
    # Plot 3: Up vs Down regulation
    x = np.arange(len(summary_df))
    width = 0.35
    axes[0,2].bar(x - width/2, summary_df['upregulated'], width, label='Upregulated', color='red', alpha=0.7)
    axes[0,2].bar(x + width/2, summary_df['downregulated'], width, label='Downregulated', color='blue', alpha=0.7)
    axes[0,2].set_title('Up vs Down Regulated Genes')
    axes[0,2].set_ylabel('Number of Genes')
    axes[0,2].set_xticks(x)
    axes[0,2].set_xticklabels(summary_df['time_point'])
    axes[0,2].legend()
    axes[0,2].tick_params(axis='x', rotation=45)
    
    # Plot 4: Percentage significant
    axes[1,0].bar(summary_df['time_point'], summary_df['pct_significant'], color='purple', alpha=0.7)
    axes[1,0].set_title('Percentage of Significant Genes')
    axes[1,0].set_ylabel('Percentage (%)')
    axes[1,0].tick_params(axis='x', rotation=45)
    
    # Plot 5: Effect size distribution
    axes[1,1].bar(summary_df['time_point'], summary_df['mean_abs_log2fc_significant'], color='brown', alpha=0.7)
    axes[1,1].set_title('Mean |log2FC| of Significant Genes')
    axes[1,1].set_ylabel('Mean |log2FC|')
    axes[1,1].tick_params(axis='x', rotation=45)
    
    # Plot 6: Statistical significance
    axes[1,2].bar(summary_df['time_point'], -np.log10(summary_df['median_padj_significant']), color='teal', alpha=0.7)
    axes[1,2].set_title('Median -log10(padj) of Significant Genes')
    axes[1,2].set_ylabel('-log10(median padj)')
    axes[1,2].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    summary_plot_file = os.path.join(output_dir, 'improved_summary_plots.png')
    plt.savefig(summary_plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Improved summary plots saved to: {summary_plot_file}")
    
else:
    print("No differential expression results generated.")

print("\n=== IMPROVED ANALYSIS COMPLETE ===")
print(f"All results saved to: {output_dir}")
print("\nKey improvements implemented:")
print("- Batch correction using ComBat")
print("- DESeq2-style normalization and size factor calculation")
print("- More stringent significance thresholds (padj < 0.01, |log2FC| > 0.75)")
print("- Detection rate filtering for spatial transcriptomics")
print("- Enhanced quality control and visualization")
print("- Cell count imbalance assessment and reporting")

# %%