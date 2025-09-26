# %%
# Suppress specific warnings
import warnings
warnings.filterwarnings('ignore', category=FutureWarning, message='.*The previous implementation of stack is deprecated.*')
warnings.filterwarnings('ignore', category=FutureWarning, message='.*The default of observed=False is deprecated.*')
warnings.filterwarnings('ignore', category=UserWarning, message='.*No data for colormapping provided.*')

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
import os
import sys

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Configure scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Define the path to the precomputed reclustered data file (produced by 1_correlation_analysis.py)
reclustered_h5ad_path = 'MERGE_SPATIAL_DATA/5_RECLUSTERING/correlation_analysis_cluster0_recluster/merged_data_with_combined_cluster0_reclustering.h5ad'

# %%
# Load the precomputed reclustered AnnData object
print(f"Loading precomputed reclustered data from: {reclustered_h5ad_path}")
adata = sc.read_h5ad(reclustered_h5ad_path)
print("Precomputed reclustered data loaded successfully.")
print(f"Data shape: {adata.shape}")
print(f"Available cluster columns: {[col for col in adata.obs.columns if 'leiden' in col]}")

# %%
# Create output directory for results
output_dir = 'MERGE_SPATIAL_DATA/marker_genes_per_sample_inclusive_cluster0_recluster'
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

print(f"Results will be saved to: {output_dir}")

# %%
# Get unique samples
samples = adata.obs['sample_id'].unique()
print(f"Found {len(samples)} samples: {samples}")

# %%
# Use precomputed reclustered data
main_cluster_column = 'leiden_combined_cluster0_reclustered'
print(f"Using precomputed reclustered data with column: {main_cluster_column}")

# %%
# Function to create UMAP plots with spatial overlay
def create_umap_spatial_plots(sample_adata, cluster_col, sample_id, output_dir, suffix=""):
    """
    Create side-by-side UMAP and spatial plots
    """
    try:
        # Create figure with subplots
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))
        
        # UMAP plot
        sc.pl.umap(sample_adata, color=cluster_col, ax=axes[0], show=False, 
                  legend_loc='on data', legend_fontsize=8, s=30)
        axes[0].set_title(f'UMAP - {sample_id}{suffix}', fontsize=12)
        
        # Spatial plot
        if 'spatial' in sample_adata.obsm.keys():
            sc.pl.spatial(sample_adata, color=cluster_col, ax=axes[1], show=False,
                         legend_loc='right margin', legend_fontsize=8, size=1.5)
            axes[1].set_title(f'Spatial - {sample_id}{suffix}', fontsize=12)
        else:
            # Fallback to scatter plot if spatial coordinates not available
            axes[1].scatter(sample_adata.obs['center_x'], sample_adata.obs['center_y'], 
                           c=sample_adata.obs[cluster_col].cat.codes, s=10, cmap='tab20')
            axes[1].set_title(f'Spatial (center coordinates) - {sample_id}{suffix}', fontsize=12)
            axes[1].set_xlabel('Center X')
            axes[1].set_ylabel('Center Y')
        
        plt.tight_layout()
        
        # Save the plot
        plot_filename = f'{sample_id}_umap_spatial{suffix.replace(" ", "_")}.png'
        plt.savefig(os.path.join(output_dir, plot_filename), dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"    UMAP/Spatial plot saved: {plot_filename}")
        
    except Exception as e:
        print(f"    Warning: Could not create UMAP/Spatial plot for {sample_id}{suffix}: {e}")
        plt.close('all')

# %%
# Analyze each sample separately
for sample_id in samples:
    try:
        print(f"\n{'='*60}")
        print(f"Processing sample: {sample_id}")
        print(f"{'='*60}")
        
        # Create sample-specific output directory
        sample_output_dir = os.path.join(output_dir, f"sample_{sample_id}")
        os.makedirs(sample_output_dir, exist_ok=True)
        
        # Filter data for current sample
        sample_adata = adata[adata.obs['sample_id'] == sample_id].copy()
        print(f"Sample {sample_id} shape: {sample_adata.shape}")
        
        # Update cluster categories to only include present clusters
        if main_cluster_column in sample_adata.obs.columns:
            sample_adata.obs[main_cluster_column] = sample_adata.obs[main_cluster_column].cat.remove_unused_categories()
        
        # Get cluster counts
        cluster_counts = sample_adata.obs[main_cluster_column].value_counts()
        
        print(f"\nUsing precomputed reclustered data for sample {sample_id}...")
        # Also get original cluster counts if available for comparison
        original_cluster_counts = sample_adata.obs['leiden_new'].value_counts() if 'leiden_new' in sample_adata.obs.columns else pd.Series()
        if len(original_cluster_counts) > 0:
            print(f"Sample {sample_id} original cluster counts:")
            for cluster, count in original_cluster_counts.sort_index().items():
                print(f"  Cluster {cluster}: {count} cells")
            
        print(f"Sample {sample_id} cluster counts ({main_cluster_column}):")
        for cluster, count in cluster_counts.sort_index().items():
            print(f"  Cluster {cluster}: {count} cells")
        
        # Check if cluster 0 was actually reclustered (has subclusters)
        has_cluster0_subclusters = any(str(cluster).startswith('0.') for cluster in cluster_counts.index)
        original_had_cluster0 = ('0' in sample_adata.obs['leiden_new'].value_counts().index) if 'leiden_new' in sample_adata.obs.columns else False
        
        # Create UMAP/Spatial plots for original clustering (if original data is available)
        if 'leiden_new' in sample_adata.obs.columns:
            print(f"Creating original clustering UMAP/Spatial plots for {sample_id}...")
            create_umap_spatial_plots(sample_adata, 'leiden_new', sample_id, sample_output_dir, "Original")
        
        # Create UMAP/Spatial plots for main clustering data
        plot_suffix = "Precomputed_Reclustered"
        print(f"Creating UMAP/Spatial plots for {sample_id} using {main_cluster_column}...")
        create_umap_spatial_plots(sample_adata, main_cluster_column, sample_id, sample_output_dir, plot_suffix)
        
        # Create subset plot for cluster 0 subclusters if they exist
        if has_cluster0_subclusters:
            print(f"Creating cluster 0 subclusters subset plots for {sample_id}...")
            
            # Filter to only cluster 0 subclusters
            cluster0_subclusters = [c for c in cluster_counts.index if str(c).startswith('0.')]
            cluster0_mask = sample_adata.obs[main_cluster_column].isin(cluster0_subclusters)
            cluster0_subset = sample_adata[cluster0_mask].copy()
            
            if len(cluster0_subset) > 0:
                # Plot cluster 0 subset
                fig, axes = plt.subplots(1, 2, figsize=(16, 7))
                
                # UMAP of cluster 0 subset with subclusters
                sc.pl.umap(cluster0_subset, color=main_cluster_column, ax=axes[0], show=False,
                          legend_loc='on data', legend_fontsize=8, s=50)
                axes[0].set_title(f'Cluster 0 Subclusters UMAP - {sample_id}\n(Precomputed)', fontsize=12)
                
                # Spatial plot of cluster 0 subset
                if 'spatial' in cluster0_subset.obsm.keys():
                    sc.pl.spatial(cluster0_subset, color=main_cluster_column, ax=axes[1], show=False,
                                 legend_loc='right margin', legend_fontsize=8, size=2.0)
                    axes[1].set_title(f'Cluster 0 Subclusters Spatial - {sample_id}\n(Precomputed)', fontsize=12)
                else:
                    # Use cluster codes for coloring
                    cluster_codes = cluster0_subset.obs[main_cluster_column].cat.codes
                    axes[1].scatter(cluster0_subset.obs['center_x'], cluster0_subset.obs['center_y'], 
                                   c=cluster_codes, s=15, cmap='tab10')
                    axes[1].set_title(f'Cluster 0 Subclusters Spatial - {sample_id}\n(Precomputed)', fontsize=12)
                    axes[1].set_xlabel('Center X')
                    axes[1].set_ylabel('Center Y')
                
                plt.tight_layout()
                subset_plot_filename = f'{sample_id}_subclusters_umap_spatial.png'
                plt.savefig(os.path.join(sample_output_dir, subset_plot_filename), dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"    Cluster 0 subclusters plot saved: {subset_plot_filename}")
        
        # Use precomputed reclustered column for analysis
        cluster_column = main_cluster_column
        final_cluster_counts = cluster_counts
        valid_clusters_for_de = final_cluster_counts[final_cluster_counts > 1].index.tolist()
        
        print(f"\nFinal clustering analysis for {sample_id}:")
        print(f"  Total clusters: {len(final_cluster_counts)}")
        print(f"  Clusters with >1 cell for DE analysis: {len(valid_clusters_for_de)} - {valid_clusters_for_de}")

        # Initialize these to empty in case DE analysis is skipped
        sample_marker_genes_df = pd.DataFrame(columns=['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges', 'group', 'sample_id'])
        top_markers_df = pd.DataFrame()
        n_top_genes = 30

        if len(valid_clusters_for_de) < 2:
            print(f"Skipping differential expression analysis for {sample_id}: insufficient clusters with >1 cell for comparison.")
            sample_adata.uns[f'rank_genes_{cluster_column}'] = {'names': np.array([]), 'scores': np.array([]), 'pvals': np.array([]), 'pvals_adj': np.array([]), 'logfoldchanges': np.array([])}
        else:
            # Perform differential expression analysis
            print(f"Performing differential expression analysis for {sample_id} using {cluster_column}...")
            
            sc.tl.rank_genes_groups(
                sample_adata, 
                cluster_column, 
                groups=valid_clusters_for_de,
                method='wilcoxon',
                key_added=f'rank_genes_{cluster_column}',
                n_genes=50
            )
            
            # Extract marker gene results
            print(f"Extracting marker gene results for {sample_id}...")
            
            marker_genes_list = []
            for group in valid_clusters_for_de: 
                group_df = sc.get.rank_genes_groups_df(sample_adata, group=group, key=f'rank_genes_{cluster_column}')
                group_df['group'] = group
                group_df['sample_id'] = sample_id
                marker_genes_list.append(group_df)
            
            if marker_genes_list:
                sample_marker_genes_df = pd.concat(marker_genes_list, ignore_index=True)
            
            # Save marker genes results
            marker_genes_filename = 'marker_genes.csv'
            sample_marker_genes_df.to_csv(os.path.join(sample_output_dir, marker_genes_filename), index=False)
            
            # Create summary of top marker genes per cluster
            print(f"Creating top marker genes summary for {sample_id}...")
            
            top_markers_per_cluster = {}
            for cluster in valid_clusters_for_de:
                cluster_markers = sample_marker_genes_df[sample_marker_genes_df['group'] == cluster].head(n_top_genes)
                top_markers_per_cluster[cluster] = cluster_markers['names'].tolist()
            
            # Create top markers summary
            if top_markers_per_cluster:
                top_markers_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in top_markers_per_cluster.items()]))
                top_markers_filename = f'top_{n_top_genes}_marker_genes.csv'
                top_markers_df.to_csv(os.path.join(sample_output_dir, top_markers_filename), index=False)
            
            # Create visualizations
            print(f"Creating visualizations for {sample_id}...")
            
            # Set figure directory for this sample
            sc.settings.figdir = sample_output_dir
            
            try:
                # Plot marker genes heatmap
                print(f"Creating heatmap for {sample_id}...")
                heatmap_ax_dict = sc.pl.rank_genes_groups_heatmap(
                    sample_adata,
                    key=f'rank_genes_{cluster_column}',
                    n_genes=5,
                    show_gene_labels=True,
                    figsize=(12, 8),
                    show=False,
                    dendrogram=False,  # Disable dendrogram to avoid category mismatch
                )
                ax = heatmap_ax_dict['heatmap_ax']
                ax.tick_params(axis='x', labelsize=4, rotation=90)
                ax.tick_params(axis='y', labelsize=6)
                fig = ax.get_figure()
                heatmap_filename = f'heatmap_{sample_id}_marker_genes.png'
                fig.savefig(os.path.join(sample_output_dir, heatmap_filename), dpi=300, bbox_inches='tight')
                plt.close('all')
                print(f"Heatmap saved: {heatmap_filename}")
                
            except Exception as e:
                print(f"Warning: Could not create heatmap for {sample_id}: {e}")
                import traceback
                traceback.print_exc()
            
            try:
                # Create dotplot for marker genes
                print(f"Creating dotplot for {sample_id}...")
                dotplot_ax_dict = sc.pl.rank_genes_groups_dotplot(
                    sample_adata,
                    key=f'rank_genes_{cluster_column}',
                    n_genes=5,
                    figsize=(14, 8),
                    show=False,
                    dendrogram=False
                )
                
                ax = None
                ax_key = None
                
                if dotplot_ax_dict:
                    if 'mainplot_ax' in dotplot_ax_dict and dotplot_ax_dict['mainplot_ax'] is not None:
                        ax = dotplot_ax_dict['mainplot_ax']
                        ax_key = 'mainplot_ax'
                    elif 'ax' in dotplot_ax_dict and dotplot_ax_dict['ax'] is not None:
                        ax = dotplot_ax_dict['ax']
                        ax_key = 'ax'
                    else:
                        for key, value in dotplot_ax_dict.items():
                            if hasattr(value, 'tick_params') and hasattr(value, 'get_figure'):
                                ax = value
                                ax_key = key
                                break
                
                if ax is not None:
                    try:
                        ax.tick_params(axis='x', labelsize=4, rotation=90)
                        ax.tick_params(axis='y', labelsize=6)
                        fig = ax.get_figure()
                        dotplot_filename = f'dotplot_{sample_id}_marker_genes.png'
                        fig.savefig(os.path.join(sample_output_dir, dotplot_filename), dpi=300, bbox_inches='tight')
                        plt.close('all')
                        print(f"Dotplot saved: {dotplot_filename} using {ax_key}")
                    except Exception as axes_error:
                        print(f"Error manipulating axes for {sample_id}: {axes_error}")
                
            except Exception as e:
                print(f"Warning: Could not create dotplot for {sample_id}: {e}")
                import traceback
                traceback.print_exc()
            
            # Create dendrogram for cluster similarity
            print(f"Creating dendrogram for {sample_id}...")
            
            # Calculate mean expression per cluster
            cluster_means = []
            cluster_labels = []
            
            for cluster in sorted(valid_clusters_for_de, key=lambda x: (len(x.split('.')), x)):
                cluster_cells = sample_adata[sample_adata.obs[cluster_column] == cluster]
                if len(cluster_cells) > 0:
                    mean_expr = np.array(cluster_cells.X.mean(axis=0)).flatten()
                    cluster_means.append(mean_expr)
                    cluster_labels.append(f"Cluster_{cluster}")
            
            if len(cluster_means) >= 2:
                cluster_means_matrix = np.array(cluster_means)
                
                # Calculate pairwise distances between clusters
                distances = pdist(cluster_means_matrix, metric='correlation')
                linkage_matrix = linkage(distances, method='ward')
                
                # Create dendrogram
                plt.figure(figsize=(12, 8))
                dendrogram(
                    linkage_matrix,
                    labels=cluster_labels,
                    leaf_rotation=45,
                    leaf_font_size=8
                )
                
                title_suffix = " (with Cluster 0 Reclustered)"
                plt.title(f'Cluster Similarity Dendrogram - {sample_id}{title_suffix}\n(Based on Gene Expression Profiles)', fontsize=11)
                plt.xlabel('Clusters', fontsize=9)
                plt.ylabel('Distance (Ward Linkage)', fontsize=9)
                plt.tight_layout()
                
                # Save dendrogram
                dendrogram_filename = f'{sample_id}_cluster_dendrogram.png'
                plt.savefig(os.path.join(sample_output_dir, dendrogram_filename), dpi=300, bbox_inches='tight')
                plt.close()
                
                # Create correlation heatmap
                correlation_matrix = np.corrcoef(cluster_means_matrix)
                
                plt.figure(figsize=(12, 8))
                mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
                sns.heatmap(
                    correlation_matrix,
                    mask=mask,
                    annot=True, annot_kws={"size": 6},
                    cmap='coolwarm',
                    center=0,
                    square=True,
                    xticklabels=cluster_labels,
                    yticklabels=cluster_labels,
                    cbar_kws={"shrink": .8}
                )
                
                plt.title(f'Cluster Correlation Heatmap - {sample_id}{title_suffix}', fontsize=11)
                plt.xticks(rotation=45, fontsize=8)
                plt.yticks(rotation=0, fontsize=8)
                plt.tight_layout()
                
                correlation_filename = f'{sample_id}_cluster_correlation.png'
                plt.savefig(os.path.join(sample_output_dir, correlation_filename), dpi=300, bbox_inches='tight')
                plt.close()
                
                # Save correlation matrix as CSV
                correlation_df = pd.DataFrame(
                    correlation_matrix,
                    index=cluster_labels,
                    columns=cluster_labels
                )
                correlation_csv_filename = f'{sample_id}_correlation_matrix.csv'
                correlation_df.to_csv(os.path.join(sample_output_dir, correlation_csv_filename))
                
                print(f"Dendrogram and correlation analysis completed for {sample_id}")
            else:
                print(f"Skipping dendrogram for {sample_id}: insufficient clusters with >1 cell for correlation analysis")
            
            # Create sample summary report
            cluster0_info = ""
            if has_cluster0_subclusters:
                n_subclusters = len([c for c in final_cluster_counts.index if str(c).startswith('0.')])
                original_cluster0_cells = original_cluster_counts.get('0', 0) if len(original_cluster_counts) > 0 else 0
                cluster0_info = f"""
Cluster 0 Reclustering (Precomputed):
- Original cluster 0 cells: {original_cluster0_cells}
- Number of subclusters created: {n_subclusters}
- Subcluster labels: {sorted([c for c in final_cluster_counts.index if str(c).startswith('0.')])}
"""
            
            sample_summary = f"""Sample: {sample_id}
================
Dataset Information:
- Total cells in sample: {sample_adata.n_obs:,}
- Total genes: {sample_adata.n_vars:,}
- Original total clusters found: {len(cluster_counts)}
- Final total clusters (after reclustering): {len(final_cluster_counts)}
- Clusters analyzed for DE (with >1 cell): {len(valid_clusters_for_de)}
{cluster0_info}
Files Generated:
1. {marker_genes_filename} - Complete marker gene results
2. {top_markers_filename} - Top marker genes per cluster
3. Various UMAP/Spatial visualization plots
4. {heatmap_filename} - Marker genes heatmap
5. {dotplot_filename} - Marker genes dotplot
6. {dendrogram_filename} - Hierarchical clustering dendrogram
7. {correlation_filename} - Cluster correlation heatmap
8. {correlation_csv_filename} - Correlation matrix data
"""
            
            with open(os.path.join(sample_output_dir, f'{sample_id}_summary.txt'), 'w') as f:
                f.write(sample_summary)
            
            print(f"Sample {sample_id} analysis completed!")
        
    except Exception as e:
        print(f"ERROR processing sample {sample_id}: {str(e)}")
        import traceback
        traceback.print_exc()
        
        # Create error log
        error_log = f"Error processing sample {sample_id}:\n{str(e)}\n\n{traceback.format_exc()}"
        error_log_path = os.path.join(sample_output_dir, f'{sample_id}_error.log')
        with open(error_log_path, 'w') as f:
            f.write(error_log)
        continue

# %%
# Create overall summary
from datetime import datetime
print(f"\n{'='*60}")
print("Creating overall summary...")
print(f"{'='*60}")

# Generate summary
print("\n" + "="*80)
print("MARKER GENES ANALYSIS SUMMARY (PRECOMPUTED RECLUSTERED)")
print("="*80)

summary_content = []
summary_content.append("MARKER GENES ANALYSIS SUMMARY (PRECOMPUTED RECLUSTERED)")
summary_content.append("="*80)
summary_content.append(f"Analysis completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
summary_content.append(f"Input file: {reclustered_h5ad_path}")
summary_content.append(f"Output directory: {output_dir}")
summary_content.append(f"Total samples processed: {len(samples)}")
summary_content.append("")
summary_content.append("CLUSTERING INFORMATION:")
summary_content.append("- Used precomputed reclustered data from correlation analysis")
summary_content.append(f"- Cluster column: {main_cluster_column}")
summary_content.append("- This includes refined subclusters for original cluster 0")

# Collect statistics from all samples
overall_stats = []
for sample_id in samples:
    sample_adata = adata[adata.obs['sample_id'] == sample_id].copy()
    
    # Get original cluster counts if available
    if 'leiden_new' in sample_adata.obs.columns:
        sample_adata.obs['leiden_new'] = sample_adata.obs['leiden_new'].cat.remove_unused_categories()
        original_cluster_counts = sample_adata.obs['leiden_new'].value_counts()
    else:
        original_cluster_counts = pd.Series()
    
    # Use precomputed reclustered data
    sample_adata.obs[main_cluster_column] = sample_adata.obs[main_cluster_column].cat.remove_unused_categories()
    final_cluster_counts = sample_adata.obs[main_cluster_column].value_counts()
    valid_clusters_for_de_summary = final_cluster_counts[final_cluster_counts > 1].index.tolist()
    
    # Check if cluster 0 was reclustered (has subclusters)
    has_cluster0_subclusters = any(str(cluster).startswith('0.') for cluster in final_cluster_counts.index)
    original_had_cluster0 = '0' in original_cluster_counts.index if len(original_cluster_counts) > 0 else False

    overall_stats.append({
        'sample_id': sample_id,
        'total_cells': sample_adata.n_obs,
        'original_total_clusters': len(original_cluster_counts) if len(original_cluster_counts) > 0 else 0,
        'final_total_clusters': len(final_cluster_counts),
        'clusters_analyzed_for_de': len(valid_clusters_for_de_summary) if len(valid_clusters_for_de_summary) >= 2 else 0,
        'had_cluster0': original_had_cluster0,
        'cluster0_reclustered': has_cluster0_subclusters,
    })

overall_stats_df = pd.DataFrame(overall_stats)
stats_filename = 'overall_sample_statistics_with_cluster0_reclustering.csv'
overall_stats_df.to_csv(os.path.join(output_dir, stats_filename), index=False)

# Create overall summary report
title = "Marker Genes and Cluster Similarity Analysis with Precomputed Cluster 0 Reclustering - Overall Summary"
data_source = "Precomputed reclustered data from previous pipeline step"
cluster_info = f"- Samples with cluster 0 subclusters (precomputed): {sum(overall_stats_df['cluster0_reclustered'])}"
output_features = """- Original UMAP/Spatial plots showing all clusters (when available)
- Precomputed reclustered UMAP/Spatial plots
- Subset UMAP/Spatial plots showing only cluster 0 subclusters (when present)
- Marker gene analysis for all final clusters (including subclusters)
- Hierarchical clustering and correlation analysis"""

overall_summary = f"""{title}
{'='*len(title)}
Analysis completed for {len(samples)} samples:
{chr(10).join([f"- {sample}" for sample in samples])}

Overall Statistics:
- Total cells across all samples: {adata.n_obs:,}
- Total genes: {adata.n_vars:,}
- Samples processed: {len(samples)}
{cluster_info}

Per-Sample Results:
{chr(10).join([f"- {row['sample_id']}: {row['total_cells']} cells, {row['original_total_clusters']} original clusters → {row['final_total_clusters']} final clusters, {row['clusters_analyzed_for_de']} analyzed for DE" for _, row in overall_stats_df.iterrows()])}

Analysis Parameters:
- Differential expression method: Wilcoxon rank-sum test
- Number of marker genes per cluster: 50
- Distance metric for clustering: Correlation distance
- Linkage method: Ward
- Minimum cells per cluster threshold for DE: >1 cell
- Data source: {data_source}
- Cluster column used: {main_cluster_column}

Output Features:
{output_features}
"""

summary_filename = 'overall_summary_with_cluster0_reclustering.txt'
with open(os.path.join(output_dir, summary_filename), 'w') as f:
    f.write(overall_summary)

print(overall_summary)
print(f"\nAnalysis completed successfully!")
print(f"Results saved to: {output_dir}")
print(f"Each sample has its own subdirectory with individual results.")

# %%