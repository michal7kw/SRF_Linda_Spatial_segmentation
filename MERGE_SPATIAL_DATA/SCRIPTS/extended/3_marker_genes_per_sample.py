# %%
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
import anndata
import os
import sys

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Configure scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Define the path to the merged AnnData file
merged_h5ad_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_extended.h5ad'

# %%
# Load the merged AnnData object
if os.path.exists(merged_h5ad_path):
    print(f"Loading merged data from: {merged_h5ad_path}")
    adata = sc.read_h5ad(merged_h5ad_path)
    print("Merged data loaded successfully.")
    print(f"Data shape: {adata.shape}")
else:
    print(f"Error: Merged data file not found at {merged_h5ad_path}")
    sys.exit()

# %%
# Create output directory for results
output_dir = 'MERGE_SPATIAL_DATA/marker_genes_per_sample_extended'
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

print(f"Results will be saved to: {output_dir}")

# %%
# Get unique samples
samples = adata.obs['sample_id'].unique()
print(f"Found {len(samples)} samples: {samples}")

# %%
# Analyze each sample separately
# NOTE: Cluster extension has already been done at the dataset level
# by create_extended_roi_dataset.py, so we use the data as-is
for sample_id in samples:
    try:
        print(f"\n{'='*60}")
        print(f"Processing sample: {sample_id}")
        print(f"{'='*60}")
        
        # Create sample-specific output directory
        sample_output_dir = os.path.join(output_dir, f"sample_{sample_id}")
        os.makedirs(sample_output_dir, exist_ok=True)
        
        # Filter data for current sample from the already-extended dataset
        sample_adata = adata[adata.obs['sample_id'] == sample_id].copy()
        print(f"Sample {sample_id} shape: {sample_adata.shape}")
        print("NOTE: Using pre-extended dataset where small clusters have already been")
        print("      augmented with cells from the same clusters across all samples.")
        
        # Use the sample data as-is (extension already done at dataset level)
        sample_adata_filtered = sample_adata
        
        # Update leiden categories to only include present clusters
        sample_adata_filtered.obs['leiden'] = sample_adata_filtered.obs['leiden'].cat.remove_unused_categories()
        
        print(f"Final data shape for analysis: {sample_adata_filtered.shape}")
        print(f"Valid clusters: {sorted(sample_adata_filtered.obs['leiden'].cat.categories.astype(int))}")
        
        # Perform differential expression analysis for this sample
        print(f"Performing differential expression analysis for {sample_id}...")
        
        sc.tl.rank_genes_groups(
            sample_adata_filtered, 
            'leiden', 
            method='wilcoxon',
            key_added='rank_genes_leiden',
            n_genes=50
        )
        
        # Extract marker gene results for this sample
        print(f"Extracting marker gene results for {sample_id}...")
        
        marker_genes_list = []
        for group in sample_adata_filtered.obs['leiden'].cat.categories:
            group_df = sc.get.rank_genes_groups_df(sample_adata_filtered, group=group, key='rank_genes_leiden')
            group_df['group'] = group
            group_df['sample_id'] = sample_id
            marker_genes_list.append(group_df)
        
        sample_marker_genes_df = pd.concat(marker_genes_list, ignore_index=True)
        
        # Save marker genes for this sample
        sample_marker_genes_df.to_csv(os.path.join(sample_output_dir, 'marker_genes.csv'), index=False)
        
        # Create summary of top marker genes per cluster for this sample
        print(f"Creating top marker genes summary for {sample_id}...")
        
        top_markers_per_cluster = {}
        n_top_genes = 30
        
        for cluster in sample_adata_filtered.obs['leiden'].cat.categories:
            cluster_markers = sample_marker_genes_df[sample_marker_genes_df['group'] == cluster].head(n_top_genes)
            top_markers_per_cluster[cluster] = cluster_markers['names'].tolist()
        
        # Save top markers summary
        top_markers_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in top_markers_per_cluster.items()]))
        top_markers_df.to_csv(os.path.join(sample_output_dir, f'top_{n_top_genes}_marker_genes.csv'), index=False)
        
        # Create visualizations for this sample
        print(f"Creating visualizations for {sample_id}...")
        
        # Set figure directory for this sample
        sc.settings.figdir = sample_output_dir
        
        try:
            # Plot marker genes heatmap
            print(f"Creating heatmap for {sample_id}...")
            heatmap_ax_dict = sc.pl.rank_genes_groups_heatmap(
                sample_adata_filtered,
                key='rank_genes_leiden',
                n_genes=5,
                show_gene_labels=True,
                figsize=(10, 7),
                show=False,
            )
            ax = heatmap_ax_dict['heatmap_ax']
            ax.tick_params(axis='x', labelsize=4, rotation=90)
            ax.tick_params(axis='y', labelsize=6)
            fig = ax.get_figure()
            fig.savefig(os.path.join(sample_output_dir, f'heatmap_{sample_id}_marker_genes_heatmap.png'), dpi=300, bbox_inches='tight')
            plt.close('all')
            print(f"Heatmap saved for {sample_id}")
            
        except Exception as e:
            print(f"Warning: Could not create heatmap for {sample_id}: {e}")
            import traceback
            traceback.print_exc()
        
        try:
            # Create dotplot for marker genes
            print(f"Creating dotplot for {sample_id}...")
            dotplot_ax_dict = sc.pl.rank_genes_groups_dotplot(
                sample_adata_filtered,
                key='rank_genes_leiden',
                n_genes=5,
                figsize=(12, 7),
                show=False,
            )
            print(f"Dotplot axes dict: {dotplot_ax_dict}")
            print(f"Available keys: {list(dotplot_ax_dict.keys()) if dotplot_ax_dict else 'None'}")
            
            # Try to find a valid axes object
            ax = None
            ax_key = None
            
            if dotplot_ax_dict:
                # First try mainplot_ax
                if 'mainplot_ax' in dotplot_ax_dict and dotplot_ax_dict['mainplot_ax'] is not None:
                    ax = dotplot_ax_dict['mainplot_ax']
                    ax_key = 'mainplot_ax'
                    print(f"Using mainplot_ax for {sample_id}")
                # Try alternative keys if mainplot_ax doesn't work
                elif 'ax' in dotplot_ax_dict and dotplot_ax_dict['ax'] is not None:
                    ax = dotplot_ax_dict['ax']
                    ax_key = 'ax'
                    print(f"Using ax for {sample_id}")
                else:
                    # Try to find any axes object in the dictionary
                    for key, value in dotplot_ax_dict.items():
                        if hasattr(value, 'tick_params') and hasattr(value, 'get_figure'):  # Check if it's an axes object
                            ax = value
                            ax_key = key
                            print(f"Using {key} as axes for {sample_id}")
                            break
            
            if ax is not None:
                try:
                    ax.tick_params(axis='x', labelsize=4, rotation=90)
                    ax.tick_params(axis='y', labelsize=6)
                    fig = ax.get_figure()
                    fig.savefig(os.path.join(sample_output_dir, f'dotplot_{sample_id}_marker_genes_dotplot.png'), dpi=300, bbox_inches='tight')
                    plt.close('all')
                    print(f"Dotplot saved for {sample_id} using {ax_key}")
                except Exception as axes_error:
                    print(f"Error manipulating axes for {sample_id}: {axes_error}")
                    # Try to save without modifications
                    try:
                        fig = ax.get_figure()
                        fig.savefig(os.path.join(sample_output_dir, f'dotplot_{sample_id}_marker_genes_dotplot.png'), dpi=300, bbox_inches='tight')
                        plt.close('all')
                        print(f"Dotplot saved for {sample_id} without tick modifications")
                    except Exception as save_error:
                        print(f"Could not save dotplot for {sample_id}: {save_error}")
            else:
                print(f"No valid axes found in dotplot result for {sample_id}")
                # Debug: print the type and content of each item
                if dotplot_ax_dict:
                    for key, value in dotplot_ax_dict.items():
                        print(f"  {key}: {type(value)} - {value}")
                
        except Exception as e:
            print(f"Warning: Could not create dotplot for {sample_id}: {e}")
            import traceback
            traceback.print_exc()
        
        # Create dendrogram for cluster similarity within this sample
        print(f"Creating dendrogram for {sample_id}...")
        
        # Calculate mean expression per cluster
        cluster_means = []
        cluster_labels = []
        
        for cluster in sorted(sample_adata_filtered.obs['leiden'].cat.categories, key=int):
            cluster_cells = sample_adata_filtered[sample_adata_filtered.obs['leiden'] == cluster]
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
            plt.figure(figsize=(10, 7))
            dendrogram_plot = dendrogram(
                linkage_matrix,
                labels=cluster_labels,
                leaf_rotation=45,
                leaf_font_size=8
            )
            
            plt.title(f'Cluster Similarity Dendrogram - {sample_id}\n(Based on Gene Expression Profiles)', fontsize=10)
            plt.xlabel('Clusters', fontsize=8)
            plt.ylabel('Distance (Ward Linkage)', fontsize=8)
            plt.tight_layout()
            
            # Save the dendrogram
            dendrogram_path = os.path.join(sample_output_dir, f'{sample_id}_cluster_dendrogram.png')
            plt.savefig(dendrogram_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            # Create correlation heatmap
            correlation_matrix = np.corrcoef(cluster_means_matrix)
            
            plt.figure(figsize=(10, 7))
            mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
            sns.heatmap(
                correlation_matrix,
                mask=mask,
                annot=True, annot_kws={"size": 8},
                cmap='coolwarm',
                center=0,
                square=True,
                xticklabels=cluster_labels,
                yticklabels=cluster_labels,
                cbar_kws={"shrink": .8}
            )
            
            plt.title(f'Cluster Correlation Heatmap - {sample_id}', fontsize=10)
            plt.xticks(rotation=45, fontsize=8)
            plt.yticks(rotation=0, fontsize=8)
            plt.tight_layout()
            
            correlation_heatmap_path = os.path.join(sample_output_dir, f'{sample_id}_cluster_correlation.png')
            plt.savefig(correlation_heatmap_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            # Save correlation matrix
            correlation_df = pd.DataFrame(
                correlation_matrix,
                index=cluster_labels,
                columns=cluster_labels
            )
            correlation_df.to_csv(os.path.join(sample_output_dir, f'{sample_id}_correlation_matrix.csv'))
            
            print(f"Dendrogram and correlation analysis completed for {sample_id}")
        else:
            print(f"Skipping dendrogram for {sample_id}: insufficient clusters")
        
        # Create sample summary report
        sample_summary = f"""Sample: {sample_id}
================

Dataset Information:
- Total cells in sample: {sample_adata.n_obs:,}
- Cells after extension: {sample_adata_filtered.n_obs:,}
- Total genes: {sample_adata_filtered.n_vars:,}
- Total clusters in sample: {len(sample_adata.obs['leiden'].cat.categories)}
- Clusters after extension: {len(sample_adata_filtered.obs['leiden'].cat.categories)}

Files Generated:
1. marker_genes.csv - Complete marker gene results
2. top_{n_top_genes}_marker_genes.csv - Top marker genes per cluster
3. heatmap_{sample_id}_marker_genes_heatmap.png - Marker genes heatmap
4. dotplot_{sample_id}_marker_genes_dotplot.png - Marker genes dotplot
5. {sample_id}_cluster_dendrogram.png - Hierarchical clustering dendrogram
6. {sample_id}_cluster_correlation.png - Cluster correlation heatmap
7. {sample_id}_correlation_matrix.csv - Correlation matrix data
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
        with open(os.path.join(sample_output_dir, f'{sample_id}_error.log'), 'w') as f:
            f.write(error_log)
        continue

# %%
# Create overall summary
print(f"\n{'='*60}")
print("Creating overall summary...")
print(f"{'='*60}")

# Collect statistics from all samples
overall_stats = []
for sample_id in samples:
    sample_adata = adata[adata.obs['sample_id'] == sample_id]
    cluster_counts = sample_adata.obs['leiden'].value_counts()
    small_clusters = cluster_counts[(cluster_counts > 0) & (cluster_counts < 50)]
    
    overall_stats.append({
        'sample_id': sample_id,
        'total_cells': sample_adata.n_obs,
        'total_clusters': len(cluster_counts),
        'small_clusters': len(small_clusters),
        'small_cluster_list': list(small_clusters.index),
        'cells_in_small_clusters': small_clusters.sum()
    })

overall_stats_df = pd.DataFrame(overall_stats)
overall_stats_df.to_csv(os.path.join(output_dir, 'overall_sample_statistics.csv'), index=False)

# Create overall summary report
overall_summary = f"""Marker Genes and Cluster Similarity Analysis - Overall Summary
============================================================

Analysis completed for {len(samples)} samples:
{chr(10).join([f"- {sample}" for sample in samples])}

Overall Statistics:
- Total cells across all samples: {adata.n_obs:,}
- Total genes: {adata.n_vars:,}
- Samples processed: {len(samples)}

Per-Sample Results:
{chr(10).join([f"- {row['sample_id']}: {row['total_cells']} cells, {row['small_clusters']} small clusters extended" for _, row in overall_stats_df.iterrows()])}

Analysis Parameters:
- Differential expression method: Wilcoxon rank-sum test
- Number of marker genes per cluster: 50
- Distance metric for clustering: Correlation distance
- Linkage method: Ward
- Small cluster threshold: < 50 cells (extended with all other cells from the same cluster)
"""

with open(os.path.join(output_dir, 'overall_analysis_summary.txt'), 'w') as f:
    f.write(overall_summary)

print(overall_summary)
print(f"\nAnalysis completed successfully!")
print(f"Results saved to: {output_dir}")
print(f"Each sample has its own subdirectory with individual results.")

# %%