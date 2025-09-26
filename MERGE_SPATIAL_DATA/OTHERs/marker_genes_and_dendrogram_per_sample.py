# %%
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
sc.settings.verbosity = 3  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Define the path to the merged AnnData file
merged_h5ad_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo.h5ad'

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
output_dir = 'MERGE_SPATIAL_DATA/marker_genes_per_sample'
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

print(f"Results will be saved to: {output_dir}")

# %%
# Get unique samples
samples = adata.obs['sample_id'].unique()
print(f"Found {len(samples)} samples: {samples}")

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
        
        # Check cluster composition and filter clusters with ≥50 cells
        cluster_counts = sample_adata.obs['leiden'].value_counts()
        valid_clusters = cluster_counts[cluster_counts >= 50]
        
        print(f"Clusters with ≥50 cells in {sample_id}:")
        for cluster, count in valid_clusters.items():
            print(f"  Cluster {cluster}: {count} cells")
        
        if len(valid_clusters) < 2:
            print(f"Skipping {sample_id}: insufficient clusters with ≥50 cells")
            continue
    
        # Filter to only include valid clusters
        valid_cluster_mask = sample_adata.obs['leiden'].isin(valid_clusters.index)
        sample_adata_filtered = sample_adata[valid_cluster_mask].copy()
        
        # Update leiden categories to only include present clusters
        sample_adata_filtered.obs['leiden'] = sample_adata_filtered.obs['leiden'].cat.remove_unused_categories()
        
        print(f"Filtered data shape: {sample_adata_filtered.shape}")
        print(f"Valid clusters: {sorted(sample_adata_filtered.obs['leiden'].cat.categories.astype(int))}")
        
        # %%
        # Perform differential expression analysis for this sample
        print(f"Performing differential expression analysis for {sample_id}...")
        
        sc.tl.rank_genes_groups(
            sample_adata_filtered,
            'leiden',
            method='wilcoxon',
            key_added='rank_genes_leiden',
            n_genes=50  # Top 50 genes per cluster
        )
        
        # %%
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
        
        # %%
        # Create summary of top marker genes per cluster for this sample
        print(f"Creating top marker genes summary for {sample_id}...")
        
        top_markers_per_cluster = {}
        n_top_genes = 10
        
        for cluster in sample_adata_filtered.obs['leiden'].cat.categories:
            cluster_markers = sample_marker_genes_df[sample_marker_genes_df['group'] == cluster].head(n_top_genes)
            top_markers_per_cluster[cluster] = cluster_markers['names'].tolist()
        
        # Save top markers summary
        top_markers_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in top_markers_per_cluster.items()]))
        top_markers_df.to_csv(os.path.join(sample_output_dir, f'top_{n_top_genes}_marker_genes.csv'), index=False)
        
        # %%
        # Create visualizations for this sample
        print(f"Creating visualizations for {sample_id}...")
        
        # Set figure directory for this sample
        sc.settings.figdir = sample_output_dir
        
        try:
            # Plot marker genes heatmap
            sc.pl.rank_genes_groups_heatmap(
                sample_adata_filtered,
                key='rank_genes_leiden',
                n_genes=5,
                show_gene_labels=True,
                figsize=(10, 6),
                save=f'_{sample_id}_marker_genes_heatmap.png',
                show=False
            )
            
            # Create dotplot for marker genes
            sc.pl.rank_genes_groups_dotplot(
                sample_adata_filtered,
                key='rank_genes_leiden',
                n_genes=5,
                figsize=(12, 6),
                save=f'_{sample_id}_marker_genes_dotplot.png',
                show=False
            )
        except Exception as e:
            print(f"Warning: Could not create scanpy plots for {sample_id}: {e}")
        
        # %%
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
            plt.figure(figsize=(10, 6))
            dendrogram_plot = dendrogram(
                linkage_matrix,
                labels=cluster_labels,
                leaf_rotation=45,
                leaf_font_size=10
            )
            
            plt.title(f'Cluster Similarity Dendrogram - {sample_id}\n(Based on Gene Expression Profiles)', fontsize=12)
            plt.xlabel('Clusters', fontsize=10)
            plt.ylabel('Distance (Ward Linkage)', fontsize=10)
            plt.tight_layout()
            
            # Save the dendrogram
            dendrogram_path = os.path.join(sample_output_dir, f'{sample_id}_cluster_dendrogram.png')
            plt.savefig(dendrogram_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            # Create correlation heatmap
            correlation_matrix = np.corrcoef(cluster_means_matrix)
            
            plt.figure(figsize=(8, 6))
            mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
            sns.heatmap(
                correlation_matrix,
                mask=mask,
                annot=True,
                cmap='coolwarm',
                center=0,
                square=True,
                xticklabels=cluster_labels,
                yticklabels=cluster_labels,
                cbar_kws={"shrink": .8}
            )
            
            plt.title(f'Cluster Correlation Heatmap - {sample_id}', fontsize=12)
            plt.xticks(rotation=45)
            plt.yticks(rotation=0)
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
        
        # %%
        # Create sample summary report
        sample_summary = f"""
Sample: {sample_id}
================

Dataset Information:
- Total cells in sample: {sample_adata.n_obs:,}
- Cells after filtering (≥50 cells/cluster): {sample_adata_filtered.n_obs:,}
- Total genes: {sample_adata_filtered.n_vars:,}
- Total clusters: {len(sample_adata.obs['leiden'].cat.categories)}
- Valid clusters (≥50 cells): {len(valid_clusters)}

Valid Clusters:
{chr(10).join([f"- Cluster {cluster}: {count} cells" for cluster, count in valid_clusters.items()])}

Files Generated:
1. marker_genes.csv - Complete marker gene results
2. top_{n_top_genes}_marker_genes.csv - Top marker genes per cluster
3. {sample_id}_marker_genes_heatmap.png - Marker genes heatmap
4. {sample_id}_marker_genes_dotplot.png - Marker genes dotplot
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
    valid_clusters = cluster_counts[cluster_counts >= 50]
    
    overall_stats.append({
        'sample_id': sample_id,
        'total_cells': sample_adata.n_obs,
        'total_clusters': len(cluster_counts),
        'valid_clusters': len(valid_clusters),
        'valid_cluster_list': list(valid_clusters.index),
        'cells_in_valid_clusters': valid_clusters.sum()
    })

overall_stats_df = pd.DataFrame(overall_stats)
overall_stats_df.to_csv(os.path.join(output_dir, 'overall_sample_statistics.csv'), index=False)

# Create overall summary report
overall_summary = f"""
Marker Genes and Cluster Similarity Analysis - Overall Summary
============================================================

Analysis completed for {len(samples)} samples:
{chr(10).join([f"- {sample}" for sample in samples])}

Overall Statistics:
- Total cells across all samples: {adata.n_obs:,}
- Total genes: {adata.n_vars:,}
- Samples processed: {len(samples)}

Per-Sample Results:
{chr(10).join([f"- {row['sample_id']}: {row['total_cells']} cells, {row['valid_clusters']} valid clusters (≥50 cells)" for _, row in overall_stats_df.iterrows()])}

Analysis Parameters:
- Differential expression method: Wilcoxon rank-sum test
- Number of marker genes per cluster: 50
- Distance metric for clustering: Correlation distance
- Linkage method: Ward
- Minimum cells per cluster threshold: 50

Note: Each sample was analyzed independently since cluster numbers 
are not correlated across samples (cluster 1 in sample A ≠ cluster 1 in sample B).
"""

with open(os.path.join(output_dir, 'overall_analysis_summary.txt'), 'w') as f:
    f.write(overall_summary)

print(overall_summary)
print(f"\nAnalysis completed successfully!")
print(f"Results saved to: {output_dir}")
print(f"Each sample has its own subdirectory with individual results.")

# %%