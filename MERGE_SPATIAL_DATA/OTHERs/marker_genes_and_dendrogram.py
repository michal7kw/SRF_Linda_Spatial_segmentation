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
output_dir = 'MERGE_SPATIAL_DATA/marker_genes_analysis'
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

print(f"Results will be saved to: {output_dir}")

# %%
# Analyze cluster composition per sample and filter clusters with ≥50 cells
print("Analyzing cluster composition per sample...")

cluster_stats = []
samples = adata.obs['sample_id'].unique()

for sample in samples:
    sample_data = adata[adata.obs['sample_id'] == sample]
    cluster_counts = sample_data.obs['leiden'].value_counts()
    
    for cluster, count in cluster_counts.items():
        cluster_stats.append({
            'sample_id': sample,
            'cluster': cluster,
            'cell_count': count,
            'meets_threshold': count >= 50
        })

cluster_stats_df = pd.DataFrame(cluster_stats)
print(f"\nCluster statistics:")
print(cluster_stats_df.head(10))

# Save cluster statistics
cluster_stats_df.to_csv(os.path.join(output_dir, 'cluster_statistics.csv'), index=False)

# %%
# Filter to only include clusters with ≥50 cells per sample
valid_clusters = cluster_stats_df[cluster_stats_df['meets_threshold']].copy()
print(f"\nClusters with ≥50 cells:")
print(f"Total valid cluster-sample combinations: {len(valid_clusters)}")
print(valid_clusters.groupby('sample_id')['cluster'].count())

# %%
# Prepare data for marker gene analysis
print("\nPreparing data for marker gene analysis...")

# Create a copy of the data for analysis
adata_analysis = adata.copy()

# Ensure leiden clusters are categorical
adata_analysis.obs['leiden'] = adata_analysis.obs['leiden'].astype('category')

print(f"Total clusters in dataset: {len(adata_analysis.obs['leiden'].cat.categories)}")
print(f"Clusters: {sorted(adata_analysis.obs['leiden'].cat.categories.astype(int))}")

# %%
# Perform differential expression analysis for marker genes
print("\nPerforming differential expression analysis...")

# Use Wilcoxon rank-sum test for marker gene identification
sc.tl.rank_genes_groups(
    adata_analysis, 
    'leiden', 
    method='wilcoxon',
    key_added='rank_genes_leiden',
    n_genes=50  # Top 50 genes per cluster
)

print("Differential expression analysis completed.")

# %%
# Extract and save marker gene results
print("Extracting marker gene results...")

# Get the results for all groups
marker_genes_list = []
for group in adata_analysis.obs['leiden'].cat.categories:
    group_df = sc.get.rank_genes_groups_df(adata_analysis, group=group, key='rank_genes_leiden')
    group_df['group'] = group  # Add group column
    marker_genes_list.append(group_df)

marker_genes_df = pd.concat(marker_genes_list, ignore_index=True)
print(f"Marker genes dataframe shape: {marker_genes_df.shape}")
print(f"Columns: {marker_genes_df.columns.tolist()}")
print(marker_genes_df.head())

# Save all marker genes
marker_genes_df.to_csv(os.path.join(output_dir, 'all_marker_genes.csv'), index=False)

# %%
# Create summary of top marker genes per cluster
print("Creating summary of top marker genes...")

top_markers_per_cluster = {}
n_top_genes = 10

for cluster in adata_analysis.obs['leiden'].cat.categories:
    cluster_markers = marker_genes_df[marker_genes_df['group'] == cluster].head(n_top_genes)
    top_markers_per_cluster[cluster] = cluster_markers['names'].tolist()

# Save top markers summary
top_markers_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in top_markers_per_cluster.items()]))
top_markers_df.to_csv(os.path.join(output_dir, f'top_{n_top_genes}_marker_genes_per_cluster.csv'), index=False)

print(f"Top {n_top_genes} marker genes per cluster:")
print(top_markers_df.head())

# %%
# Plot marker genes heatmap
print("Creating marker genes heatmap...")

sc.pl.rank_genes_groups_heatmap(
    adata_analysis,
    key='rank_genes_leiden',
    n_genes=5,  # Top 5 genes per cluster for visualization
    show_gene_labels=True,
    figsize=(12, 8),
    save='_marker_genes_heatmap.png',
    show=False
)

# %%
# Create dotplot for marker genes
print("Creating marker genes dotplot...")

sc.pl.rank_genes_groups_dotplot(
    adata_analysis,
    key='rank_genes_leiden',
    n_genes=5,  # Top 5 genes per cluster
    figsize=(15, 8),
    save='_marker_genes_dotplot.png',
    show=False
)

# %%
# Prepare data for dendrogram analysis
print("\nPreparing data for cluster similarity dendrogram...")

# Calculate mean expression per cluster
cluster_means = []
cluster_labels = []

for cluster in sorted(adata_analysis.obs['leiden'].cat.categories, key=int):
    cluster_cells = adata_analysis[adata_analysis.obs['leiden'] == cluster]
    if len(cluster_cells) > 0:  # Only include clusters with cells
        mean_expr = np.array(cluster_cells.X.mean(axis=0)).flatten()
        cluster_means.append(mean_expr)
        cluster_labels.append(f"Cluster_{cluster}")

cluster_means_matrix = np.array(cluster_means)
print(f"Cluster means matrix shape: {cluster_means_matrix.shape}")

# %%
# Calculate pairwise distances between clusters
print("Calculating pairwise distances between clusters...")

# Use correlation distance (1 - correlation) for similarity
distances = pdist(cluster_means_matrix, metric='correlation')

# Perform hierarchical clustering
linkage_matrix = linkage(distances, method='ward')

# %%
# Create and save dendrogram
print("Creating cluster similarity dendrogram...")

plt.figure(figsize=(12, 8))
dendrogram_plot = dendrogram(
    linkage_matrix,
    labels=cluster_labels,
    leaf_rotation=45,
    leaf_font_size=10
)

plt.title('Cluster Similarity Dendrogram\n(Based on Gene Expression Profiles)', fontsize=14)
plt.xlabel('Clusters', fontsize=12)
plt.ylabel('Distance (Ward Linkage)', fontsize=12)
plt.tight_layout()

# Save the dendrogram
dendrogram_path = os.path.join(output_dir, 'cluster_similarity_dendrogram.png')
plt.savefig(dendrogram_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"Dendrogram saved to: {dendrogram_path}")

# %%
# Create a more detailed dendrogram with sample information
print("Creating detailed dendrogram with cluster statistics...")

# Calculate cluster statistics for annotation
cluster_info = []
for cluster in sorted(adata_analysis.obs['leiden'].cat.categories, key=int):
    cluster_data = adata_analysis[adata_analysis.obs['leiden'] == cluster]
    total_cells = len(cluster_data)
    
    # Count cells per sample
    sample_counts = cluster_data.obs['sample_id'].value_counts()
    
    cluster_info.append({
        'cluster': cluster,
        'total_cells': total_cells,
        'samples': dict(sample_counts)
    })

cluster_info_df = pd.DataFrame(cluster_info)
cluster_info_df.to_csv(os.path.join(output_dir, 'cluster_detailed_info.csv'), index=False)

# %%
# Create correlation heatmap between clusters
print("Creating cluster correlation heatmap...")

# Calculate correlation matrix
correlation_matrix = np.corrcoef(cluster_means_matrix)

plt.figure(figsize=(10, 8))
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

plt.title('Cluster Correlation Heatmap\n(Gene Expression Similarity)', fontsize=14)
plt.xticks(rotation=45)
plt.yticks(rotation=0)
plt.tight_layout()

correlation_heatmap_path = os.path.join(output_dir, 'cluster_correlation_heatmap.png')
plt.savefig(correlation_heatmap_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"Correlation heatmap saved to: {correlation_heatmap_path}")

# %%
# Save correlation matrix as CSV
correlation_df = pd.DataFrame(
    correlation_matrix,
    index=cluster_labels,
    columns=cluster_labels
)
correlation_df.to_csv(os.path.join(output_dir, 'cluster_correlation_matrix.csv'))

# %%
# Create summary report
print("Creating summary report...")

summary_report = f"""
Marker Genes and Cluster Similarity Analysis Report
==================================================

Dataset Information:
- Total cells: {adata.n_obs:,}
- Total genes: {adata.n_vars:,}
- Total samples: {len(adata.obs['sample_id'].unique())}
- Total clusters: {len(adata.obs['leiden'].cat.categories)}

Cluster Analysis:
- Clusters with ≥50 cells (any sample): {len(valid_clusters['cluster'].unique())}
- Total valid cluster-sample combinations: {len(valid_clusters)}

Files Generated:
1. cluster_statistics.csv - Cell counts per cluster per sample
2. all_marker_genes.csv - Complete marker gene results
3. top_{n_top_genes}_marker_genes_per_cluster.csv - Top marker genes summary
4. cluster_detailed_info.csv - Detailed cluster information
5. cluster_correlation_matrix.csv - Cluster correlation matrix
6. marker_genes_heatmap.png - Marker genes heatmap
7. marker_genes_dotplot.png - Marker genes dotplot
8. cluster_similarity_dendrogram.png - Hierarchical clustering dendrogram
9. cluster_correlation_heatmap.png - Cluster correlation heatmap

Analysis Parameters:
- Differential expression method: Wilcoxon rank-sum test
- Number of marker genes per cluster: 50
- Distance metric for clustering: Correlation distance
- Linkage method: Ward
- Minimum cells per cluster threshold: 50
"""

with open(os.path.join(output_dir, 'analysis_summary.txt'), 'w') as f:
    f.write(summary_report)

print(summary_report)
print(f"\nAnalysis completed successfully!")
print(f"All results saved to: {output_dir}")

# %%