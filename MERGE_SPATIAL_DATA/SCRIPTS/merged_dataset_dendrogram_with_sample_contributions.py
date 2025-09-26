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
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Configure scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Define the path to the merged AnnData file
merged_h5ad_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'

print("="*80)
print("MERGED DATASET DENDROGRAM WITH SAMPLE CONTRIBUTIONS")
print("="*80)
print("This script creates a hierarchical clustering dendrogram of all clusters")
print("from the merged dataset and shows sample contribution percentages.")
print("="*80)

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
output_dir = 'MERGE_SPATIAL_DATA/merged_dataset_dendrogram_analysis'
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

print(f"Results will be saved to: {output_dir}")

# %%
# Analyze dataset structure
print(f"\n{'='*60}")
print("DATASET STRUCTURE ANALYSIS")
print(f"{'='*60}")

# Get unique samples and clusters
samples = adata.obs['sample_id'].unique()
clusters = adata.obs['leiden'].unique()

print(f"Found {len(samples)} samples: {list(samples)}")
print(f"Found {len(clusters)} clusters: {list(sorted(clusters, key=int))}")
print(f"Total cells: {adata.n_obs:,}")
print(f"Total genes: {adata.n_vars:,}")

# %%
# Clean up cluster categories to only include present clusters
adata.obs['leiden'] = adata.obs['leiden'].cat.remove_unused_categories()
valid_clusters = sorted(adata.obs['leiden'].cat.categories, key=int)

print(f"\nValid clusters for analysis: {valid_clusters}")

# %%
# Calculate sample contribution percentages for each cluster
print(f"\n{'='*60}")
print("CALCULATING SAMPLE CONTRIBUTIONS PER CLUSTER")
print(f"{'='*60}")

sample_contributions = {}
cluster_sample_counts = {}
cluster_total_counts = {}

for cluster in valid_clusters:
    cluster_mask = adata.obs['leiden'] == cluster
    cluster_cells = adata[cluster_mask]
    
    # Count cells from each sample in this cluster
    sample_counts = cluster_cells.obs['sample_id'].value_counts()
    total_cells_in_cluster = len(cluster_cells)
    
    # Calculate percentages
    sample_percentages = (sample_counts / total_cells_in_cluster * 100).round(2)
    
    sample_contributions[cluster] = sample_percentages.to_dict()
    cluster_sample_counts[cluster] = sample_counts.to_dict()
    cluster_total_counts[cluster] = total_cells_in_cluster
    
    print(f"Cluster {cluster}: {total_cells_in_cluster} cells")
    for sample_id, percentage in sample_percentages.items():
        count = sample_counts[sample_id]
        print(f"  {sample_id}: {count} cells ({percentage:.1f}%)")

# %%
# Create sample contributions DataFrame for easy export
contributions_data = []
for cluster in valid_clusters:
    for sample_id in samples:
        count = cluster_sample_counts[cluster].get(sample_id, 0)
        percentage = sample_contributions[cluster].get(sample_id, 0.0)
        contributions_data.append({
            'cluster': cluster,
            'sample_id': sample_id,
            'cell_count': count,
            'percentage': percentage,
            'total_cluster_cells': cluster_total_counts[cluster]
        })

contributions_df = pd.DataFrame(contributions_data)
contributions_df.to_csv(os.path.join(output_dir, 'cluster_sample_contributions.csv'), index=False)
print(f"\nSample contributions saved to: cluster_sample_contributions.csv")

# %%
# Calculate mean expression per cluster for hierarchical clustering
print(f"\n{'='*60}")
print("CALCULATING CLUSTER MEAN EXPRESSIONS")
print(f"{'='*60}")

cluster_means = []
cluster_labels = []

for cluster in valid_clusters:
    cluster_cells = adata[adata.obs['leiden'] == cluster]
    if len(cluster_cells) > 0:
        mean_expr = np.array(cluster_cells.X.mean(axis=0)).flatten()
        cluster_means.append(mean_expr)
        cluster_labels.append(f"Cluster_{cluster}\n({len(cluster_cells)} cells)")
        print(f"Processed Cluster {cluster}: {len(cluster_cells)} cells")

print(f"Calculated mean expressions for {len(cluster_means)} clusters")

# %%
# Perform hierarchical clustering
print(f"\n{'='*60}")
print("PERFORMING HIERARCHICAL CLUSTERING")
print(f"{'='*60}")

if len(cluster_means) >= 2:
    cluster_means_matrix = np.array(cluster_means)
    
    # Calculate pairwise distances between clusters
    print("Calculating pairwise distances using correlation metric...")
    distances = pdist(cluster_means_matrix, metric='correlation')
    
    print("Performing Ward linkage...")
    linkage_matrix = linkage(distances, method='ward')
    
    print("Hierarchical clustering completed successfully!")
else:
    print("Error: Need at least 2 clusters for hierarchical clustering")
    sys.exit()

# %%
# Create color palette for samples
print(f"\n{'='*60}")
print("PREPARING VISUALIZATION")
print(f"{'='*60}")

# Generate distinct colors for each sample
sample_colors = plt.cm.Set3(np.linspace(0, 1, len(samples)))
sample_color_map = dict(zip(samples, sample_colors))

print("Sample color mapping:")
for sample, color in sample_color_map.items():
    print(f"  {sample}: RGB{tuple(color[:3])}")

# %%
# Create the main dendrogram plot with sample contribution bars
print(f"\n{'='*60}")
print("CREATING DENDROGRAM WITH SAMPLE CONTRIBUTIONS")
print(f"{'='*60}")

# Create figure with specific layout for dendrogram and sample bars
fig = plt.figure(figsize=(20, 12))

# Create main dendrogram subplot
ax_dend = plt.subplot2grid((3, 4), (0, 0), colspan=3, rowspan=2)

# Create the dendrogram
dend_result = dendrogram(
    linkage_matrix,
    labels=cluster_labels,
    leaf_rotation=90,
    leaf_font_size=10,
    ax=ax_dend
)

ax_dend.set_title('Cluster Similarity Dendrogram - Full Merged Dataset\n(Based on Gene Expression Profiles)', 
                  fontsize=16, fontweight='bold')
ax_dend.set_xlabel('Clusters', fontsize=12)
ax_dend.set_ylabel('Distance (Ward Linkage)', fontsize=12)

# Get the order of clusters from the dendrogram
cluster_order = [label.replace('Cluster_', '').split('\n')[0] for label in dend_result['ivl']]

# %%
# Create sample contribution stacked bars next to dendrogram
ax_bars = plt.subplot2grid((3, 4), (0, 3), rowspan=2)

# Prepare data for stacked bars
bar_width = 0.8
x_positions = range(len(cluster_order))

# Create stacked bars showing sample contributions
bottom_values = np.zeros(len(cluster_order))

for i, sample_id in enumerate(samples):
    percentages = []
    for cluster in cluster_order:
        percentage = sample_contributions[cluster].get(sample_id, 0.0)
        percentages.append(percentage)
    
    ax_bars.bar(x_positions, percentages, bar_width, 
                bottom=bottom_values, label=sample_id, 
                color=sample_color_map[sample_id], alpha=0.8)
    
    bottom_values += np.array(percentages)

ax_bars.set_title('Sample Contributions\nper Cluster (%)', fontsize=12, fontweight='bold')
ax_bars.set_xlabel('Clusters', fontsize=10)
ax_bars.set_ylabel('Percentage (%)', fontsize=10)
ax_bars.set_xticks(x_positions)
ax_bars.set_xticklabels([f'C{c}' for c in cluster_order], rotation=90)
ax_bars.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
ax_bars.set_ylim(0, 100)

# Add grid for better readability
ax_bars.grid(True, alpha=0.3, axis='y')

# %%
# Create sample contribution heatmap
ax_heatmap = plt.subplot2grid((3, 4), (2, 0), colspan=4)

# Prepare heatmap data
heatmap_data = []
for cluster in cluster_order:
    row = []
    for sample_id in samples:
        percentage = sample_contributions[cluster].get(sample_id, 0.0)
        row.append(percentage)
    heatmap_data.append(row)

heatmap_df = pd.DataFrame(heatmap_data, 
                         index=[f'Cluster_{c}' for c in cluster_order], 
                         columns=samples)

# Create heatmap
sns.heatmap(heatmap_df, annot=True, fmt='.1f', cmap='YlOrRd', 
            cbar_kws={'label': 'Percentage (%)'}, ax=ax_heatmap,
            annot_kws={'size': 8})

ax_heatmap.set_title('Sample Contribution Heatmap per Cluster', fontsize=12, fontweight='bold')
ax_heatmap.set_xlabel('Samples', fontsize=10)
ax_heatmap.set_ylabel('Clusters (Dendrogram Order)', fontsize=10)

plt.tight_layout()

# Save the complete figure
output_path = os.path.join(output_dir, 'merged_dataset_dendrogram_with_sample_contributions.png')
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"Main combined plot saved to: {output_path}")

# %%
# Save individual subplot figures
print(f"\n{'='*60}")
print("SAVING INDIVIDUAL SUBPLOT FIGURES")
print(f"{'='*60}")

# 1. Save standalone dendrogram
fig_dend = plt.figure(figsize=(14, 8))
dend_result_standalone = dendrogram(
    linkage_matrix,
    labels=cluster_labels,
    leaf_rotation=90,
    leaf_font_size=12
)
plt.title('Cluster Similarity Dendrogram - Full Merged Dataset\n(Based on Gene Expression Profiles)',
          fontsize=16, fontweight='bold')
plt.xlabel('Clusters', fontsize=14)
plt.ylabel('Distance (Ward Linkage)', fontsize=14)
plt.tight_layout()

dendrogram_path = os.path.join(output_dir, 'standalone_dendrogram.png')
plt.savefig(dendrogram_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Standalone dendrogram saved to: {dendrogram_path}")

# 2. Save standalone sample contribution bars
fig_bars = plt.figure(figsize=(12, 8))
ax_bars_standalone = fig_bars.add_subplot(111)

# Recreate the stacked bars
cluster_order_standalone = [label.replace('Cluster_', '').split('\n')[0] for label in dend_result_standalone['ivl']]
x_positions = range(len(cluster_order_standalone))
bottom_values = np.zeros(len(cluster_order_standalone))

for i, sample_id in enumerate(samples):
    percentages = []
    for cluster in cluster_order_standalone:
        percentage = sample_contributions[cluster].get(sample_id, 0.0)
        percentages.append(percentage)
    
    ax_bars_standalone.bar(x_positions, percentages, 0.8,
                          bottom=bottom_values, label=sample_id,
                          color=sample_color_map[sample_id], alpha=0.8)
    
    bottom_values += np.array(percentages)

ax_bars_standalone.set_title('Sample Contributions per Cluster (%)', fontsize=16, fontweight='bold')
ax_bars_standalone.set_xlabel('Clusters (Dendrogram Order)', fontsize=14)
ax_bars_standalone.set_ylabel('Percentage (%)', fontsize=14)
ax_bars_standalone.set_xticks(x_positions)
ax_bars_standalone.set_xticklabels([f'Cluster {c}' for c in cluster_order_standalone], rotation=90)
ax_bars_standalone.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
ax_bars_standalone.set_ylim(0, 100)
ax_bars_standalone.grid(True, alpha=0.3, axis='y')
plt.tight_layout()

bars_path = os.path.join(output_dir, 'standalone_sample_contribution_bars.png')
plt.savefig(bars_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Standalone sample contribution bars saved to: {bars_path}")

# 3. Save standalone heatmap
fig_heatmap = plt.figure(figsize=(14, 10))
ax_heatmap_standalone = fig_heatmap.add_subplot(111)

# Recreate heatmap data in dendrogram order
heatmap_data_standalone = []
for cluster in cluster_order_standalone:
    row = []
    for sample_id in samples:
        percentage = sample_contributions[cluster].get(sample_id, 0.0)
        row.append(percentage)
    heatmap_data_standalone.append(row)

heatmap_df_standalone = pd.DataFrame(heatmap_data_standalone,
                                    index=[f'Cluster {c}' for c in cluster_order_standalone],
                                    columns=samples)

sns.heatmap(heatmap_df_standalone, annot=True, fmt='.1f', cmap='YlOrRd',
            cbar_kws={'label': 'Percentage (%)'}, ax=ax_heatmap_standalone,
            annot_kws={'size': 10})

ax_heatmap_standalone.set_title('Sample Contribution Heatmap per Cluster', fontsize=16, fontweight='bold')
ax_heatmap_standalone.set_xlabel('Samples', fontsize=14)
ax_heatmap_standalone.set_ylabel('Clusters (Dendrogram Order)', fontsize=14)
ax_heatmap_standalone.tick_params(axis='x', labelrotation=90)
plt.tight_layout()

heatmap_path = os.path.join(output_dir, 'standalone_sample_contribution_heatmap.png')
plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Standalone sample contribution heatmap saved to: {heatmap_path}")

# 4. Save UMAP plots colored by clusters and samples
fig_umap = plt.figure(figsize=(16, 6))

# Get UMAP coordinates
umap_coords = adata.obsm['X_umap']

# Left subplot: UMAP colored by clusters
ax1 = fig_umap.add_subplot(121)
for i, cluster in enumerate(valid_clusters):
    cluster_mask = adata.obs['leiden'] == cluster
    if cluster_mask.sum() > 0:
        # Use consistent colors with the rest of the analysis
        color = plt.cm.tab20(int(cluster) % 20)  # Cycle through colors for many clusters
        ax1.scatter(umap_coords[cluster_mask, 0], umap_coords[cluster_mask, 1],
                   c=[color], s=0.5, alpha=0.7, label=f'Cluster {cluster}')

ax1.set_xlabel('UMAP1', fontsize=12)
ax1.set_ylabel('UMAP2', fontsize=12)
ax1.set_title('UMAP Colored by Leiden Clusters', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)

# Create legend for clusters with smaller text
handles, labels = ax1.get_legend_handles_labels()
# Split legend into two columns if too many clusters
if len(handles) > 12:
    ax1.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left',
              fontsize=12, ncol=2, markerscale=12)
else:
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12, markerscale=12)

# Right subplot: UMAP colored by samples
ax2 = fig_umap.add_subplot(122)
for i, sample_id in enumerate(samples):
    sample_mask = adata.obs['sample_id'] == sample_id
    if sample_mask.sum() > 0:
        color = sample_color_map[sample_id]
        ax2.scatter(umap_coords[sample_mask, 0], umap_coords[sample_mask, 1],
                   c=[color], s=0.5, alpha=0.7, label=sample_id)

ax2.set_xlabel('UMAP1', fontsize=12)
ax2.set_ylabel('UMAP2', fontsize=12)
ax2.set_title('UMAP Colored by Sample ID', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12, markerscale=12)

plt.tight_layout()

umap_path = os.path.join(output_dir, 'umap_clusters_and_samples.png')
plt.savefig(umap_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"UMAP plots (clusters and samples) saved to: {umap_path}")

# 5. Save individual UMAP plots
# Individual UMAP by clusters
fig_umap_clusters = plt.figure(figsize=(10, 8))
ax_clusters = fig_umap_clusters.add_subplot(111)

for i, cluster in enumerate(valid_clusters):
    cluster_mask = adata.obs['leiden'] == cluster
    if cluster_mask.sum() > 0:
        color = plt.cm.tab20(int(cluster) % 20)
        ax_clusters.scatter(umap_coords[cluster_mask, 0], umap_coords[cluster_mask, 1],
                           c=[color], s=1, alpha=0.7, label=f'Cluster {cluster}')

ax_clusters.set_xlabel('UMAP1', fontsize=14)
ax_clusters.set_ylabel('UMAP2', fontsize=14)
ax_clusters.set_title('UMAP Colored by Leiden Clusters', fontsize=16, fontweight='bold')
ax_clusters.grid(True, alpha=0.3)

# Create legend
handles, labels = ax_clusters.get_legend_handles_labels()
if len(handles) > 12:
    ax_clusters.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left',
                      fontsize=10, ncol=2, markerscale=8)
else:
    ax_clusters.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, markerscale=8)

plt.tight_layout()

umap_clusters_path = os.path.join(output_dir, 'standalone_umap_by_clusters.png')
plt.savefig(umap_clusters_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Standalone UMAP by clusters saved to: {umap_clusters_path}")

# Individual UMAP by samples
fig_umap_samples = plt.figure(figsize=(10, 8))
ax_samples = fig_umap_samples.add_subplot(111)

for i, sample_id in enumerate(samples):
    sample_mask = adata.obs['sample_id'] == sample_id
    if sample_mask.sum() > 0:
        color = sample_color_map[sample_id]
        ax_samples.scatter(umap_coords[sample_mask, 0], umap_coords[sample_mask, 1],
                          c=[color], s=1, alpha=0.7, label=sample_id)

ax_samples.set_xlabel('UMAP1', fontsize=14)
ax_samples.set_ylabel('UMAP2', fontsize=14)
ax_samples.set_title('UMAP Colored by Sample ID', fontsize=16, fontweight='bold')
ax_samples.grid(True, alpha=0.3)
ax_samples.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, markerscale=8)

plt.tight_layout()

umap_samples_path = os.path.join(output_dir, 'standalone_umap_by_samples.png')
plt.savefig(umap_samples_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Standalone UMAP by samples saved to: {umap_samples_path}")

# %%
# Create separate detailed sample contribution analysis
print(f"\n{'='*60}")
print("CREATING DETAILED SAMPLE CONTRIBUTION ANALYSIS")
print(f"{'='*60}")

# Create summary statistics
summary_stats = {}
for sample_id in samples:
    sample_stats = {
        'total_cells': len(adata[adata.obs['sample_id'] == sample_id]),
        'clusters_present': len([c for c in valid_clusters if sample_contributions[c].get(sample_id, 0) > 0]),
        'avg_contribution_when_present': np.mean([sample_contributions[c].get(sample_id, 0) 
                                                 for c in valid_clusters 
                                                 if sample_contributions[c].get(sample_id, 0) > 0])
    }
    summary_stats[sample_id] = sample_stats

# Create sample summary DataFrame
sample_summary_data = []
for sample_id, stats in summary_stats.items():
    sample_summary_data.append({
        'sample_id': sample_id,
        'total_cells': stats['total_cells'],
        'clusters_present_in': stats['clusters_present'],
        'avg_contribution_when_present': round(stats['avg_contribution_when_present'], 2)
    })

sample_summary_df = pd.DataFrame(sample_summary_data)
sample_summary_df.to_csv(os.path.join(output_dir, 'sample_summary_statistics.csv'), index=False)

# %%
# Create cluster summary statistics
cluster_summary_data = []
for cluster in valid_clusters:
    cluster_info = {
        'cluster': cluster,
        'total_cells': cluster_total_counts[cluster],
        'num_samples_contributing': len([s for s in samples if sample_contributions[cluster].get(s, 0) > 0]),
        'dominant_sample': max(sample_contributions[cluster].items(), key=lambda x: x[1])[0],
        'dominant_sample_percentage': max(sample_contributions[cluster].values()),
        'sample_diversity_index': len([s for s in samples if sample_contributions[cluster].get(s, 0) > 5.0])  # samples with >5% contribution
    }
    cluster_summary_data.append(cluster_info)

cluster_summary_df = pd.DataFrame(cluster_summary_data)
cluster_summary_df.to_csv(os.path.join(output_dir, 'cluster_summary_statistics.csv'), index=False)

# %%
# Save linkage matrix and dendrogram order for reproducibility
linkage_df = pd.DataFrame(linkage_matrix, columns=['cluster1', 'cluster2', 'distance', 'sample_count'])
linkage_df.to_csv(os.path.join(output_dir, 'hierarchical_clustering_linkage_matrix.csv'), index=False)

dendrogram_order_df = pd.DataFrame({'position': range(len(cluster_order)), 'cluster': cluster_order})
dendrogram_order_df.to_csv(os.path.join(output_dir, 'dendrogram_cluster_order.csv'), index=False)

# %%
# Generate comprehensive analysis report
print(f"\n{'='*60}")
print("GENERATING ANALYSIS REPORT")
print(f"{'='*60}")

report_content = f"""MERGED DATASET DENDROGRAM WITH SAMPLE CONTRIBUTIONS - ANALYSIS REPORT
===============================================================================

DATASET OVERVIEW:
- Total cells analyzed: {adata.n_obs:,}
- Total genes: {adata.n_vars:,}
- Number of samples: {len(samples)}
- Number of clusters: {len(valid_clusters)}
- Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

SAMPLES INCLUDED:
{chr(10).join([f"- {sample}: {len(adata[adata.obs['sample_id'] == sample])} cells" for sample in samples])}

CLUSTERS ANALYZED:
{chr(10).join([f"- Cluster {cluster}: {cluster_total_counts[cluster]} cells" for cluster in valid_clusters])}

HIERARCHICAL CLUSTERING PARAMETERS:
- Distance metric: Correlation distance
- Linkage method: Ward linkage
- Number of clusters in dendrogram: {len(valid_clusters)}

SAMPLE CONTRIBUTION SUMMARY:
{chr(10).join([f"- {row['sample_id']}: {row['total_cells']} total cells, present in {row['clusters_present_in']} clusters, avg contribution {row['avg_contribution_when_present']:.1f}%" for _, row in sample_summary_df.iterrows()])}

CLUSTER DIVERSITY ANALYSIS:
{chr(10).join([f"- Cluster {row['cluster']}: {row['total_cells']} cells, {row['num_samples_contributing']} samples contributing, dominated by {row['dominant_sample']} ({row['dominant_sample_percentage']:.1f}%)" for _, row in cluster_summary_df.iterrows()])}

FILES GENERATED:
1. merged_dataset_dendrogram_with_sample_contributions.png - Main combined visualization
2. standalone_dendrogram.png - Individual dendrogram plot
3. standalone_sample_contribution_bars.png - Individual sample contribution bars
4. standalone_sample_contribution_heatmap.png - Individual sample contribution heatmap
5. umap_clusters_and_samples.png - UMAP plots (clusters and samples side-by-side)
6. standalone_umap_by_clusters.png - Individual UMAP colored by clusters
7. standalone_umap_by_samples.png - Individual UMAP colored by samples
8. cluster_sample_contributions.csv - Detailed contribution data per cluster
9. sample_summary_statistics.csv - Summary statistics per sample
10. cluster_summary_statistics.csv - Summary statistics per cluster
11. hierarchical_clustering_linkage_matrix.csv - Linkage matrix for reproducibility
12. dendrogram_cluster_order.csv - Cluster order in dendrogram
13. analysis_report.txt - This comprehensive report

ANALYSIS METHODOLOGY:
1. Loaded merged dataset with all samples and clusters
2. Calculated mean gene expression per cluster across entire dataset
3. Performed hierarchical clustering using correlation distance and Ward linkage
4. Computed sample contribution percentages for each cluster
5. Created dendrogram visualization with color-coded sample contribution bars
6. Generated sample contribution heatmap showing cluster composition
7. Produced summary statistics and detailed data exports

INTERPRETATION GUIDELINES:
- The dendrogram shows similarity relationships between clusters based on gene expression
- Sample contribution bars show what percentage of each cluster comes from each sample
- Clusters with high diversity (multiple samples contributing) may represent common cell types
- Clusters dominated by single samples may represent sample-specific populations
- The correlation distance metric emphasizes gene expression pattern similarity
"""

# Save the report
report_path = os.path.join(output_dir, 'analysis_report.txt')
with open(report_path, 'w') as f:
    f.write(report_content)

print(f"Analysis report saved to: {report_path}")

# %%
# Print final summary
print(f"\n{'='*80}")
print("ANALYSIS COMPLETED SUCCESSFULLY!")
print(f"{'='*80}")
print(f"Results saved to: {output_dir}")
print(f"\nKey outputs:")
print(f"- Main combined visualization: merged_dataset_dendrogram_with_sample_contributions.png")
print(f"- Individual dendrogram: standalone_dendrogram.png")
print(f"- Individual sample bars: standalone_sample_contribution_bars.png")
print(f"- Individual heatmap: standalone_sample_contribution_heatmap.png")
print(f"- UMAP plots (combined): umap_clusters_and_samples.png")
print(f"- UMAP by clusters: standalone_umap_by_clusters.png")
print(f"- UMAP by samples: standalone_umap_by_samples.png")
print(f"- Contribution data: cluster_sample_contributions.csv")
print(f"- Analysis report: analysis_report.txt")
print(f"\nThe dendrogram shows hierarchical relationships between all {len(valid_clusters)} clusters")
print(f"from the merged dataset, with sample contribution percentages visualized")
print(f"as color-coded stacked bars.")
print(f"{'='*80}")