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
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Define the path to the merged AnnData file
merged_h5ad_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'

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
# Create output directory for combined results
output_dir = 'MERGE_SPATIAL_DATA/correlation_analysis_cluster0_recluster'
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

print(f"Combined results will be saved to: {output_dir}")

# %%
# Function to perform reclustering on cluster 0 cells for the entire dataset
def recluster_all_cluster0_cells(adata, cluster_col='leiden_new', cluster_to_recluster='0'):
    """
    Recluster all cluster 0 cells across all samples
    """
    print("Performing combined cluster 0 reclustering across all samples...")
    
    # Extract all cluster 0 cells from all samples
    cluster0_cells = adata[adata.obs[cluster_col] == cluster_to_recluster].copy()
    
    if len(cluster0_cells) < 10:
        print(f"  Insufficient cluster 0 cells for reclustering: {len(cluster0_cells)} cells")
        return adata, {}
    
    print(f"  Reclustering {len(cluster0_cells)} cluster 0 cells from all samples")
    
    # Perform neighborhood graph and clustering on all cluster 0 cells
    sc.pp.neighbors(cluster0_cells, n_neighbors=min(30, len(cluster0_cells)-1), n_pcs=40)
    sc.tl.leiden(cluster0_cells, resolution=0.8, key_added='leiden_combined_recluster')
    
    # Generate UMAP embedding for all cluster 0 cells
    sc.tl.umap(cluster0_cells)
    
    # Create new cluster labels for the full dataset
    new_cluster_labels = adata.obs[cluster_col].astype(str).copy()
    
    # Map reclustered cells to new labels
    recluster_map = {}
    for idx, cell in enumerate(cluster0_cells.obs_names):
        old_label = cluster0_cells.obs.loc[cell, 'leiden_combined_recluster']
        new_label = f"0.{old_label}"
        new_cluster_labels.loc[cell] = new_label
        if new_label not in recluster_map:
            recluster_map[new_label] = {'count': 0, 'samples': set()}
        recluster_map[new_label]['count'] += 1
        recluster_map[new_label]['samples'].add(cluster0_cells.obs.loc[cell, 'sample_id'])
    
    # Add new cluster column to the full dataset
    adata.obs['leiden_combined_cluster0_reclustered'] = new_cluster_labels
    adata.obs['leiden_combined_cluster0_reclustered'] = adata.obs['leiden_combined_cluster0_reclustered'].astype('category')
    
    return adata, {'cluster0_subset': cluster0_cells, 'recluster_map': recluster_map}

# %%
# Perform combined cluster 0 reclustering
adata_with_combined_recluster, combined_recluster_results = recluster_all_cluster0_cells(adata)

# %%
# Use the reclustered labels
cluster_column = 'leiden_combined_cluster0_reclustered'
print(f"Using combined reclustered labels: {cluster_column}")

# Print reclustering summary
print("\nCombined Cluster 0 Reclustering Summary:")
for subcluster, info in combined_recluster_results['recluster_map'].items():
    samples_list = sorted(list(info['samples']))
    print(f"  {subcluster}: {info['count']} cells from {len(info['samples'])} samples ({', '.join(samples_list)})")

# Update cluster counts
cluster_counts = adata_with_combined_recluster.obs[cluster_column].value_counts()
valid_clusters_for_de = cluster_counts[cluster_counts > 1].index.tolist()

print(f"\nCombined dataset analysis:")
print(f"  Total cells: {adata_with_combined_recluster.n_obs:,}")
print(f"  Total clusters: {len(cluster_counts)}")
print(f"  Clusters with >1 cell for DE analysis: {len(valid_clusters_for_de)}")

# %%
# Create combined UMAP plots
print("Creating combined UMAP visualizations...")

# Combined UMAP plots - only the two UMAP plots side by side
fig, axes = plt.subplots(1, 2, figsize=(20, 8))

# UMAP colored by final clusters
sc.pl.umap(adata_with_combined_recluster, color=cluster_column, ax=axes[0], show=False, 
          legend_loc='right margin', legend_fontsize=8, s=15)
axes[0].set_title('Combined UMAP - All Samples\n(With Cluster 0 Reclustered)', fontsize=14)

# UMAP colored by sample
sc.pl.umap(adata_with_combined_recluster, color='sample_id', ax=axes[1], show=False,
          legend_loc='right margin', legend_fontsize=8, s=15)
axes[1].set_title('Combined UMAP - Colored by Sample', fontsize=14)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'UMAP_full.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

# %%
# Create cluster 0 subset combined analysis
# print("Creating combined cluster 0 subset analysis...")

# cluster0_combined = combined_recluster_results['cluster0_subset']

# # Combined cluster 0 subset plots
# fig, axes = plt.subplots(2, 1, figsize=(20, 16))

# # UMAP of cluster 0 subset with combined reclusters
# sc.pl.umap(cluster0_combined, color='leiden_combined_recluster', ax=axes[0], show=False,
#           legend_loc='on data', legend_fontsize=8, s=30)
# axes[0].set_title('Cluster 0 Combined UMAP\n(All Samples Reclustered Together)', fontsize=14)

# # UMAP of cluster 0 subset colored by sample
# sc.pl.umap(cluster0_combined, color='sample_id', ax=axes[1], show=False,
#           legend_loc='right margin', legend_fontsize=8, s=30)
# axes[1].set_title('Cluster 0 Combined UMAP\n(Colored by Sample)', fontsize=14)

# plt.tight_layout()
# plt.savefig(os.path.join(output_dir, 'UMAP_cluster0.png'), 
#             dpi=300, bbox_inches='tight')
# plt.close()

# %%
# Perform combined differential expression analysis
if len(valid_clusters_for_de) >= 2:
    print("Performing combined differential expression analysis...")
    
    sc.tl.rank_genes_groups(
        adata_with_combined_recluster, 
        cluster_column, 
        groups=valid_clusters_for_de,
        method='wilcoxon',
        key_added=f'rank_genes_combined_{cluster_column}',
        n_genes=100
    )
    
    # Extract combined marker gene results
    print("Extracting combined marker gene results...")
    
    combined_marker_genes_list = []
    for group in valid_clusters_for_de: 
        group_df = sc.get.rank_genes_groups_df(adata_with_combined_recluster, 
                                             group=group, 
                                             key=f'rank_genes_combined_{cluster_column}')
        group_df['group'] = group
        group_df['analysis_type'] = 'combined'
        combined_marker_genes_list.append(group_df)
    
    if combined_marker_genes_list:
        combined_marker_genes_df = pd.concat(combined_marker_genes_list, ignore_index=True)
        
        # Save combined marker genes
        marker_genes_filename = 'marker_genes.csv'
        combined_marker_genes_df.to_csv(os.path.join(output_dir, marker_genes_filename), index=False)
        
        # Create summary of top marker genes per cluster
        print("Creating combined top marker genes summary...")
        
        n_top_genes = 50
        top_markers_per_cluster = {}
        for cluster in valid_clusters_for_de:
            cluster_markers = combined_marker_genes_df[combined_marker_genes_df['group'] == cluster].head(n_top_genes)
            top_markers_per_cluster[cluster] = cluster_markers['names'].tolist()
        
        # Save top markers summary
        if top_markers_per_cluster:
            top_markers_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in top_markers_per_cluster.items()]))
            top_markers_filename = f'top_{n_top_genes}_marker_genes.csv'
            top_markers_df.to_csv(os.path.join(output_dir, top_markers_filename), index=False)
    
    # %%
    # Create combined visualizations
    print("Creating combined visualizations...")
    
    try:
        # Combined marker genes heatmap
        print("Creating combined heatmap...")
        heatmap_ax_dict = sc.pl.rank_genes_groups_heatmap(
            adata_with_combined_recluster,
            key=f'rank_genes_combined_{cluster_column}',
            n_genes=5,
            show_gene_labels=True,
            figsize=(16, 12),
            show=False,
        )
        ax = heatmap_ax_dict['heatmap_ax']
        ax.tick_params(axis='x', labelsize=6, rotation=90)
        ax.tick_params(axis='y', labelsize=8)
        fig = ax.get_figure()
        heatmap_filename = 'heatmap_marker_genes.png'
        fig.savefig(os.path.join(output_dir, heatmap_filename), dpi=300, bbox_inches='tight')
        plt.close('all')
        print(f"Combined heatmap saved: {heatmap_filename}")
        
    except Exception as e:
        print(f"Warning: Could not create combined heatmap: {e}")
    
    try:
        # Combined dotplot
        print("Creating combined dotplot...")
        dotplot_ax_dict = sc.pl.rank_genes_groups_dotplot(
            adata_with_combined_recluster,
            key=f'rank_genes_combined_{cluster_column}',
            n_genes=5,
            figsize=(18, 12),
            show=False,
        )
        
        ax = None
        if dotplot_ax_dict:
            if 'mainplot_ax' in dotplot_ax_dict and dotplot_ax_dict['mainplot_ax'] is not None:
                ax = dotplot_ax_dict['mainplot_ax']
            elif 'ax' in dotplot_ax_dict and dotplot_ax_dict['ax'] is not None:
                ax = dotplot_ax_dict['ax']
            else:
                for key, value in dotplot_ax_dict.items():
                    if hasattr(value, 'tick_params') and hasattr(value, 'get_figure'):
                        ax = value
                        break
        
        if ax is not None:
            ax.tick_params(axis='x', labelsize=6, rotation=90)
            ax.tick_params(axis='y', labelsize=8)
            fig = ax.get_figure()
            dotplot_filename = 'dotplot_marker_genes.png'
            fig.savefig(os.path.join(output_dir, dotplot_filename), dpi=300, bbox_inches='tight')
            plt.close('all')
            print(f"Combined dotplot saved: {dotplot_filename}")
        
    except Exception as e:
        print(f"Warning: Could not create combined dotplot: {e}")
    
    # %%
    # Create combined dendrogram and correlation analysis
    print("Creating combined dendrogram...")
    
    # Calculate mean expression per cluster across all samples
    cluster_means = []
    cluster_labels = []
    
    for cluster in sorted(valid_clusters_for_de, key=lambda x: (len(x.split('.')), x)):
        cluster_cells = adata_with_combined_recluster[adata_with_combined_recluster.obs[cluster_column] == cluster]
        if len(cluster_cells) > 0:
            mean_expr = np.array(cluster_cells.X.mean(axis=0)).flatten()
            cluster_means.append(mean_expr)
            cluster_labels.append(f"Cluster_{cluster}")
    
    if len(cluster_means) >= 2:
        cluster_means_matrix = np.array(cluster_means)
        
        # Calculate pairwise distances between clusters
        distances = pdist(cluster_means_matrix, metric='correlation')
        linkage_matrix = linkage(distances, method='ward')
        
        # Create combined dendrogram
        plt.figure(figsize=(16, 10))
        dendrogram(
            linkage_matrix,
            labels=cluster_labels,
            leaf_rotation=45,
            leaf_font_size=8
        )
        
        plt.title('Combined Cluster Similarity Dendrogram - All Samples (with Combined Cluster 0 Reclustering)\n(Based on Gene Expression Profiles)', fontsize=12)
        plt.xlabel('Clusters', fontsize=10)
        plt.ylabel('Distance (Ward Linkage)', fontsize=10)
        plt.tight_layout()
        
        # Save the combined dendrogram
        dendrogram_filename = 'dendrogram.png'
        plt.savefig(os.path.join(output_dir, dendrogram_filename), dpi=300, bbox_inches='tight')
        plt.close()
        
        # Create combined correlation heatmap
        correlation_matrix = np.corrcoef(cluster_means_matrix)
        
        plt.figure(figsize=(16, 12))
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
        
        plt.title('Combined Cluster Correlation Heatmap - All Samples (with Combined Cluster 0 Reclustering)', fontsize=12)
        plt.xticks(rotation=45, fontsize=8)
        plt.yticks(rotation=0, fontsize=8)
        plt.tight_layout()
        
        correlation_filename = 'clusters_correlation.png'
        plt.savefig(os.path.join(output_dir, correlation_filename), dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save combined correlation matrix
        correlation_df = pd.DataFrame(
            correlation_matrix,
            index=cluster_labels,
            columns=cluster_labels
        )
        correlation_csv_filename = 'clusters_correlation.csv'
        correlation_df.to_csv(os.path.join(output_dir, correlation_csv_filename))
        
        print(f"Combined dendrogram and correlation analysis completed")
    
    else:
        print("Skipping combined dendrogram: insufficient clusters for correlation analysis")

# %%
# Create sample composition analysis
print("Creating sample composition analysis...")

# Analyze cluster composition across samples
composition_data = []
for cluster in sorted(adata_with_combined_recluster.obs[cluster_column].cat.categories, 
                     key=lambda x: (len(x.split('.')), x)):
    cluster_cells = adata_with_combined_recluster[adata_with_combined_recluster.obs[cluster_column] == cluster]
    sample_counts = cluster_cells.obs['sample_id'].value_counts()
    
    for sample in adata_with_combined_recluster.obs['sample_id'].cat.categories:
        count = sample_counts.get(sample, 0)
        total_in_sample = len(adata_with_combined_recluster[adata_with_combined_recluster.obs['sample_id'] == sample])
        percentage = (count / total_in_sample * 100) if total_in_sample > 0 else 0
        
        composition_data.append({
            'cluster': cluster,
            'sample': sample,
            'cell_count': count,
            'percentage': percentage,
            'total_cells_in_sample': total_in_sample
        })

composition_df = pd.DataFrame(composition_data)
composition_df.to_csv(os.path.join(output_dir, 'cluster_composition_across_samples.csv'), index=False)

# Create composition heatmap
pivot_counts = composition_df.pivot(index='cluster', columns='sample', values='cell_count').fillna(0)
pivot_percentages = composition_df.pivot(index='cluster', columns='sample', values='percentage').fillna(0)

# Cell count heatmap
plt.figure(figsize=(14, 20))
sns.heatmap(pivot_counts, annot=True, fmt='g', cmap='viridis', cbar_kws={'label': 'Cell Count'})
plt.title('Cluster Cell Count Distribution Across Samples', fontsize=14)
plt.xlabel('Sample', fontsize=12)
plt.ylabel('Cluster', fontsize=12)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'cluster_cell_count_heatmap.png'), dpi=300, bbox_inches='tight')
plt.close()

# Percentage heatmap
plt.figure(figsize=(14, 20))
sns.heatmap(pivot_percentages, annot=True, fmt='.1f', cmap='plasma', cbar_kws={'label': 'Percentage (%)'})
plt.title('Cluster Percentage Distribution Across Samples', fontsize=14)
plt.xlabel('Sample', fontsize=12)
plt.ylabel('Cluster', fontsize=12)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'cluster_percentage_heatmap.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Sample composition analysis completed")

# %%
# Create comprehensive summary report
print("Creating comprehensive summary report...")

# Collect detailed statistics
samples = adata_with_combined_recluster.obs['sample_id'].unique()
total_cluster0_cells = len(adata_with_combined_recluster[adata_with_combined_recluster.obs['leiden_new'] == '0'])

# Count cluster 0 subclusters if reclustering was performed
cluster0_subclusters = []
if combined_recluster_results:
    cluster0_subclusters = [c for c in adata_with_combined_recluster.obs[cluster_column].cat.categories if c.startswith('0.')]

combined_summary = f"""Combined Analysis with Cluster 0 Reclustering - Comprehensive Summary
================================================================================

Dataset Overview:
- Total cells across all samples: {adata_with_combined_recluster.n_obs:,}
- Total genes: {adata_with_combined_recluster.n_vars:,}
- Number of samples: {len(samples)}
- Samples analyzed: {', '.join(sorted(samples))}

Original Clustering:
- Original cluster 0 cells: {total_cluster0_cells:,}

Combined Cluster 0 Reclustering Results:
{"- Reclustering performed: Yes" if combined_recluster_results else "- Reclustering performed: No"}
{f"- Number of subclusters created: {len(cluster0_subclusters)}" if combined_recluster_results else ""}
{f"- Subcluster labels: {sorted(cluster0_subclusters)}" if combined_recluster_results else ""}

Final Combined Analysis:
- Total final clusters: {len(cluster_counts)}
- Clusters analyzed for differential expression: {len(valid_clusters_for_de)}
- Analysis type: {"Combined reclustered" if combined_recluster_results else "Original clustering"}

Files Generated:
1. combined_umap_overview.png - Combined UMAP plots side by side
2. combined_spatial_overview.png - Combined spatial plots side by side
{f"2. combined_cluster0_subset_analysis.png - Cluster 0 subset analysis plots" if combined_recluster_results else ""}
3. {marker_genes_filename} - Complete combined marker gene results
4. {top_markers_filename} - Top marker genes per cluster (combined analysis)
5. {heatmap_filename} - Combined marker genes heatmap
6. {dotplot_filename} - Combined marker genes dotplot
7. {dendrogram_filename} - Combined hierarchical clustering dendrogram
8. {correlation_filename} - Combined cluster correlation heatmap
9. {correlation_csv_filename} - Combined correlation matrix data
10. cluster_composition_across_samples.csv - Cluster composition analysis
11. cluster_cell_count_heatmap.png - Cell count distribution heatmap
12. cluster_percentage_heatmap.png - Percentage distribution heatmap

Analysis Parameters:
- Combined differential expression method: Wilcoxon rank-sum test
- Number of marker genes per cluster: 100 (combined analysis)
- Distance metric for clustering: Correlation distance
- Linkage method: Ward
- Minimum cells per cluster threshold for DE: >1 cell
- Combined cluster 0 reclustering resolution: 0.8
- Neighborhood size for reclustering: 30

Key Insights:
- This analysis performs cluster 0 reclustering on ALL cluster 0 cells combined across samples
- Enables identification of shared and sample-specific subclusters within cluster 0
- Provides comprehensive cross-sample comparisons and marker gene analysis
- Sample composition analysis reveals cluster distribution patterns across samples

Notes:
- Individual sample analyses are available in the per-sample output directories
- This combined analysis complements the individual sample analyses
- Cluster 0 subclusters identified here represent shared patterns across all samples
"""

with open(os.path.join(output_dir, 'analysis_summary.txt'), 'w') as f:
    f.write(combined_summary)

print(combined_summary)
print(f"\nCombined analysis completed successfully!")
print(f"Results saved to: {output_dir}")

# Save the updated AnnData object with combined reclustering
if combined_recluster_results:
    print("Saving updated AnnData object with combined reclustering...")
    adata_with_combined_recluster.write_h5ad(os.path.join(output_dir, 'merged_data_with_combined_cluster0_reclustering.h5ad'))
    print("Updated AnnData object saved")

print("\n" + "="*80)
print("COMBINED ANALYSIS WITH CLUSTER 0 RECLUSTERING COMPLETED SUCCESSFULLY!")
print("="*80)
print(f"\nReclustering Summary:")
print(f"- Original Cluster 0 cells: {total_cluster0_cells}")
print(f"- New subclusters created: {len(cluster0_subclusters)}")
print(f"- Cluster labels updated: {len(combined_recluster_results['cluster0_subset'])} cells")
print(f"\nAll outputs saved to: {output_dir}")
print("\nKey files generated:")
print(f"- Marker genes: {marker_genes_filename}")
print(f"- Top markers: {top_markers_filename}")
print(f"- UMAP plots: UMAP_full.png")
print(f"- Heatmap: heatmap_marker_genes.png")
print(f"- Dotplot: combined_dotplot_marker_genes_with_cluster0_reclustered.png")
print(f"- Dendrogram: combined_cluster_dendrogram_with_cluster0_reclustered.png")
print(f"- Correlation heatmap: combined_cluster_correlation_with_cluster0_reclustered.png")
print(f"- Correlation matrix: combined_correlation_matrix_with_cluster0_reclustered.csv")

# %%