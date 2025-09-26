# %%
import scanpy as sc
import squidpy as sq
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Define the path to the merged AnnData file and gene sets
merged_h5ad_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
gene_sets_path = 'GENE_LISTS/output_combined/import_ready.csv'

# %%
# Load the merged AnnData object
if os.path.exists(merged_h5ad_path):
    print(f"Loading merged data from: {merged_h5ad_path}")
    adata = sc.read_h5ad(merged_h5ad_path)
    print("Merged data loaded successfully.")
else:
    print(f"Error: Merged data file not found at {merged_h5ad_path}")
    sys.exit()

# %%
# Load gene sets
print(f"Loading gene sets from: {gene_sets_path}")
gene_sets_df = pd.read_csv(gene_sets_path, header=None, names=['gene_set', 'gene', 'col3', 'col4'])

# Group genes by gene set
gene_sets_dict = {}
for gene_set, group in gene_sets_df.groupby('gene_set'):
    gene_sets_dict[gene_set] = group['gene'].tolist()

print(f"Loaded {len(gene_sets_dict)} gene sets:")
for gene_set, genes in gene_sets_dict.items():
    print(f"  - {gene_set}: {len(genes)} genes")

# %%
# Create spatial coordinates
adata.obsm['spatial'] = np.c_[adata.obs['center_x'], adata.obs['center_y']]

# %%
# Create spatial metadata in .uns
adata.uns['spatial'] = {}
for sample_id in adata.obs['sample_id'].unique():
    adata.uns['spatial'][sample_id] = {'images': {}, 'scalefactors': {'spot_diameter_fullres': 1.0, 'hires': 1.0}}

# %%
# Calculate gene set scores
print("Calculating gene set scores...")

for gene_set_name, genes in gene_sets_dict.items():
    print(f"  Scoring gene set: {gene_set_name}")
    
    # Filter genes to only those present in the dataset
    genes_in_data = [gene for gene in genes if gene in adata.var.index]
    genes_missing = [gene for gene in genes if gene not in adata.var.index]
    
    if len(genes_missing) > 0:
        print(f"    Warning: {len(genes_missing)} genes not found in dataset")
    
    if len(genes_in_data) == 0:
        print(f"    Error: No genes from {gene_set_name} found in dataset. Skipping.")
        continue
    
    print(f"    Using {len(genes_in_data)}/{len(genes)} genes from {gene_set_name}")
    
    # Calculate gene set score using scanpy
    sc.tl.score_genes(adata, gene_list=genes_in_data, 
                     score_name=f"{gene_set_name}_score", use_raw=False)

# %%
# Define output directory for plots
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/spatial_plots'
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

# %%
# Create spatial plots for each gene set and sample
samples = adata.obs['sample_id'].unique()

for gene_set_name in gene_sets_dict.keys():
    score_column = f"{gene_set_name}_score"
    
    # Skip if score column doesn't exist (gene set was skipped during scoring)
    if score_column not in adata.obs.columns:
        print(f"Skipping {gene_set_name} - no scores calculated")
        continue
    
    print(f"Creating spatial plots for gene set: {gene_set_name}")
    
    # Clean gene set name for filename
    clean_name = gene_set_name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
    
    # Plot for each sample individually
    for sample_id in samples:
        sample_adata = adata[adata.obs['sample_id'] == sample_id].copy()
        
        if len(sample_adata) == 0:
            continue
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Get spatial coordinates and scores
        spatial_coords = sample_adata.obsm['spatial']
        scores = sample_adata.obs[score_column]
        
        # Create scatter plot with red colormap
        scatter = ax.scatter(spatial_coords[:, 0], spatial_coords[:, 1], 
                           c=scores, cmap='Reds', s=5, alpha=0.8, 
                           vmin=scores.quantile(0.01), vmax=scores.quantile(0.99))
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
        cbar.set_label(f'{gene_set_name} Score', rotation=270, labelpad=20)
        
        # Set labels and title
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_title(f'{gene_set_name} - Sample {sample_id}\n'
                   f'Score range: {scores.min():.3f} to {scores.max():.3f}')
        ax.set_aspect('equal')
        
        # Save the plot
        plot_filename = f"{clean_name}_sample_{sample_id}.png"
        fig.savefig(os.path.join(output_dir, plot_filename), bbox_inches='tight', dpi=300)
        plt.close(fig)
        
        print(f"  Saved: {plot_filename}")

# %%
# Create summary statistics
print("\nCreating summary report...")
summary_data = []

for gene_set_name in gene_sets_dict.keys():
    score_column = f"{gene_set_name}_score"
    
    if score_column not in adata.obs.columns:
        continue
    
    scores = adata.obs[score_column]
    
    # Overall statistics
    overall_stats = {
        'Gene_Set': gene_set_name,
        'Sample': 'All',
        'N_Cells': len(scores),
        'Mean_Score': scores.mean(),
        'Std_Score': scores.std(),
        'Min_Score': scores.min(),
        'Max_Score': scores.max(),
        'N_Genes_in_Set': len(gene_sets_dict[gene_set_name]),
        'N_Genes_Found': len([g for g in gene_sets_dict[gene_set_name] if g in adata.var.index])
    }
    summary_data.append(overall_stats)
    
    # Per-sample statistics
    for sample_id in adata.obs['sample_id'].unique():
        sample_scores = adata.obs[adata.obs['sample_id'] == sample_id][score_column]
        
        sample_stats = {
            'Gene_Set': gene_set_name,
            'Sample': sample_id,
            'N_Cells': len(sample_scores),
            'Mean_Score': sample_scores.mean(),
            'Std_Score': sample_scores.std(),
            'Min_Score': sample_scores.min(),
            'Max_Score': sample_scores.max(),
            'N_Genes_in_Set': len(gene_sets_dict[gene_set_name]),
            'N_Genes_Found': len([g for g in gene_sets_dict[gene_set_name] if g in adata.var.index])
        }
        summary_data.append(sample_stats)

# Convert to DataFrame and save
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(os.path.join(output_dir, 'gene_set_scores_summary.csv'), index=False)

print(f"Summary saved to: {os.path.join(output_dir, 'gene_set_scores_summary.csv')}")

# %%
print("\n" + "="*60)
print("GENE SET SCORING AND SPATIAL VISUALIZATION COMPLETED")
print("="*60)
print(f"Results saved to: {output_dir}")
print(f"- Spatial plots for each gene set and sample")
print(f"- Summary statistics: gene_set_scores_summary.csv")
print("="*60)

# %%