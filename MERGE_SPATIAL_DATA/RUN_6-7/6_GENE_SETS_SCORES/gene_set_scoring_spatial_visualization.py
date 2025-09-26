#!/usr/bin/env python3
"""
Gene Set Scoring and Spatial Visualization

This script performs gene set scoring at single cell level and creates spatial 
visualizations showing gene set scores as red color intensity across spatial 
transcriptomics samples.
"""

import scanpy as sc
import squidpy as sq
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Define paths
merged_h5ad_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
gene_sets_path = 'GENE_LISTS/output_combined/import_ready.csv'
output_dir = 'MERGE_SPATIAL_DATA/6_GENE_SETS_SCORES/plots'

# Create output directory
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

def load_gene_sets(gene_sets_path):
    """
    Load gene sets from CSV file.
    
    Parameters:
    -----------
    gene_sets_path : str
        Path to the gene sets CSV file
    
    Returns:
    --------
    dict
        Dictionary with gene set names as keys and lists of genes as values
    """
    print(f"Loading gene sets from: {gene_sets_path}")
    
    # Read the CSV file
    gene_sets_df = pd.read_csv(gene_sets_path, header=None, 
                              names=['gene_set', 'gene', 'col3', 'col4'])
    
    # Group genes by gene set
    gene_sets_dict = {}
    for gene_set, group in gene_sets_df.groupby('gene_set'):
        gene_sets_dict[gene_set] = group['gene'].tolist()
    
    print(f"Loaded {len(gene_sets_dict)} gene sets:")
    for gene_set, genes in gene_sets_dict.items():
        print(f"  - {gene_set}: {len(genes)} genes")
    
    return gene_sets_dict

def calculate_gene_set_scores(adata, gene_sets_dict):
    """
    Calculate gene set scores using scanpy's score_genes function.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object containing the gene expression data
    gene_sets_dict : dict
        Dictionary with gene set names as keys and lists of genes as values
    
    Returns:
    --------
    AnnData
        Updated AnnData object with gene set scores in .obs
    """
    print("Calculating gene set scores...")
    
    # Make a copy to avoid modifying the original
    adata_copy = adata.copy()
    
    # Calculate scores for each gene set
    for gene_set_name, genes in gene_sets_dict.items():
        print(f"  Scoring gene set: {gene_set_name}")
        
        # Filter genes to only those present in the dataset
        genes_in_data = [gene for gene in genes if gene in adata_copy.var.index]
        genes_missing = [gene for gene in genes if gene not in adata_copy.var.index]
        
        if len(genes_missing) > 0:
            print(f"    Warning: {len(genes_missing)} genes not found in dataset: {genes_missing[:5]}{'...' if len(genes_missing) > 5 else ''}")
        
        if len(genes_in_data) == 0:
            print(f"    Error: No genes from {gene_set_name} found in dataset. Skipping.")
            continue
        
        print(f"    Using {len(genes_in_data)}/{len(genes)} genes from {gene_set_name}")
        
        # Calculate gene set score using scanpy
        # The score will be stored in adata.obs with the gene set name
        sc.tl.score_genes(adata_copy, gene_list=genes_in_data, 
                         score_name=f"{gene_set_name}_score", use_raw=False)
    
    return adata_copy

def create_spatial_plots(adata, gene_sets_dict, output_dir):
    """
    Create spatial plots showing gene set scores as red color intensity.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names as keys and lists of genes as values
    output_dir : str
        Directory to save the plots
    """
    print("Creating spatial plots...")
    
    # Create spatial coordinates if not present
    if 'spatial' not in adata.obsm:
        adata.obsm['spatial'] = np.c_[adata.obs['center_x'], adata.obs['center_y']]
    
    # Create spatial metadata in .uns if not present
    if 'spatial' not in adata.uns:
        adata.uns['spatial'] = {}
        for sample_id in adata.obs['sample_id'].unique():
            adata.uns['spatial'][sample_id] = {
                'images': {}, 
                'scalefactors': {'spot_diameter_fullres': 1.0, 'hires': 1.0}
            }
    
    # Get unique samples
    samples = adata.obs['sample_id'].unique()
    
    # Create plots for each gene set and sample
    for gene_set_name in gene_sets_dict.keys():
        score_column = f"{gene_set_name}_score"
        
        # Skip if score column doesn't exist (gene set was skipped during scoring)
        if score_column not in adata.obs.columns:
            print(f"  Skipping {gene_set_name} - no scores calculated")
            continue
        
        print(f"  Plotting gene set: {gene_set_name}")
        
        # Create a subdirectory for this gene set
        gene_set_dir = os.path.join(output_dir, gene_set_name.replace(' ', '_').replace('/', '_'))
        os.makedirs(gene_set_dir, exist_ok=True)
        
        # Plot for all samples combined
        try:
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Use squidpy for spatial plotting with red colormap
            sq.pl.spatial_scatter(adata, color=score_column, 
                                shape=None, 
                                library_key='sample_id',
                                cmap='Reds',
                                ax=ax,
                                title=f"{gene_set_name} - All Samples")
            
            plt.savefig(os.path.join(gene_set_dir, f'{gene_set_name.replace(" ", "_")}_all_samples.png'), 
                       bbox_inches='tight', dpi=300)
            plt.close()
            
        except Exception as e:
            print(f"    Error plotting all samples for {gene_set_name}: {e}")
        
        # Plot for each sample individually
        for sample_id in samples:
            try:
                sample_adata = adata[adata.obs['sample_id'] == sample_id].copy()
                
                if len(sample_adata) == 0:
                    continue
                
                fig, ax = plt.subplots(figsize=(10, 8))
                
                # Get spatial coordinates
                spatial_coords = sample_adata.obsm['spatial']
                scores = sample_adata.obs[score_column]
                
                # Create scatter plot with red colormap
                scatter = ax.scatter(spatial_coords[:, 0], spatial_coords[:, 1], 
                                   c=scores, cmap='Reds', s=3, alpha=0.8)
                
                # Add colorbar
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label(f'{gene_set_name} Score', rotation=270, labelpad=15)
                
                # Set labels and title
                ax.set_xlabel('X coordinate')
                ax.set_ylabel('Y coordinate')
                ax.set_title(f'{gene_set_name} - Sample {sample_id}\n'
                           f'Score range: {scores.min():.3f} to {scores.max():.3f}')
                ax.set_aspect('equal')
                
                # Save the plot
                plt.savefig(os.path.join(gene_set_dir, f'{gene_set_name.replace(" ", "_")}_sample_{sample_id}.png'), 
                           bbox_inches='tight', dpi=300)
                plt.close()
                
            except Exception as e:
                print(f"    Error plotting sample {sample_id} for {gene_set_name}: {e}")

def create_summary_report(adata, gene_sets_dict, output_dir):
    """
    Create a summary report of gene set scores.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object with gene set scores
    gene_sets_dict : dict
        Dictionary with gene set names as keys and lists of genes as values
    output_dir : str
        Directory to save the report
    """
    print("Creating summary report...")
    
    # Create summary statistics
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
    
    # Create a formatted text report
    with open(os.path.join(output_dir, 'gene_set_scores_report.txt'), 'w') as f:
        f.write("Gene Set Scoring Summary Report\n")
        f.write("=" * 50 + "\n\n")
        
        f.write(f"Total gene sets analyzed: {len(gene_sets_dict)}\n")
        f.write(f"Total cells: {len(adata)}\n")
        f.write(f"Total samples: {len(adata.obs['sample_id'].unique())}\n\n")
        
        # Gene set overview
        f.write("Gene Sets Overview:\n")
        f.write("-" * 20 + "\n")
        for gene_set_name, genes in gene_sets_dict.items():
            score_column = f"{gene_set_name}_score"
            if score_column in adata.obs.columns:
                genes_found = len([g for g in genes if g in adata.var.index])
                f.write(f"{gene_set_name}: {genes_found}/{len(genes)} genes found\n")
            else:
                f.write(f"{gene_set_name}: SKIPPED (no genes found)\n")
        
        f.write(f"\nDetailed statistics saved to: gene_set_scores_summary.csv\n")
        f.write(f"Spatial plots saved to individual gene set directories\n")

def main():
    """Main function to run the gene set scoring and visualization pipeline."""
    
    print("Gene Set Scoring and Spatial Visualization Pipeline")
    print("=" * 60)
    
    # Load the merged AnnData object
    if os.path.exists(merged_h5ad_path):
        print(f"Loading merged data from: {merged_h5ad_path}")
        adata = sc.read_h5ad(merged_h5ad_path)
        print(f"Loaded data: {adata.n_obs} cells x {adata.n_vars} genes")
    else:
        print(f"Error: Merged data file not found at {merged_h5ad_path}")
        sys.exit(1)
    
    # Load gene sets
    if os.path.exists(gene_sets_path):
        gene_sets_dict = load_gene_sets(gene_sets_path)
    else:
        print(f"Error: Gene sets file not found at {gene_sets_path}")
        sys.exit(1)
    
    # Calculate gene set scores
    adata_scored = calculate_gene_set_scores(adata, gene_sets_dict)
    
    # Create spatial plots
    create_spatial_plots(adata_scored, gene_sets_dict, output_dir)
    
    # Create summary report
    create_summary_report(adata_scored, gene_sets_dict, output_dir)
    
    print("\n" + "=" * 60)
    print("Pipeline completed successfully!")
    print(f"Results saved to: {output_dir}")
    print("=" * 60)

if __name__ == "__main__":
    main()