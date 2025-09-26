# %%
import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import warnings
from celltypist import models, annotate
import seaborn as sns
import random

# %%
# --- Configuration ---

# Set seeds for reproducibility
random_seed = 0
np.random.seed(random_seed)
random.seed(random_seed)
sc.settings.seed = random_seed

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Define input and output paths
INPUT_H5AD_FILE = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
OUTPUT_DIR = 'MERGE_SPATIAL_DATA/celltypist_annotation'
OUTPUT_H5AD_FILE = os.path.join(OUTPUT_DIR, 'merged_spatial_hippo_annotated.h5ad')
MODEL_DIR = os.path.join(project_root, "MERGE_SPATIAL_DATA", "models")

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(MODEL_DIR, exist_ok=True)

# CellTypist parameters
MAJORITY_VOTING = True
PROB_THRESHOLD = 0.5
LEIDEN_KEY = 'leiden_new'  # Based on your merge script

# Model configuration - focusing on Dentate Gyrus
MODEL_NAME = "Mouse_Dentate_Gyrus"
MODEL_INFO = {
    "description": "Cell types from the dentate gyrus in perinatal, juvenile, and adult mice",
    "cell_types": 24,
    "version": "v1",
    "reference": "https://doi.org/10.1038/s41593-017-0056-2"
}

print("="*60)
print("CELLTYPIST ANNOTATION FOR SPATIAL HIPPOCAMPUS DATA")
print("="*60)
print(f"Model: {MODEL_NAME}")
print(f"Description: {MODEL_INFO['description']}")
print(f"Input file: {INPUT_H5AD_FILE}")
print(f"Output directory: {OUTPUT_DIR}")
print("="*60)

# %% [markdown]
# ## Load Merged Spatial Data

# %%
print("\nLoading merged spatial dataset...")
if not os.path.exists(INPUT_H5AD_FILE):
    raise FileNotFoundError(f"Input file not found: {INPUT_H5AD_FILE}. Please run merge_hippo_inclusive.py first.")

adata = sc.read_h5ad(INPUT_H5AD_FILE)
print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes from {INPUT_H5AD_FILE}")

# Display basic information about the dataset
print(f"\nDataset overview:")
print(f"  - Number of cells: {adata.n_obs}")
print(f"  - Number of genes: {adata.n_vars}")
print(f"  - Samples: {adata.obs['sample_id'].unique()}")
print(f"  - Time points: {adata.obs['time_point'].unique()}")
print(f"  - Types: {adata.obs['type'].unique()}")

# %%
# Basic check for expected structure
if 'X_umap' not in adata.obsm:
    warnings.warn("UMAP coordinates ('X_umap') not found. Plots might fail.")
if LEIDEN_KEY not in adata.obs:
    warnings.warn(f"Leiden clustering key '{LEIDEN_KEY}' not found. Using available leiden key.")
    # Find available leiden keys
    leiden_keys = [col for col in adata.obs.columns if col.startswith('leiden')]
    if leiden_keys:
        LEIDEN_KEY = leiden_keys[0]
        print(f"Using leiden key: {LEIDEN_KEY}")
    else:
        warnings.warn("No leiden clustering found in data.")

# Check for spatial coordinates
if 'spatial' not in adata.obsm:
    print("Creating spatial coordinates from center_x and center_y...")
    adata.obsm['spatial'] = np.c_[adata.obs['center_x'], adata.obs['center_y']]

# %% [markdown]
# ## Prepare Data for CellTypist

# %%
# Create a copy for CellTypist analysis
adata_ct = adata.copy()

print("\nPreparing data for CellTypist analysis...")

# CellTypist expects log1p normalized data with a target sum (e.g., 10,000).
# Check if we have raw counts available
if adata_ct.raw is not None:
    print("Using raw counts from .raw slot for currently filtered genes.")
    # Use raw counts for the current genes in adata_ct.var_names
    adata_ct.X = adata_ct.raw[:, adata_ct.var_names].X.copy()
else:
    print("No raw counts found in .raw slot. Using current .X matrix.")
    # Assume current X contains counts or suitable data for normalization
    pass

# Check if data appears to be already normalized/processed
x_max = adata_ct.X.max()
x_mean_sum = np.mean(adata_ct.X.sum(axis=1))

print(f"Data statistics before preprocessing:")
print(f"  - Maximum value: {x_max:.2f}")
print(f"  - Mean sum per cell: {x_mean_sum:.2f}")

# If data appears to be already log-transformed (max value < 15), we might need to exponentiate
if x_max < 15 and x_mean_sum < 50:
    print("Data appears to be log-transformed. Converting back to counts...")
    adata_ct.X = np.expm1(adata_ct.X)

# 1. Normalize total counts per cell to 10,000
print("Normalizing total counts to 10,000 per cell...")
sc.pp.normalize_total(adata_ct, target_sum=10000)

# 2. Log1p transform the data
print("Applying log1p transformation...")
sc.pp.log1p(adata_ct)

# Verification step
print("Verifying data format for CellTypist...")
x_max_final = adata_ct.X.max()
if x_max_final > 20:
    warnings.warn(f"Maximum value in .X ({x_max_final:.2f}) seems high for log1p data.")

print(f"Final data statistics:")
print(f"  - Maximum value after preprocessing: {x_max_final:.2f}")
print(f"  - Mean sum per cell after preprocessing: {np.mean(adata_ct.X.sum(axis=1)):.2f}")
print("Data preparation for CellTypist complete.")

# %% [markdown]
# ## Run CellTypist Annotation

# %%
print(f"\nRunning CellTypist annotation with {MODEL_NAME}")
print(f"Parameters: majority_voting={MAJORITY_VOTING}, prob_threshold={PROB_THRESHOLD}")

# Set up model path
model_path = os.path.join(MODEL_DIR, f"{MODEL_NAME}.pkl")

# Download or load model
if not os.path.exists(model_path):
    print(f"Model file not found at {model_path}. Attempting to download...")
    try:
        models.celltypist_models(force_update=False)
        model = models.Model.load(MODEL_NAME)
        model.write(model_path)
        print(f"Model downloaded and saved to {model_path}")
    except Exception as e:
        raise FileNotFoundError(f"Model file not found at {model_path} and download failed: {e}")
else:
    print(f"Loading model from {model_path}")

# Run annotation
print("Running CellTypist annotation...")
try:
    predictions = annotate(
        adata_ct,
        model=model_path,
        majority_voting=MAJORITY_VOTING,
        mode='prob match',
        p_thres=PROB_THRESHOLD
    )

    # Add annotations with DG prefix (Dentate Gyrus)
    prefix = 'DG_'
    predictions.to_adata(adata_ct, prefix=prefix)
    print(f"Added CellTypist predictions with prefix '{prefix}' to adata.obs")

    # Display available annotation columns
    annotation_cols = [col for col in adata_ct.obs.columns if col.startswith(prefix)]
    print(f"Added annotation columns: {annotation_cols}")

except Exception as e:
    print(f"ERROR running CellTypist for model {MODEL_NAME}: {e}")
    sys.exit(1)

# %% [markdown]
# ## Generate Visualizations

# %%
print("\nGenerating visualizations...")
plot_subdir = os.path.join(OUTPUT_DIR, "plots")
os.makedirs(plot_subdir, exist_ok=True)
sc.settings.figdir = plot_subdir

# Set up annotation column names
cell_type_col = f'{prefix}majority_voting'
conf_score_col = f'{prefix}conf_score'

# 1. UMAP plots
if 'X_umap' in adata_ct.obsm:
    print("Generating UMAP plots...")
    
    # Cell type annotation UMAP
    if cell_type_col in adata_ct.obs.columns:
        sc.pl.umap(adata_ct, color=cell_type_col, 
                  title=f"Cell Type Annotation ({MODEL_NAME})",
                  legend_loc='on data', 
                  legend_fontsize='xx-small',
                  save='_cell_types_umap.png', show=False)
        print("Generated UMAP plot for cell types.")
    
    # Confidence score UMAP
    if conf_score_col in adata_ct.obs.columns:
        sc.pl.umap(adata_ct, color=conf_score_col, 
                  title=f"Annotation Confidence Score ({MODEL_NAME})",
                  cmap='viridis',
                  save='_confidence_score_umap.png', show=False)
        print("Generated UMAP plot for confidence scores.")
    
    # Original leiden clusters for comparison
    if LEIDEN_KEY in adata_ct.obs.columns:
        sc.pl.umap(adata_ct, color=LEIDEN_KEY,
                  title="Original Leiden Clusters",
                  save='_leiden_clusters_umap.png', show=False)
        print("Generated UMAP plot for original leiden clusters.")

# 2. Spatial plots using squidpy
if 'spatial' in adata_ct.obsm:
    print("Generating spatial plots...")
    
    # Cell type annotation spatial plot
    if cell_type_col in adata_ct.obs.columns:
        sq.pl.spatial_scatter(adata_ct, color=cell_type_col, shape=None,
                             library_key='sample_id',
                             save='cell_types_spatial.png')
        print("Generated spatial plot for cell types.")
    
    # Confidence score spatial plot
    if conf_score_col in adata_ct.obs.columns:
        sq.pl.spatial_scatter(adata_ct, color=conf_score_col, shape=None,
                             library_key='sample_id', 
                             save='confidence_score_spatial.png')
        print("Generated spatial plot for confidence scores.")
    
    # Original leiden clusters spatial plot
    if LEIDEN_KEY in adata_ct.obs.columns:
        sq.pl.spatial_scatter(adata_ct, color=LEIDEN_KEY, shape=None,
                             library_key='sample_id',
                             save='leiden_clusters_spatial.png')
        print("Generated spatial plot for original leiden clusters.")

# %% [markdown]
# ## Generate Summary Statistics

# %%
print("\nGenerating summary statistics...")

# Cell type summary
if cell_type_col in adata_ct.obs.columns:
    cell_type_counts = adata_ct.obs[cell_type_col].value_counts()
    cell_type_df = pd.DataFrame({
        'cell_type': cell_type_counts.index,
        'cell_count': cell_type_counts.values,
        'percentage': (cell_type_counts.values / cell_type_counts.sum() * 100).round(2)
    })
    cell_type_df = cell_type_df.sort_values('cell_count', ascending=False).reset_index(drop=True)
    
    # Save overall summary
    summary_file = os.path.join(OUTPUT_DIR, f"{MODEL_NAME}_cell_type_summary.csv")
    cell_type_df.to_csv(summary_file, index=False)
    print(f"Saved cell type summary to {summary_file}")
    
    # Print summary
    print(f"\nCell type annotation summary:")
    print(cell_type_df.head(10))

# Summary by sample
if cell_type_col in adata_ct.obs.columns and 'sample_id' in adata_ct.obs.columns:
    sample_summary = adata_ct.obs.groupby(['sample_id', cell_type_col]).size().unstack(fill_value=0)
    sample_summary_pct = sample_summary.div(sample_summary.sum(axis=1), axis=0) * 100
    
    # Save sample summaries
    sample_counts_file = os.path.join(OUTPUT_DIR, f"{MODEL_NAME}_cell_type_by_sample_counts.csv")
    sample_pct_file = os.path.join(OUTPUT_DIR, f"{MODEL_NAME}_cell_type_by_sample_percentages.csv")
    
    sample_summary.to_csv(sample_counts_file)
    sample_summary_pct.to_csv(sample_pct_file)
    
    print(f"Saved sample-wise cell type counts to {sample_counts_file}")
    print(f"Saved sample-wise cell type percentages to {sample_pct_file}")

# Confidence score summary
if conf_score_col in adata_ct.obs.columns:
    conf_stats = adata_ct.obs[conf_score_col].describe()
    print(f"\nConfidence score statistics:")
    print(conf_stats)
    
    # Save confidence stats
    conf_stats_file = os.path.join(OUTPUT_DIR, f"{MODEL_NAME}_confidence_statistics.csv")
    conf_stats.to_csv(conf_stats_file, header=['confidence_score'])
    print(f"Saved confidence statistics to {conf_stats_file}")

# %% [markdown]
# ## Compare with Original Clustering

# %%
print("\nComparing CellTypist annotations with original Leiden clustering...")

if LEIDEN_KEY in adata_ct.obs.columns and cell_type_col in adata_ct.obs.columns:
    # Create contingency table
    cluster_comparison = pd.crosstab(
        adata_ct.obs[LEIDEN_KEY],
        adata_ct.obs[cell_type_col],
        margins=True
    )
    
    # Save comparison table
    comparison_file = os.path.join(OUTPUT_DIR, f"{MODEL_NAME}_vs_leiden_comparison.csv")
    cluster_comparison.to_csv(comparison_file)
    print(f"Saved cluster comparison table to {comparison_file}")
    
    # Generate comparison heatmap
    try:
        plt.figure(figsize=(15, 10))
        # Remove margins for heatmap
        heatmap_data = cluster_comparison.iloc[:-1, :-1]  
        sns.heatmap(heatmap_data, cmap='YlOrRd', annot=True, fmt='g', linewidths=.5)
        plt.title(f'Cell Type Annotation vs Leiden Clustering\n({MODEL_NAME})')
        plt.xlabel('CellTypist Annotation')
        plt.ylabel('Leiden Cluster')
        plt.tight_layout()
        
        heatmap_file = os.path.join(plot_subdir, "cluster_comparison_heatmap.png")
        plt.savefig(heatmap_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Generated comparison heatmap: {heatmap_file}")
    except Exception as e:
        print(f"Could not generate comparison heatmap: {e}")

# %% [markdown]
# ## Save Annotated Data

# %%
print(f"\nSaving annotated data to {OUTPUT_H5AD_FILE}")
try:
    os.makedirs(os.path.dirname(OUTPUT_H5AD_FILE), exist_ok=True)
    adata_ct.write(OUTPUT_H5AD_FILE)
    print("Successfully saved annotated AnnData object.")
except Exception as e:
    print(f"Error saving annotated AnnData object: {e}")

# %% [markdown]
# ## Final Summary

# %%
print(f"\n{'='*60}")
print("CELLTYPIST ANNOTATION COMPLETED")
print(f"{'='*60}")
print(f"Input file: {INPUT_H5AD_FILE}")
print(f"Output file: {OUTPUT_H5AD_FILE}")
print(f"Number of cells: {adata_ct.n_obs}")
print(f"Number of genes: {adata_ct.n_vars}")
print(f"Model used: {MODEL_NAME}")
print(f"Leiden clustering key: {LEIDEN_KEY}")

if cell_type_col in adata_ct.obs.columns:
    n_cell_types = adata_ct.obs[cell_type_col].nunique()
    print(f"Number of cell types identified: {n_cell_types}")
    
    if conf_score_col in adata_ct.obs.columns:
        mean_confidence = adata_ct.obs[conf_score_col].mean()
        print(f"Mean confidence score: {mean_confidence:.3f}")

print(f"\nOutput files:")
print(f"  - Annotated data: {OUTPUT_H5AD_FILE}")
print(f"  - Plots directory: {plot_subdir}")
print(f"  - Summary files: {OUTPUT_DIR}")
print(f"{'='*60}")
print("Analysis complete!")

# %%