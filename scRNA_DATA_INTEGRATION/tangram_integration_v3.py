# %%
import os
import scanpy as sc
import tangram as tg
import pandas as pd
import numpy as np

# Set working directory
os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation")

# %%
# Define file paths
spatial_data_path = "DATA/p30-E165/R4/data.h5ad"
sc_data_path = "DATA/RNA_counts/Emx1_Ctrl/Emx1_Ctrl_processed.h5ad"
output_dir = "DATA_INTEGRATION"
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, "p30_R4_tangram_integrated_v3.h5ad")

# %%
# 1. Load Data
print("Loading spatial data...")
adata_sp = sc.read_h5ad(spatial_data_path)
print("Spatial data loaded:")
print(adata_sp)

print("\nLoading single-cell data...")
adata_sc = sc.read_h5ad(sc_data_path)
print("Single-cell data loaded:")
print(adata_sc)

# %%
# Ensure spatial coordinates are in .obsm
if 'spatial' not in adata_sp.obsm_keys():
    if 'center_x' in adata_sp.obs.columns and 'center_y' in adata_sp.obs.columns:
        adata_sp.obsm['spatial'] = adata_sp.obs[['center_x', 'center_y']].values
        print("Added spatial coordinates to adata_sp.obsm['spatial']")
    else:
        raise KeyError("Spatial coordinates ('center_x', 'center_y') not found in adata_sp.obs")
# %%
# 2. Pre-process Data for Tangram
print("\nPreprocessing data for Tangram...")

# Find common genes between the two datasets.
common_genes = list(set(adata_sc.var_names) & set(adata_sp.var_names))
print(f"Found {len(common_genes)} common genes.")

# --- MERSCOPE-specific gene quality assessment ---
# As per Tangram recommendations for image-based data, we assess gene quality.
# Problematic genes (e.g., with low detection rate) are excluded from training.
print("\nAssessing MERSCOPE gene quality to select high-confidence genes...")

# Calculate detection rate for each gene in the spatial data
detection_rate = np.array((adata_sp[:, common_genes].X > 0).mean(axis=0)).flatten()
detection_df = pd.DataFrame({'gene': common_genes, 'detection_rate': detection_rate})

# Flag genes with very low detection rate (e.g., < 0.1%) as problematic
low_detection_threshold = 0.001
problem_genes = set(detection_df[detection_df['detection_rate'] < low_detection_threshold]['gene'])
print(f"Identified {len(problem_genes)} problematic genes with detection rate < {low_detection_threshold}.")

# --- Strategic Gene Selection for Training ---
# For robust mapping, we train on highly variable genes among the high-confidence candidates.
print("\nSelecting highly variable genes from high-confidence set for training...")
training_candidates = [g for g in common_genes if g not in problem_genes]

# Find highly variable genes among the candidates.
adata_sc_training_candidates = adata_sc[:, training_candidates].copy()
sc.pp.highly_variable_genes(adata_sc_training_candidates, n_top_genes=1000, flavor='seurat')
training_genes = adata_sc_training_candidates.var_names[adata_sc_training_candidates.var.highly_variable].tolist()

# For best results, ensure key cell-type markers are included in the training set.
# Example:
# key_markers = {'Slc17a7', 'Gad1', 'Mog', 'Pecam1'}
# training_genes = list(set(training_genes) | (key_markers & set(common_genes)))

print(f"Using {len(training_genes)} high-confidence variable genes for training.")

# %%
# Pre-process AnnDatas for Tangram.
# This modifies the AnnDatas in-place, so we use copies for mapping.
print("\nPreparing data for Tangram mapping...")
# Subset the data to only include training genes before passing to Tangram
adata_sc_map = adata_sc[:, training_genes].copy()
adata_sp_map = adata_sp[:, training_genes].copy()
# Now that the data is subset, we can let pp_adatas find all genes for training
# and compute density priors.
tg.pp_adatas(adata_sc_map, adata_sp_map, genes=None)

# %%
# 3. Map single cells to space
print("\nMapping single cells to space...")
# Check for GPU availability, otherwise use CPU
try:
    import torch
    if torch.cuda.is_available():
        device = "cuda:0"
        print("Using GPU for mapping.")
    else:
        device = "cpu"
        print("Using CPU for mapping.")
except ImportError:
    device = "cpu"
    print("PyTorch not found. Using CPU for mapping.")

ad_map = tg.map_cells_to_space(
    adata_sc=adata_sc_map,
    adata_sp=adata_sp_map,
    device=device,
    mode='cells',
    density_prior='uniform',  # MERSCOPE has single cells, so uniform prior is appropriate
    num_epochs=1000,
)

print("Mapping complete.")
print("Resulting mapping AnnData:")
print(ad_map)

# %%
# 4. Project gene expression for correction and expansion
print("\nProjecting gene expression onto spatial data...")
# Project ALL genes from the original scRNA-seq data.
# This achieves two goals:
# 1. Correction: The expression of common genes is corrected based on the mapping.
# 2. Expansion: The expression of genes in the scRNA-seq reference but not in the
#    MERSCOPE panel are imputed onto the spatial map.
adata_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_sc)
print("Gene projection complete.")
print("Projected gene expression AnnData (with corrected/expanded genes):")
print(adata_ge)

# %%
# 5. Validation (Optional)
# The Tangram paper recommends cross-validation to ensure the mapping is robust.
# This involves holding out a gene, re-mapping, predicting, and checking the correlation.
# This can be computationally intensive.

# print("\nRunning cross-validation...")
# correlations = []
# for gene in training_genes[:10]: # Example: validate first 10 genes
#     print(f"Validating on gene: {gene}")
#     # Leave out gene
#     train_without = [g for g in training_genes if g != gene]
    
#     # Re-map
#     ad_map_cv = tg.map_cells_to_space(
#         adata_sc=adata_sc_map,
#         adata_sp=adata_sp_map,
#         device=device,
#         mode='cells',
#         density_prior='uniform',
#         num_epochs=500, # Fewer epochs for speed
#         cv_train_genes=train_without,
#     )
    
#     # Predict and compare
#     pred_cv = tg.project_genes(ad_map_cv, adata_sc=adata_sc_map)
#     true_expr = adata_sp_map[:, gene].X.squeeze()
#     pred_expr = pred_cv[:, gene].X.squeeze()
    
#     # To calculate spatial correlation, you might need a library like squidpy
#     # import squidpy as sq
#     # corr, pvalue = sq.metrics.spatial_corr(adata_sp_map, gene_x=gene, gene_y=gene, layer_x=None, layer_y='tangram_pred')
#     # correlations.append(corr)

# print(f"Cross-validation correlations: {correlations}")

# Biological validation should also be performed, e.g., checking if known
# layer markers show expected spatial patterns in the projected data.

# %%
# 6. Save the integrated data
print(f"\nSaving integrated AnnData to {output_path}")
adata_ge.write(output_path)

print("\nTangram integration script finished successfully!")
print(f"The final integrated data is saved at: {output_path}")