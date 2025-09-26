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
output_path = os.path.join(output_dir, "p30_R4_tangram_integrated_v4.h5ad")

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

# --- Step 2a: Ensure Consistent Cell Type Annotations (Crucial for Matched Samples) ---
# For matched samples, ensuring annotations are consistent is key.
# This assumes 'cell_type' is in `adata_sc.obs`. We can project these labels to the spatial data.
# If adata_sp lacks annotations, uncomment the following lines:
# print("\nAnnotating spatial data with cell types from single-cell reference...")
# tg.annotate_cell_types(
#     adata_sp=adata_sp,
#     adata_sc=adata_sc,
#     annotation='cell_type' # The column in adata_sc.obs with cell type labels
# )
# assert 'cell_type' in adata_sp.obs.columns

# --- Step 2b: Quality Control using Matched Data ---
# With matched data, we can identify problematic genes by comparing their expression
# profiles across platforms. Genes with low correlation may indicate technical issues.
print("\nIdentifying problematic genes by correlating expression across platforms...")
# This step requires consistent cell type annotations (see Step 2a).
if 'cell_type' not in adata_sc.obs.columns or 'cell_type' not in adata_sp.obs.columns:
    raise ValueError("Both AnnData objects must have a 'cell_type' column in .obs for correlation analysis.")

# Calculate average expression per cell type for each platform
sc_avg_by_type = adata_sc[:, common_genes].to_df().groupby(adata_sc.obs['cell_type']).mean()
sp_avg_by_type = adata_sp[:, common_genes].to_df().groupby(adata_sp.obs['cell_type']).mean()

# Align cell types and genes
common_cell_types = list(set(sc_avg_by_type.index) & set(sp_avg_by_type.index))
sc_avg_aligned = sc_avg_by_type.loc[common_cell_types]
sp_avg_aligned = sp_avg_by_type.loc[common_cell_types]

# Calculate correlation for each gene
gene_correlations = sc_avg_aligned.corrwith(sp_avg_aligned, axis=0)
correlation_threshold = 0.5
problem_genes = set(gene_correlations[gene_correlations < correlation_threshold].index)
print(f"Identified {len(problem_genes)} problematic genes with correlation < {correlation_threshold}.")

# --- Step 2c: Strategic Gene Selection for Training ---
# For robust mapping, we train on highly variable genes among the high-confidence candidates.
print("\nSelecting highly variable genes from high-confidence set for training...")
training_candidates = [g for g in common_genes if g not in problem_genes]

# Find highly variable genes among the candidates.
adata_sc_training_candidates = adata_sc[:, training_candidates].copy()
sc.pp.highly_variable_genes(adata_sc_training_candidates, n_top_genes=1000, flavor='seurat')
training_genes = adata_sc_training_candidates.var_names[adata_sc_training_candidates.var.highly_variable].tolist()
print(f"Using {len(training_genes)} high-confidence variable genes for training.")

# %%
# Pre-process AnnDatas for Tangram.
# This modifies the AnnDatas in-place, so we use copies for mapping.
print("\nPreparing data for Tangram mapping...")
adata_sc_map = adata_sc.copy()
adata_sp_map = adata_sp.copy()
tg.pp_adatas(adata_sc_map, adata_sp_map, genes=common_genes)

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
    train_genes=training_genes,
    filter_gene_counts=True
)

# --- Advanced Option: Cell Type-Aware Mapping ---
# For higher accuracy with matched samples, you can map each cell type individually.
# This requires uncommenting and running the following loop instead of the block above.
#
# print("\nPerforming cell type-aware mapping...")
# ad_map = None # Initialize the map object
# for cell_type in adata_sc.obs['cell_type'].unique():
#     print(f"Mapping cell type: {cell_type}")
#     # Create subsets for the current cell type
#     sc_subset = adata_sc_map[adata_sc_map.obs['cell_type'] == cell_type].copy()
#     sp_subset = adata_sp_map[adata_sp_map.obs['cell_type'] == cell_type].copy()
#
#     # Map the subset
#     ad_map_subset = tg.map_cells_to_space(
#         adata_sc=sc_subset,
#         adata_sp=sp_subset,
#         device=device,
#         mode='cells',
#         density_prior='uniform',
#         num_epochs=500,
#         train_genes=training_genes
#     )
#     # To-do: Combine the individual mapping results into a single `ad_map` object.
#     # This typically involves careful concatenation of the mapping matrices.
#     # For now, we proceed with the combined mapping.

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
# 5. Validation using Matched-Sample Advantage (Optional)
# With matched samples, we can perform more powerful validation.
print("\nPerforming validation using matched sample properties...")

# --- 5a. Cell Type Composition Validation ---
# Compare the frequency of cell types between the two datasets. They should be similar.
# print("Validating cell type composition...")
# sp_freq = ad_map.obs['cell_type'].value_counts(normalize=True)
# sc_freq = adata_sc.obs['cell_type'].value_counts(normalize=True)
# common_types = list(set(sp_freq.index) & set(sc_freq.index))
# composition_corr = np.corrcoef(sp_freq[common_types], sc_freq[common_types])[0, 1]
# print(f"Correlation of cell type frequencies: {composition_corr:.4f}")
# if composition_corr < 0.8:
#     print("Warning: Cell type composition correlation is low. Review annotations or mapping.")

# --- 5b. Gene Expression Pattern Validation ---
# Compare expression profiles for high-quality genes.
# print("\nValidating gene expression patterns...")
# validation_genes = [g for g in common_genes if g not in problem_genes][:20] # Check top 20 genes
# sc_means = adata_sc[:, validation_genes].to_df().groupby(adata_sc.obs['cell_type']).mean()
# projected_means = adata_ge[:, validation_genes].to_df().groupby(ad_map.obs['cell_type']).mean()
#
# # Align and calculate correlation
# common_types_val = list(set(sc_means.index) & set(projected_means.index))
# gene_pattern_corr = sc_means.loc[common_types_val].corrwith(projected_means.loc[common_types_val], axis=0)
# print(f"Average gene expression pattern correlation: {gene_pattern_corr.mean():.4f}")
# print("Correlations for individual genes:")
# print(gene_pattern_corr)

# %%
# 6. Save the integrated data
print(f"\nSaving integrated AnnData to {output_path}")
# The `adata_ge` object contains the spatial locations from MERSCOPE with the
# corrected and expanded gene expression matrix from the scRNA-seq reference.
adata_ge.write(output_path)

print("\nTangram integration script finished successfully!")
print(f"The final integrated data is saved at: {output_path}")