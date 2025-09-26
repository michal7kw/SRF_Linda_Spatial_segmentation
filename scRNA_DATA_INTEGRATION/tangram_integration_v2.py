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
output_path = os.path.join(output_dir, "p30_R4_tangram_integrated_v2.h5ad")

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

# --- Step 2a: Quality Control & Gene Selection ---
# The script will now automatically choose the QC method based on data availability.
# Method 1: If cell type annotations are present (supervised)
# Method 2: If annotations are absent (unsupervised)

problem_genes = set()
annotations_present = 'cell_type' in adata_sc.obs.columns and 'cell_type' in adata_sp.obs.columns

if annotations_present:
    print("\nAnnotations found. Using supervised QC for problematic genes...")
    # With matched data, we can identify problematic genes by comparing their expression
    # profiles across platforms. Genes with low correlation may indicate technical issues.
    sc_avg_by_type = adata_sc[:, common_genes].to_df().groupby(adata_sc.obs['cell_type']).mean()
    sp_avg_by_type = adata_sp[:, common_genes].to_df().groupby(adata_sp.obs['cell_type']).mean()

    common_cell_types = list(set(sc_avg_by_type.index) & set(sp_avg_by_type.index))
    sc_avg_aligned = sc_avg_by_type.loc[common_cell_types]
    sp_avg_aligned = sp_avg_by_type.loc[common_cell_types]

    gene_correlations = sc_avg_aligned.corrwith(sp_avg_aligned, axis=0)
    correlation_threshold = 0.5
    problem_genes = set(gene_correlations[gene_correlations < correlation_threshold].index)
    print(f"Identified {len(problem_genes)} problematic genes with correlation < {correlation_threshold}.")

else:
    print("\nNo cell type annotations found. Using unsupervised QC for problematic genes...")
    # Compare expression distributions between platforms to find genes with major discrepancies.
    from scipy.stats import ks_2samp
    for gene in common_genes:
        sc_expr = adata_sc[:, gene].X.flatten()
        sp_expr = adata_sp[:, gene].X.flatten()
        
        # Check for near-complete dropout in spatial data or loss of variation
        if np.mean(sp_expr) < 0.05 * np.mean(sc_expr) or np.std(sp_expr) < 0.1 * np.std(sc_expr):
            problem_genes.add(gene)
        
        # Use Kolmogorov-Smirnov test to check if distributions are significantly different
        # This helps find genes that behave differently, even if means are similar.
        if len(sc_expr[sc_expr > 0]) > 10 and len(sp_expr[sp_expr > 0]) > 10:
            ks_stat, p_value = ks_2samp(sc_expr[sc_expr > 0], sp_expr[sp_expr > 0])
            if p_value < 0.001:
                problem_genes.add(gene)
    print(f"Identified {len(problem_genes)} problematic genes based on expression distribution.")

# --- Step 2b: Strategic Gene Selection for Training ---
# For robust mapping, we train on highly variable genes among the high-confidence candidates.
# This data-driven approach works well with or without annotations.
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
# The validation method will be chosen based on the presence of annotations.
print("\nPerforming validation...")

if annotations_present:
    print("Using supervised validation based on cell type annotations...")
    # --- 5a. Cell Type Composition Validation ---
    # Compare the frequency of cell types between the two datasets. They should be similar.
    # sp_freq = ad_map.obs['cell_type'].value_counts(normalize=True)
    # sc_freq = adata_sc.obs['cell_type'].value_counts(normalize=True)
    # common_types = list(set(sp_freq.index) & set(sc_freq.index))
    # composition_corr = np.corrcoef(sp_freq[common_types], sc_freq[common_types])[0, 1]
    # print(f"Correlation of cell type frequencies: {composition_corr:.4f}")

    # --- 5b. Gene Expression Pattern Validation ---
    # validation_genes = [g for g in common_genes if g not in problem_genes][:20]
    # sc_means = adata_sc[:, validation_genes].to_df().groupby(adata_sc.obs['cell_type']).mean()
    # projected_means = adata_ge[:, validation_genes].to_df().groupby(ad_map.obs['cell_type']).mean()
    # common_types_val = list(set(sc_means.index) & set(projected_means.index))
    # gene_pattern_corr = sc_means.loc[common_types_val].corrwith(projected_means.loc[common_types_val], axis=0)
    # print(f"Average gene expression pattern correlation: {gene_pattern_corr.mean():.4f}")

else:
    print("Using unsupervised validation based on spatial coherence...")
    # Check if the projected gene expression shows spatial patterns.
    # Good mappings should result in spatially coherent expression for many genes.
    # This requires a function to calculate spatial autocorrelation, e.g., Moran's I.
    # For simplicity, we will check a few known markers if available, or top HVGs.
    # Note: A full implementation requires a spatial analysis library like `squidpy`.
    
    # Example:
    # import squidpy as sq
    # sq.gr.spatial_neighbors(adata_ge, coord_type="generic")
    # sq.gr.spatial_autocorr(adata_ge, mode="moran", genes=training_genes[:50])
    # moran_scores = adata_ge.uns['moranI']['I']
    # print(f"Average Moran's I for top 50 training genes: {np.nanmean(moran_scores):.4f}")
    # A positive score indicates that the mapping preserved or created spatial structure.
    pass # Placeholder for unsupervised validation

# %%
# 6. Save the integrated data
print(f"\nSaving integrated AnnData to {output_path}")
# The `adata_ge` object contains the spatial locations from MERSCOPE with the
# corrected and expanded gene expression matrix from the scRNA-seq reference.
adata_ge.write(output_path)

print("\nTangram integration script finished successfully!")
print(f"The final integrated data is saved at: {output_path}")