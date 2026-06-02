#!/usr/bin/env python3
"""
Convert h5ad file to RDS format for R/Seurat analysis
Exports each sample separately as RDS files
"""

import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np

# Set project root
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Configuration
MERGED_H5AD_PATH = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
OUTPUT_DIR = 'MERGE_SPATIAL_DATA/RDS_FILES'
CLUSTER_COLUMN = 'leiden_new'

# Sample list
SAMPLES = [
    'e16_R1_mut',
    'e16_R2_ctrl',
    'p0_R5p_mut',
    'p0_R6_ctrl',
    'p7_R1_mut',
    'p7_R2_ctrl',
    'p30_R3_mut',
    'p30_R4p_ctrl'
]

print("="*80)
print("Converting h5ad to RDS for R/Seurat")
print("="*80)
print()

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Output directory: {OUTPUT_DIR}\n")

# Load merged spatial data
print(f"Loading merged data from: {MERGED_H5AD_PATH}")
adata = sc.read_h5ad(MERGED_H5AD_PATH)
print(f"Loaded data shape: {adata.shape}\n")

# Export each sample
for sample_name in SAMPLES:
    print(f"Processing: {sample_name}")

    # Subset data
    sample_adata = adata[adata.obs['sample_id'] == sample_name].copy()
    print(f"  Cells: {sample_adata.n_obs}")

    # Export expression matrix (genes x cells)
    expr_path = os.path.join(OUTPUT_DIR, f'{sample_name}_expression.csv')
    expr_df = pd.DataFrame(
        sample_adata.X.toarray() if hasattr(sample_adata.X, 'toarray') else sample_adata.X,
        index=sample_adata.obs_names,
        columns=sample_adata.var_names
    ).T  # Transpose to genes x cells
    expr_df.to_csv(expr_path)
    print(f"  Saved expression: {expr_path}")

    # Export metadata
    meta_path = os.path.join(OUTPUT_DIR, f'{sample_name}_metadata.csv')
    sample_adata.obs.to_csv(meta_path)
    print(f"  Saved metadata: {meta_path}")

    print()

print("="*80)
print("COMPLETED")
print("="*80)
print(f"Files saved to: {OUTPUT_DIR}")
print("\nYou can now load these in R using:")
print("  expr <- read.csv('sample_expression.csv', row.names=1)")
print("  meta <- read.csv('sample_metadata.csv', row.names=1)")
print("  seurat_obj <- CreateSeuratObject(counts=expr, meta.data=meta)")
