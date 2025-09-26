#!/usr/bin/env python3

import os
import sys
import scanpy as sc
import pandas as pd

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Load one sample to check data structure
h5ad_path = 'DATA/p0-p7/R1/data.h5ad'

print(f"Loading data from: {h5ad_path}")
adata = sc.read_h5ad(h5ad_path)

print(f"\nData shape: {adata.shape}")
print(f"Observations (cells): {adata.n_obs}")
print(f"Variables (genes): {adata.n_vars}")

print("\n=== OBSM KEYS (multidimensional observations) ===")
print(adata.obsm.keys())

print("\n=== OBS COLUMNS (cell metadata) ===")
print(adata.obs.columns.tolist())

print("\n=== VAR COLUMNS (gene metadata) ===")
print(adata.var.columns.tolist())

# Check if spatial coordinates are in obs columns
if 'center_x' in adata.obs.columns and 'center_y' in adata.obs.columns:
    print("\n=== SPATIAL COORDINATES FOUND IN OBS ===")
    print(f"center_x: {adata.obs['center_x'].head()}")
    print(f"center_y: {adata.obs['center_y'].head()}")
    print(f"center_x range: {adata.obs['center_x'].min()} - {adata.obs['center_x'].max()}")
    print(f"center_y range: {adata.obs['center_y'].min()} - {adata.obs['center_y'].max()}")

# Check if there are any spatial-related keys in obsm
for key in adata.obsm.keys():
    print(f"\nobsm['{key}'] shape: {adata.obsm[key].shape}")
    if adata.obsm[key].shape[1] == 2:
        print(f"  Could be spatial coordinates: {adata.obsm[key][:5]}")

print("\n=== SAMPLE OF FIRST FEW ROWS ===")
print(adata.obs.head())