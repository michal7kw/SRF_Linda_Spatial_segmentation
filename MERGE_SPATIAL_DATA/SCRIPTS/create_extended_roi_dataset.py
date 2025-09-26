#!/usr/bin/env python3

import os
import sys
import scanpy as sc
import pandas as pd
import anndata
import numpy as np

# Set project root and add it to the system path
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

print("="*80)
print("CREATING EXTENDED ROI DATASET")
print("Extending ROI selection by adding cells from small clusters")
print("="*80)

# Define the data structure (same as merge scripts)
data_map = {
    'p0-p7': {
        'R1': {'time_point': 'p7', 'region': 'R1', 'type': 'mut'},
        'R2': {'time_point': 'p7', 'region': 'R2', 'type': 'ctrl'},
        'R4': {'time_point': 'p30', 'region': 'R4p', 'type': 'ctrl'},
        'R5': {'time_point': 'p0', 'region': 'R5p', 'type': 'mut'},
        'R6': {'time_point': 'p0', 'region': 'R6', 'type': 'ctrl'}
    },
    'p30-E165': {
        'R1': {'time_point': 'e16', 'region': 'R1', 'type': 'mut'},
        'R2': {'time_point': 'e16', 'region': 'R2', 'type': 'ctrl'},
        'R3': {'time_point': 'p30', 'region': 'R3', 'type': 'mut'}
    }
}

# Configuration
MIN_CLUSTER_SIZE = 50
base_data_dir = 'DATA'

# Step 1: Load the ROI-selected merged data
roi_data_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_all_cells.h5ad'

if not os.path.exists(roi_data_path):
    print(f"ERROR: ROI merged data not found at {roi_data_path}")
    print("Please run merge_hippo_extended.py first with KEEP_ALL_ROI_CELLS = True")
    sys.exit(1)

print(f"Loading ROI-selected merged data from: {roi_data_path}")
roi_adata = sc.read_h5ad(roi_data_path)
print(f"ROI data shape: {roi_adata.shape}")

# Step 2: Check if clustering has been performed
if 'leiden' not in roi_adata.obs.columns and 'leiden_new' not in roi_adata.obs.columns:
    print("ERROR: No leiden clustering found in ROI data.")
    print("Clustering should be performed during the merge step.")
    sys.exit(1)

# Use leiden_new if available, otherwise leiden
cluster_col = 'leiden_new' if 'leiden_new' in roi_adata.obs.columns else 'leiden'
print(f"Using clustering from column: {cluster_col}")

# Step 3: Analyze cluster sizes PER SAMPLE in ROI selection
print(f"\nAnalyzing cluster sizes PER SAMPLE in ROI selection (threshold: {MIN_CLUSTER_SIZE} cells)...")

# Find clusters that are small in ANY sample
clusters_needing_extension = set()
sample_cluster_analysis = {}

for sample_id in roi_adata.obs['sample_id'].unique():
    sample_data = roi_adata[roi_adata.obs['sample_id'] == sample_id]
    sample_cluster_counts = sample_data.obs[cluster_col].value_counts()
    
    small_in_sample = sample_cluster_counts[sample_cluster_counts < MIN_CLUSTER_SIZE].index.tolist()
    large_in_sample = sample_cluster_counts[sample_cluster_counts >= MIN_CLUSTER_SIZE].index.tolist()
    
    clusters_needing_extension.update(small_in_sample)
    
    sample_cluster_analysis[sample_id] = {
        'small_clusters': small_in_sample,
        'large_clusters': large_in_sample,
        'cluster_counts': sample_cluster_counts.to_dict()
    }
    
    print(f"\nSample {sample_id}:")
    for cluster_id in sorted(sample_cluster_counts.index, key=int):
        count = sample_cluster_counts[cluster_id]
        status = "SMALL" if count < MIN_CLUSTER_SIZE else "LARGE"
        print(f"  Cluster {cluster_id}: {count} cells ({status})")

small_clusters = list(clusters_needing_extension)
print(f"\nClusters requiring extension (< {MIN_CLUSTER_SIZE} cells in ANY sample): {len(small_clusters)}")
print(f"Small clusters: {sorted(small_clusters, key=int)}")

if len(small_clusters) == 0:
    print("No clusters require extension. All clusters have >= 50 cells in all samples.")
    print("Copying ROI data as extended dataset...")
    extended_output_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_extended_rois.h5ad'
    roi_adata.write(extended_output_path)
    print(f"Extended dataset saved to: {extended_output_path}")
    sys.exit(0)

# Step 4: Load original data to get additional cells from small clusters
print(f"\nLoading original data to extract cells from small clusters...")
additional_cells_list = []
total_additional_cells = 0

for time_group, regions in data_map.items():
    for region, metadata in regions.items():
        print(f"\nProcessing: {time_group}/{region}")
        
        # Load original full data
        original_path = os.path.join(base_data_dir, time_group, region, 'data.h5ad')
        if not os.path.exists(original_path):
            print(f"  WARNING: Original data not found at {original_path}, skipping...")
            continue
            
        print(f"  Loading original data from: {original_path}")
        original_adata = sc.read_h5ad(original_path)
        
        # Add metadata to match ROI data structure
        sample_id = f"{metadata['time_point']}_{metadata['region']}_{metadata['type']}"
        original_adata.obs['time_point'] = metadata['time_point']
        original_adata.obs['region'] = metadata['region']
        original_adata.obs['type'] = metadata['type']
        original_adata.obs['sample_id'] = sample_id
        
        print(f"  Original data shape: {original_adata.shape}")
        
        # Check if clustering information exists in original data
        if cluster_col not in original_adata.obs.columns:
            print(f"  WARNING: No {cluster_col} clustering in original data, skipping...")
            continue
            
        # Extract cells from small clusters that are NOT already in ROI selection
        original_small_cluster_cells = original_adata[original_adata.obs[cluster_col].isin(small_clusters)]
        
        if original_small_cluster_cells.n_obs == 0:
            print(f"  No cells from small clusters found in original data")
            continue
            
        # Remove cells that are already in ROI selection (avoid duplicates)
        roi_cell_ids = set(roi_adata.obs.index)
        original_cell_ids = set(original_small_cluster_cells.obs.index)
        new_cell_ids = original_cell_ids - roi_cell_ids
        
        if len(new_cell_ids) == 0:
            print(f"  All small cluster cells already in ROI selection")
            continue
            
        # Get only the new cells
        new_cells_mask = original_small_cluster_cells.obs.index.isin(list(new_cell_ids))
        additional_cells = original_small_cluster_cells[new_cells_mask].copy()
        
        print(f"  Found {additional_cells.n_obs} additional cells from small clusters")
        
        # Add roi_assignment for consistency
        additional_cells.obs['roi_assignment'] = 'extended_from_small_cluster'
        
        if additional_cells.n_obs > 0:
            additional_cells_list.append(additional_cells)
            total_additional_cells += additional_cells.n_obs

print(f"\nSUMMARY:")
print(f"  Original ROI selection: {roi_adata.n_obs:,} cells")
print(f"  Additional cells from small clusters: {total_additional_cells:,} cells")
print(f"  Total extended dataset: {roi_adata.n_obs + total_additional_cells:,} cells")

# Step 5: Merge ROI data with additional cells
if additional_cells_list:
    print(f"\nMerging ROI data with {len(additional_cells_list)} additional cell groups...")
    
    # Combine all additional cells
    all_additional_cells = anndata.concat(additional_cells_list, join='outer', index_unique='-')
    print(f"Combined additional cells shape: {all_additional_cells.shape}")
    
    # Merge with ROI data
    print("Merging ROI data with additional cells...")
    extended_adata = anndata.concat([roi_adata, all_additional_cells], join='outer', index_unique='-')
    
    print(f"Extended data shape: {extended_adata.shape}")
    
else:
    print("No additional cells to add.")
    extended_adata = roi_adata.copy()

# Step 6: Verify cluster extension
print(f"\nVerifying cluster extension...")
extended_cluster_counts = extended_adata.obs[cluster_col].value_counts().sort_index()
original_total_counts = roi_adata.obs[cluster_col].value_counts().sort_index()

print(f"\nOverall cluster sizes after extension:")
for cluster_id in sorted(extended_cluster_counts.index, key=int):
    original_count = original_total_counts.get(cluster_id, 0)
    extended_count = extended_cluster_counts[cluster_id]
    added = extended_count - original_count
    status = "EXTENDED" if str(cluster_id) in [str(c) for c in small_clusters] else "UNCHANGED"
    print(f"  Cluster {cluster_id}: {original_count} → {extended_count} cells (+{added}) [{status}]")

# Show per-sample impact
print(f"\nPer-sample extension summary:")
for sample_id, analysis in sample_cluster_analysis.items():
    print(f"\nSample {sample_id}:")
    small_in_sample = analysis['small_clusters']
    if small_in_sample:
        print(f"  Small clusters extended: {sorted(small_in_sample, key=int)}")
        for cluster_id in small_in_sample:
            original_count = analysis['cluster_counts'][cluster_id]
            # Count how many cells this cluster has in extended dataset
            cluster_cells_extended = extended_adata[extended_adata.obs[cluster_col] == cluster_id]
            total_extended_count = cluster_cells_extended.n_obs
            print(f"    Cluster {cluster_id}: {original_count} → {total_extended_count} cells available for analysis")
    else:
        print(f"  No small clusters requiring extension")

# Step 7: Save extended dataset
extended_output_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_extended.h5ad'
print(f"\nSaving extended ROI dataset to: {extended_output_path}")
extended_adata.write(extended_output_path)

# Step 8: Final summary
print("\n" + "="*80)
print("EXTENSION COMPLETE")
print("="*80)
print(f"Original ROI cells: {roi_adata.n_obs:,}")
print(f"Additional cells added: {total_additional_cells:,}")
print(f"Final extended dataset: {extended_adata.n_obs:,} cells")
print(f"Clusters extended: {len(small_clusters)}")
print(f"Extended dataset saved to: {extended_output_path}")
print("\nThis extended dataset should now be used for all subsequent analyses.")
print("="*80)