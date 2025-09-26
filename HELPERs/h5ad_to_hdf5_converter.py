#!/usr/bin/env python3
"""
Convert AnnData .h5ad file to HDF5 .hdf5 format

This script loads an AnnData object from a .h5ad file and converts it to 
a standard HDF5 format that can be used with other tools and platforms.
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import h5py
import os
import warnings
from scipy import sparse
import argparse

# Suppress warnings
warnings.filterwarnings('ignore', category=FutureWarning)

def convert_h5ad_to_hdf5(h5ad_file_path, output_hdf5_path=None, compression='gzip'):
    """
    Convert AnnData .h5ad file to HDF5 .hdf5 format
    
    Parameters:
    -----------
    h5ad_file_path : str
        Path to the input .h5ad file
    output_hdf5_path : str, optional
        Path for the output .hdf5 file. If None, will use input filename with .hdf5 extension
    compression : str, default 'gzip'
        Compression method for HDF5 datasets
    """
    
    # Load the AnnData object
    print(f"Loading AnnData file: {h5ad_file_path}")
    try:
        adata = sc.read_h5ad(h5ad_file_path)
        print(f"Successfully loaded AnnData object:")
        print(f"  - Shape: {adata.n_obs} cells × {adata.n_vars} genes")
        print(f"  - Observations (obs): {list(adata.obs.columns)}")
        print(f"  - Variables (var): {list(adata.var.columns)}")
        print(f"  - Unstructured (uns): {list(adata.uns.keys())}")
        print(f"  - Obsm keys: {list(adata.obsm.keys())}")
        print(f"  - Varm keys: {list(adata.varm.keys())}")
        print(f"  - Obsp keys: {list(adata.obsp.keys())}")
    except Exception as e:
        print(f"Error loading AnnData file: {e}")
        return False
    
    # Determine output file path
    if output_hdf5_path is None:
        output_hdf5_path = h5ad_file_path.replace('.h5ad', '.hdf5')
    
    print(f"\nConverting to HDF5 format: {output_hdf5_path}")
    
    try:
        with h5py.File(output_hdf5_path, 'w') as f:
            # Create main groups
            obs_group = f.create_group('obs')
            var_group = f.create_group('var')
            obsm_group = f.create_group('obsm')
            varm_group = f.create_group('varm')
            obsp_group = f.create_group('obsp')
            uns_group = f.create_group('uns')
            
            # Save the main expression matrix (X)
            print("  Saving expression matrix (X)...")
            if sparse.issparse(adata.X):
                # Save sparse matrix in CSR format
                x_group = f.create_group('X')
                x_group.create_dataset('data', data=adata.X.data, compression=compression)
                x_group.create_dataset('indices', data=adata.X.indices, compression=compression)
                x_group.create_dataset('indptr', data=adata.X.indptr, compression=compression)
                x_group.create_dataset('shape', data=adata.X.shape)
                x_group.attrs['format'] = 'csr'
                x_group.attrs['dtype'] = str(adata.X.dtype)
            else:
                # Save dense matrix
                f.create_dataset('X', data=adata.X, compression=compression)
            
            # Save cell/observation metadata (obs)
            print("  Saving observation metadata (obs)...")
            for col in adata.obs.columns:
                data = adata.obs[col]
                if data.dtype == 'object' or data.dtype.name == 'category':
                    # Handle categorical/string data
                    if hasattr(data, 'cat'):
                        # Categorical data
                        obs_group.create_dataset(f'{col}_codes', data=data.cat.codes, compression=compression)
                        obs_group.create_dataset(f'{col}_categories', 
                                               data=[str(cat).encode('utf-8') for cat in data.cat.categories],
                                               compression=compression)
                    else:
                        # String data
                        obs_group.create_dataset(col, 
                                               data=[str(val).encode('utf-8') for val in data],
                                               compression=compression)
                else:
                    # Numeric data
                    obs_group.create_dataset(col, data=data.values, compression=compression)
            
            # Save cell/observation indices
            obs_group.create_dataset('_index', 
                                   data=[str(idx).encode('utf-8') for idx in adata.obs.index],
                                   compression=compression)
            
            # Save gene/variable metadata (var)
            print("  Saving variable metadata (var)...")
            for col in adata.var.columns:
                data = adata.var[col]
                if data.dtype == 'object' or data.dtype.name == 'category':
                    if hasattr(data, 'cat'):
                        var_group.create_dataset(f'{col}_codes', data=data.cat.codes, compression=compression)
                        var_group.create_dataset(f'{col}_categories', 
                                               data=[str(cat).encode('utf-8') for cat in data.cat.categories],
                                               compression=compression)
                    else:
                        var_group.create_dataset(col, 
                                               data=[str(val).encode('utf-8') for val in data],
                                               compression=compression)
                else:
                    var_group.create_dataset(col, data=data.values, compression=compression)
            
            # Save gene/variable indices
            var_group.create_dataset('_index', 
                                   data=[str(idx).encode('utf-8') for idx in adata.var.index],
                                   compression=compression)
            
            # Save obsm (multi-dimensional observations, e.g., PCA, UMAP)
            print("  Saving obsm (multi-dimensional observations)...")
            for key in adata.obsm.keys():
                obsm_group.create_dataset(key, data=adata.obsm[key], compression=compression)
            
            # Save varm (multi-dimensional variables)
            print("  Saving varm (multi-dimensional variables)...")
            for key in adata.varm.keys():
                varm_group.create_dataset(key, data=adata.varm[key], compression=compression)
            
            # Save obsp (pairwise observations, e.g., distances, connectivities)
            print("  Saving obsp (pairwise observations)...")
            for key in adata.obsp.keys():
                if sparse.issparse(adata.obsp[key]):
                    obsp_key_group = obsp_group.create_group(key)
                    obsp_key_group.create_dataset('data', data=adata.obsp[key].data, compression=compression)
                    obsp_key_group.create_dataset('indices', data=adata.obsp[key].indices, compression=compression)
                    obsp_key_group.create_dataset('indptr', data=adata.obsp[key].indptr, compression=compression)
                    obsp_key_group.create_dataset('shape', data=adata.obsp[key].shape)
                    obsp_key_group.attrs['format'] = 'csr'
                else:
                    obsp_group.create_dataset(key, data=adata.obsp[key], compression=compression)
            
            # Save uns (unstructured annotations)
            print("  Saving uns (unstructured annotations)...")
            def save_uns_recursive(group, data, path=""):
                for key, value in data.items():
                    current_path = f"{path}/{key}" if path else key
                    try:
                        if isinstance(value, dict):
                            subgroup = group.create_group(key)
                            save_uns_recursive(subgroup, value, current_path)
                        elif isinstance(value, (list, tuple)):
                            if len(value) > 0 and isinstance(value[0], str):
                                group.create_dataset(key, 
                                                   data=[str(v).encode('utf-8') for v in value],
                                                   compression=compression)
                            else:
                                group.create_dataset(key, data=value, compression=compression)
                        elif isinstance(value, np.ndarray):
                            if value.dtype.kind in ['U', 'S', 'O']:  # String types
                                group.create_dataset(key, 
                                                   data=[str(v).encode('utf-8') for v in value.flatten()],
                                                   compression=compression)
                            else:
                                group.create_dataset(key, data=value, compression=compression)
                        elif isinstance(value, (str, int, float, bool, np.integer, np.floating)):
                            if isinstance(value, str):
                                group.create_dataset(key, data=value.encode('utf-8'))
                            else:
                                group.create_dataset(key, data=value)
                        elif hasattr(value, 'toarray'):  # Sparse matrix
                            sparse_group = group.create_group(key)
                            sparse_group.create_dataset('data', data=value.data, compression=compression)
                            sparse_group.create_dataset('indices', data=value.indices, compression=compression)
                            sparse_group.create_dataset('indptr', data=value.indptr, compression=compression)
                            sparse_group.create_dataset('shape', data=value.shape)
                            sparse_group.attrs['format'] = 'csr'
                        else:
                            print(f"    Warning: Skipping uns['{current_path}'] - unsupported type: {type(value)}")
                    except Exception as e:
                        print(f"    Warning: Could not save uns['{current_path}']: {e}")
            
            save_uns_recursive(uns_group, adata.uns)
            
            # Add metadata about the conversion
            f.attrs['converted_from'] = 'h5ad'
            f.attrs['original_file'] = os.path.basename(h5ad_file_path)
            f.attrs['n_obs'] = adata.n_obs
            f.attrs['n_vars'] = adata.n_vars
            f.attrs['conversion_tool'] = 'h5ad_to_hdf5_converter.py'
            
        print(f"\nSuccessfully converted to HDF5 format!")
        print(f"Output file: {output_hdf5_path}")
        
        # Verify the output file
        file_size = os.path.getsize(output_hdf5_path) / (1024**2)  # Size in MB
        print(f"File size: {file_size:.2f} MB")
        
        return True
        
    except Exception as e:
        print(f"Error during conversion: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Convert AnnData .h5ad file to HDF5 .hdf5 format')
    parser.add_argument('input_file', help='Path to input .h5ad file')
    parser.add_argument('-o', '--output', help='Path to output .hdf5 file (optional)')
    parser.add_argument('-c', '--compression', default='gzip', 
                       choices=['gzip', 'lzf', 'szip'], 
                       help='Compression method (default: gzip)')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found!")
        return 1
    
    # Perform conversion
    success = convert_h5ad_to_hdf5(args.input_file, args.output, args.compression)
    
    return 0 if success else 1

if __name__ == "__main__":
    # Example usage for the specific file mentioned
    # h5ad_file = '202503071102_SESSA-p30-E165_VMSC10702/region_R1/202503071102_SESSA-p30-E165_VMSC10702_region_R1_clustered.h5ad'
    # h5ad_file = '202503071102_SESSA-p30-E165_VMSC10702/region_R2/202503071102_SESSA-p30-E165_VMSC10702_region_R2_clustered.h5ad'
    # h5ad_file = '202503071102_SESSA-p30-E165_VMSC10702/region_R3/202503071102_SESSA-p30-E165_VMSC10702_region_R3_clustered.h5ad'
    # h5ad_file = '202503071102_SESSA-p30-E165_VMSC10702/region_R4/202503071102_SESSA-p30-E165_VMSC10702_region_R4_clustered.h5ad'
    
    # h5ad_file = "202504111150_Sessa-p0-p7_VMSC10702/R4/202504111150_Sessa-p0-p7_VMSC10702_region_R1_clustered.h5ad"
    # h5ad_file = "202504111150_Sessa-p0-p7_VMSC10702/R4/202504111150_Sessa-p0-p7_VMSC10702_region_R2_clustered.h5ad"
    # h5ad_file = "202504111150_Sessa-p0-p7_VMSC10702/R4/202504111150_Sessa-p0-p7_VMSC10702_region_R3_clustered.h5ad"
    h5ad_file = "202504111150_Sessa-p0-p7_VMSC10702/R4/202504111150_Sessa-p0-p7_VMSC10702_region_R4_clustered.h5ad"

    if os.path.exists(h5ad_file):
        print("Converting file...")
        convert_h5ad_to_hdf5(h5ad_file)
    else:
        print(f"File not found: {h5ad_file}")
        print("Please run with the correct file path or use command line arguments:")
        print("python h5ad_to_hdf5_converter.py <input_file.h5ad> [-o output_file.hdf5]") 