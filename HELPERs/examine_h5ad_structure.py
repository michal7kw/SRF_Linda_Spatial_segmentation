import h5py
import numpy as np
import pandas as pd
try:
    import anndata
    HAS_ANNDATA = True
except ImportError:
    HAS_ANNDATA = False

def examine_h5ad_structure(file_path):
    """
    Examine the structure of an H5AD file (HDF5 format for AnnData objects)
    
    Parameters:
    file_path (str): Path to the H5AD file
    """
    print(f"Examining structure of: {file_path}")
    print("=" * 50)
    
    # First, examine the HDF5 structure
    with h5py.File(file_path, 'r') as f:
        print("\nHDF5 File Structure:")
        print("-" * 30)
        
        def print_structure(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f"Dataset: {name}")
                print(f"  Shape: {obj.shape}")
                print(f"  Dtype: {obj.dtype}")
                if len(obj.shape) > 0 and obj.shape[0] < 10:  # Print small arrays
                    print(f"  Values: {obj[:]}")
            elif isinstance(obj, h5py.Group):
                print(f"Group: {name}")
                for key in obj.keys():
                    print(f"  Key: {key}")
        
        f.visititems(print_structure)
    
    # If anndata is available, try to load as AnnData object
    if HAS_ANNDATA:
        try:
            print("\nAnnData Object Information:")
            print("-" * 30)
            adata = anndata.read_h5ad(file_path)
            
            adata
            print(f"Shape: {adata.shape} (obs × var)")
            print(f"Number of observations: {adata.n_obs}")
            print(f"Number of variables: {adata.n_vars}")
            
            print("\nObservations (obs) info:")
            print(f"  Shape: {adata.obs.shape}")
            print(f"  Columns: {list(adata.obs.columns)}")
            
            print("\nVariables (var) info:")
            print(f"  Shape: {adata.var.shape}")
            print(f"  Columns: {list(adata.var.columns)}")
            
            print("\nLayers:")
            for layer in adata.layers:
                print(f"  {layer}: {adata.layers[layer].shape}")
            
            print("\nObsm:")
            for key in adata.obsm:
                print(f"  {key}: {adata.obsm[key].shape}")
            
            print("\nVarm:")
            for key in adata.varm:
                print(f"  {key}: {adata.varm[key].shape}")
            
            print("\nUns:")
            for key in adata.uns:
                print(f"  {key}: {type(adata.uns[key])}")
            
        except Exception as e:
            print(f"\nError loading as AnnData: {str(e)}")
    else:
        print("\nAnnData library not available. Install with: pip install anndata")

if __name__ == "__main__":
    file_path = "MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad"
    examine_h5ad_structure(file_path)