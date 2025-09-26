# %% [markdown]
# # MERSCOPE Analysis

# %%
# Import necessary libraries
import scanpy as sc
import anndata as ad
import pandas as pd
import geopandas as gpd # For .parquet file with geometries
import matplotlib.pyplot as plt
import matplotlib.image as mpimg # For displaying images
import seaborn as sns
import numpy as np # For calculations if needed
import os
import squidpy
import warnings
import scipy.sparse

# Suppress FutureWarning messages
warnings.filterwarnings('ignore', category=FutureWarning)

# %%
os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation")

# Define file paths
base_path = 'DATA/p30-E165/R4'
base_path_integrated = 'DATA_INTEGRATION'
h5ad_file = os.path.join(base_path_integrated, 'p30_R4_tangram_integrated.h5ad')
original_h5ad_file = os.path.join(base_path, 'data.h5ad')

# %% [markdown]
# ## 1. Data Loading

# %%
adata_integrated = None
try:
    adata_integrated = sc.read_h5ad(h5ad_file)
    print(f"Successfully loaded Integrated AnnData file: {h5ad_file}")
    print(adata_integrated)
except FileNotFoundError:
    print(f"Integrated AnnData file not found: {h5ad_file}.")
except Exception as e:
    print(f"Error loading Integrated AnnData file {h5ad_file}: {e}.")

# %%
adata_original = None
try:
    adata_original = sc.read_h5ad(original_h5ad_file)
    print(f"Successfully loaded Raw Spatial AnnData file: {original_h5ad_file}")
    print(adata_original)
except FileNotFoundError:
    print(f"Raw spatial AnnData file not found: {original_h5ad_file}.")
except Exception as e:
    print(f"Error loading raw spatial AnnData file {original_h5ad_file}: {e}.")


# %%
print(adata_integrated.obs.index[:5])
print(adata_integrated.obs.volume.head())
print(adata_integrated.obs.center_x.head())
print(adata_integrated.obs.leiden.head())


# %%
print(f"adata_integrated.X data type: {adata_integrated.X.dtype}")
if scipy.sparse.issparse(adata_integrated.X):
    print(adata_integrated.X[:5, :5].toarray())
    if adata_integrated.X.nnz > 0: # nnz is number of stored_elements
        print(f"Min non-zero value in adata_integrated.X: {adata_integrated.X.data.min()}")
        print(f"Max non-zero value in adata_integrated.X: {adata_integrated.X.data.max()}")
else: # Dense array
    print(adata_integrated.X[:5, :5])
    print(f"Min value in adata_integrated.X: {adata_integrated.X.min()}")
    print(f"Max value in adata_integrated.X: {adata_integrated.X.max()}")

# %%
print(f"adata_original.X data type: {adata_original.X.dtype}")
if scipy.sparse.issparse(adata_original.X):
    print(adata_original.X[:5, :5].toarray())
    if adata_original.X.nnz > 0:
        print(f"Min non-zero value in adata_original.X: {adata_original.X.data.min()}")
        print(f"Max non-zero value in adata_original.X: {adata_original.X.data.max()}")
else: # Dense array
    print(adata_original.X[:5, :5])
    print(f"Min value in adata_original.X: {adata_original.X.min()}")
    print(f"Max value in adata_original.X: {adata_original.X.max()}")

# %% [markdown]
# # 2. Exploratory Data Analysis

# %%
print(f"Number of cells: {adata_integrated.n_obs}")
print(f"Number of genes: {adata_integrated.n_vars}")


# %%
coll_to_summary = 'leiden'

if adata_integrated.obs[coll_to_summary].nunique() < 30 and adata_integrated.obs[coll_to_summary].nunique() > 1:
    plt.figure(figsize=(8, max(4, adata_integrated.obs[coll_to_summary].nunique() * 0.3)))
    sns.countplot(y=adata_integrated.obs[coll_to_summary], order = adata_integrated.obs[coll_to_summary].value_counts(dropna=False).index)
    plt.title(f"Cell Counts by {coll_to_summary}")
    plt.xlabel("Number of Cells")
    plt.ylabel(coll_to_summary)
    plt.tight_layout()
    plt.show()

# %% [markdown]
# ## 3. Visualization

# %%
print(adata_integrated.obsm_keys())
print(list(adata_integrated.obs.columns))

# %%
adata_integrated.obsm['spatial'] = adata_integrated.obs[['center_x', 'center_y']].to_numpy()
spatial_coords_available = True
print(adata_integrated.obsm_keys())

# %%
features = ['leiden', 'rna_count_based_density']
for color_by in features:
    plt.figure(figsize=(10, 10))
    sc.pl.spatial(adata_integrated, color=color_by, spot_size=30, show=False, frameon=False)
    plt.title(f"Spatial Plot of Cells (Colored by {color_by if color_by else 'default'})")
    plt.show()

# %%
print(list(adata_integrated.var.columns))
print(list(adata_integrated.uns))

# %%
adata_integrated.var['highly_variable'].head()

# %%
hvg_genes = adata_integrated.var_names[adata_integrated.var['highly_variable']].tolist()
if len(hvg_genes) > 0:
    num_genes_to_plot = min(len(hvg_genes), 3)
    genes_to_plot = hvg_genes[:num_genes_to_plot]
    print(f"Plotting spatial expression for HVGs: {genes_to_plot}")
    sc.pl.spatial(adata_integrated, color=genes_to_plot, spot_size=30, show=True, frameon=False, ncols=num_genes_to_plot)
else:
    print("No highly variable genes found after computation.")


# %%
sc.tl.pca(adata_integrated)

# %%
sc.pp.neighbors(adata_integrated)
sc.tl.umap(adata_integrated)

# %%
sc.pl.umap(adata_integrated, color="leiden", legend_loc='on data')


# %% [markdown]
# ## 4. Plotting specific genes

# %%
adata_integrated.var_names[:20]

# %%
# Top DEGs from Cluster_Mature GC_vs_Rest_up_significant.csv
genes_to_plot_custom = ['ppfia2', 'rfx3', 'dgkh', 'glis3']
existing_genes = [gene for gene in genes_to_plot_custom if gene in adata_integrated.var_names]

if len(existing_genes) > 0:
    print(f"Plotting spatial expression for custom genes: {existing_genes}")
    sc.pl.spatial(adata_integrated, color=existing_genes, spot_size=30, show=True, frameon=False, ncols=2, cmap='Reds')
else:
    print(f"None of the specified genes were found in the AnnData object: {genes_to_plot_custom}")

# %%
# Top DEGs from Cluster_Mature Cluster_CTX Cluster_CTX_vs_Rest_up_significant.csv
genes_to_plot_custom = ["cpne7", "epha6", "hs6st3", "gm10754"]
existing_genes = [gene for gene in genes_to_plot_custom if gene in adata_integrated.var_names]

if len(existing_genes) > 0:
    print(f"Plotting spatial expression for custom genes: {existing_genes}")
    sc.pl.spatial(adata_integrated, color=existing_genes, spot_size=30, show=True, frameon=False, ncols=2, cmap='Reds')
else:
    print(f"None of the specified genes were found in the AnnData object: {genes_to_plot_custom}")
