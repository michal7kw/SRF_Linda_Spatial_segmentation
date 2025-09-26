# %% [markdown]
# # MERSCOPE Analysis

# %% [markdown]
# This notebook performs an analysis of MERSCOPE data for region R2, focusing on data loading, exploratory data analysis (EDA), and visualization.

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

# Suppress FutureWarning messages
warnings.filterwarnings('ignore', category=FutureWarning)

# %%
os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation")

# Define file paths
base_path = 'DATA/p30-E165/R3'
h5ad_file = os.path.join(base_path, 'data.h5ad')

# %% [markdown]
# ## 1. Data Loading
# 
# We will prioritize loading the AnnData file (`.h5ad`) as it should contain the most comprehensive data. If needed, we will supplement with other files.

# %%
cell_by_gene_file = os.path.join(base_path, 'cell_by_gene.csv')
cell_metadata_file = os.path.join(base_path, 'cell_metadata.csv')
differentially_expressed_genes_file = os.path.join(base_path, 'differentially_expressed_genes.csv')
# cell_categories_file = os.path.join(base_path, 'cell_categories.csv')
# cell_numeric_categories_file = os.path.join(base_path, 'cell_numeric_categories.csv')
cell_boundaries_file = os.path.join(base_path, 'cell_boundaries.parquet')
detected_transcripts_file = os.path.join(base_path, 'detected_transcripts.parquet')

summary_image_file = os.path.join(base_path, 'summary.png')

# Attempt to load the .h5ad file
adata = None
try:
    adata = sc.read_h5ad(h5ad_file)
    print(f"Successfully loaded AnnData file: {h5ad_file}")
    print(adata)
except FileNotFoundError:
    print(f"AnnData file not found: {h5ad_file}. Will attempt to load individual files.")
except Exception as e:
    print(f"Error loading AnnData file {h5ad_file}: {e}. Will attempt to load individual files.")

# %%
print(adata.obs.index[:5])
print(adata.obs.volume.head())
print(adata.obs.center_x.head())
print(adata.obs.leiden.head())

# %% [markdown]
# # Load supplementary files if .h5ad loading failed

# %%
cell_categories_df = None
cell_boundaries_gdf = None

if adata is None:
    print("\nAttempting to load data from individual CSV files...")
    try:
        # Load gene expression data
        counts_df = pd.read_csv(cell_by_gene_file, index_col=0) # Assuming first column is cell ID
        print(f"Loaded {cell_by_gene_file}: {counts_df.shape[0]} cells, {counts_df.shape[1]} genes")
        
        # Load cell metadata
        metadata_df = pd.read_csv(cell_metadata_file, index_col=0) # Assuming first column is cell ID
        print(f"Loaded {cell_metadata_file}: {metadata_df.shape[0]} cells, {metadata_df.shape[1]} metadata columns")
        
        # Align indices (important!)
        common_cells = counts_df.index.intersection(metadata_df.index)
        counts_df = counts_df.loc[common_cells]
        metadata_df = metadata_df.loc[common_cells]
        print(f"Found {len(common_cells)} common cells between counts and metadata.")

        if len(common_cells) == 0:
            raise ValueError("No common cells found between cell_by_gene.csv and cell_metadata.csv. Cannot create AnnData object.")

        # Create AnnData object
        adata = ad.AnnData(X=counts_df.values, obs=metadata_df, var=pd.DataFrame(index=counts_df.columns))
        adata.X = adata.X.astype('float32') # Ensure X is float for scanpy operations
        print("Successfully created AnnData object from CSV files.")
        print(adata)
        
    except FileNotFoundError as e:
        print(f"Error: A required CSV file was not found: {e}. Cannot proceed with manual AnnData creation.")
        adata = None # Ensure adata is None if creation fails
    except ValueError as e:
        print(f"Error creating AnnData object: {e}")
        adata = None
    except Exception as e:
        print(f"An unexpected error occurred while loading CSV files: {e}")
        adata = None

# %% [markdown]
# # Load additional data

# %% [markdown]
# ## cell boundaries

# %%
cell_boundaries_gdf = gpd.read_parquet(cell_boundaries_file)
print(f"Loaded {cell_boundaries_file}. Shape: {cell_boundaries_gdf.shape}")

# %%
cell_boundaries_gdf.head()

# %%
cell_boundaries_gdf = cell_boundaries_gdf.set_index('EntityID', drop=False)
cell_boundaries_gdf.head()

# %%
adata.obs.index = adata.obs.index.astype(str)
cell_boundaries_gdf.index = cell_boundaries_gdf.index.astype(str)

print(adata.obs.index[:5])
print(cell_boundaries_gdf.index[:5])

# %%
common_cells_boundaries = adata.obs.index.intersection(cell_boundaries_gdf.index)
common_cells_boundaries[:5]

# %%
adata.uns['cell_boundaries_gdf'] = cell_boundaries_gdf.loc[common_cells_boundaries]

# %% [markdown]
# ## differentially expressed genes

# %%
degs_df = pd.read_csv(differentially_expressed_genes_file, index_col=0) # Assuming first column is cell ID

# %%
degs_df.head()

# %% [markdown]
# ## detected transcripts

# %%
det_trans_df = pd.read_parquet(detected_transcripts_file)

# %%
det_trans_df.head()

# %% [markdown]
# ## 2. Exploratory Data Analysis (EDA)
# 
# Basic statistics and distributions of the data.

# %%
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")

# %%
print("Calculating QC metrics (total_counts, n_genes_by_counts)...")
sc.pp.calculate_qc_metrics(adata, inplace=True)

# %%
fig, axes = plt.subplots(1, 2, figsize=(16, 5))

# Plot for total_counts
sns.histplot(adata.obs['total_counts'], bins=50, kde=False, ax=axes[0])
axes[0].set_xlabel("Total transcripts per cell")
axes[0].set_ylabel("Number of cells")
axes[0].set_title("Distribution of Total Transcripts per Cell")

# Plot for n_genes_by_counts
sns.histplot(adata.obs['n_genes_by_counts'], bins=50, kde=False, ax=axes[1])
axes[1].set_xlabel("Number of genes detected per cell")
axes[1].set_ylabel("Number of cells")
axes[1].set_title("Distribution of Genes Detected per Cell")

plt.tight_layout()
plt.show()

# %%
coll_to_summary = 'leiden'

# print(f"\n--- {coll_to_summary} ---")
# print(adata.obs[coll_to_summary].value_counts(dropna=False))

if adata.obs[coll_to_summary].nunique() < 30 and adata.obs[coll_to_summary].nunique() > 1:
    plt.figure(figsize=(8, max(4, adata.obs[coll_to_summary].nunique() * 0.3)))
    sns.countplot(y=adata.obs[coll_to_summary], order = adata.obs[coll_to_summary].value_counts(dropna=False).index)
    plt.title(f"Cell Counts by {coll_to_summary}")
    plt.xlabel("Number of Cells")
    plt.ylabel(coll_to_summary)
    plt.tight_layout()
    plt.show()

# %% [markdown]
# ## 3. Visualization

# %%
img = mpimg.imread(summary_image_file)
plt.figure(figsize=(10,10))
plt.imshow(img)
plt.axis('off') # Turn off axis numbers and ticks
plt.title("Experiment Summary Image (region_R2/summary.png)")
plt.show()

# %%
print(adata.obsm_keys())
print(list(adata.obs.columns))

# %%
adata.obsm['spatial'] = adata.obs[['center_x', 'center_y']].to_numpy()
spatial_coords_available = True

# %%
adata

# %%
features = ['total_counts', 'n_genes_by_counts', 'leiden']
for color_by in features:
    plt.figure(figsize=(10, 10))
    sc.pl.spatial(adata, color=color_by, spot_size=30, show=False, frameon=False)
    plt.title(f"Spatial Plot of Cells (Colored by {color_by if color_by else 'default'})")
    plt.show()

# %%
print(list(adata.var.columns))
print(list(adata.uns))

# %%
# Check the log1p metadata to understand the issue
print("Checking adata.uns['log1p'] contents:")
if 'log1p' in adata.uns_keys():
    print(adata.uns['log1p'])
    print(f"Keys in log1p: {list(adata.uns['log1p'].keys())}")
else:
    print("No log1p key found in adata.uns")

# %%
# Fix the log1p metadata issue by adding the missing 'base' key or using a workaround
if 'log1p' in adata.uns_keys() and 'base' not in adata.uns['log1p']:
    print("Adding missing 'base' key to log1p metadata...")
    # Set base to None (natural log) which is the default for scanpy
    adata.uns['log1p']['base'] = None
    print("Fixed log1p metadata:", adata.uns['log1p'])

# Alternative approach: use flavor='cell_ranger' which doesn't require log1p metadata
print("Computing highly variable genes with cell_ranger flavor to avoid log1p issues...")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, flavor='cell_ranger')

# %%
import numpy as np
import scipy.sparse

print(f"adata.X data type: {adata.X.dtype}")

if scipy.sparse.issparse(adata.X):
    print(adata.X[:5, :5].toarray())
    if adata.X.nnz > 0: # nnz is number of stored_elements
        print(f"Min non-zero value in adata.X: {adata.X.data.min()}")
        print(f"Max non-zero value in adata.X: {adata.X.data.max()}")
        has_nan_sparse = np.isnan(adata.X.data).any()
        has_inf_sparse = np.isinf(adata.X.data).any()
        print(f"adata.X.data contains NaNs: {has_nan_sparse}")
        print(f"adata.X.data contains Infs: {has_inf_sparse}")

else: # Dense array
    print(adata.X[:5, :5])
    print(f"Min value in adata.X: {adata.X.min()}")
    print(f"Max value in adata.X: {adata.X.max()}")
    has_nan_dense = np.isnan(adata.X).any()
    has_inf_dense = np.isinf(adata.X).any()
    print(f"adata.X contains NaNs: {has_nan_dense}")
    print(f"adata.X contains Infs: {has_inf_dense}")

# %%
hvg_genes = adata.var_names[adata.var['highly_variable']].tolist()
if len(hvg_genes) > 0:
    num_genes_to_plot = min(len(hvg_genes), 3)
    genes_to_plot = hvg_genes[:num_genes_to_plot]
    print(f"Plotting spatial expression for HVGs: {genes_to_plot}")
    sc.pl.spatial(adata, color=genes_to_plot, spot_size=30, show=True, frameon=False, ncols=num_genes_to_plot)
else:
    print("No highly variable genes found after computation.")




