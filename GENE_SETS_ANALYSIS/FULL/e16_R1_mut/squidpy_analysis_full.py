# %% [markdown]
# # Squidpy analysis of Vizgen data

# %% [markdown]
# ## 1. Imports

# %%
import numpy as np
import pandas as pd
from pathlib import Path
import shutil
import os
from matplotlib import pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq
import geopandas as gpd

# %% [markdown]
# ## 2. Setup

# %%
os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation")

# %% [markdown]
# ## 3. Load Data

# %%
data_dir = "./DATA/p30-E165/R1"
# Create a directory to save the results
results_dir = os.path.join(data_dir, "analysis_results")
os.makedirs(results_dir, exist_ok=True)

# %%
adata = sq.read.vizgen(
    path=data_dir,
    counts_file="cell_by_gene.csv",
    meta_file="cell_metadata.csv",
)
print("Data loaded:")
print(adata)

# %% [markdown]
# ### Get Library ID

# %%
# Get the library id from the anndata object
try:
    library_id = list(adata.uns['spatial'].keys())[0]
    print(f"Using library_id: {library_id}")
except (KeyError, IndexError):
    print("Could not automatically determine library_id. Spatial plots may fail.")
    library_id = None

# %% [markdown]
# ### Load Cell Boundaries

# %%
# Load the cell boundaries and add them to the anndata object
try:
    segmentation_path = os.path.join(data_dir, "cell_boundaries.parquet")
    boundaries = gpd.read_parquet(segmentation_path)
    adata.uns['spatial'][library_id]['segmentations'] = boundaries
    print("Successfully loaded cell boundaries.")
except Exception as e:
    print(f"Could not load cell boundaries: {e}")

# %% [markdown]
# ## 4. Pre-processing and QC
# 
# We perform standard pre-processing and quality control steps.

# %%
adata.var_names_make_unique()
sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100, 200, 300), inplace=True)

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
axs[0].set_title("Total transcripts per cell")
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1])
axs[1].set_title("Unique transcripts per cell")
sns.histplot(adata.obs.groupby("fov").sum()["total_counts"], kde=False, ax=axs[2])
axs[2].set_title("Transcripts per FOV")
sns.histplot(adata.obs["volume"], kde=False, ax=axs[3])
axs[3].set_title("Volume of segmented cells")
fig.tight_layout()
plt.savefig(os.path.join(results_dir, "qc_metrics_distribution.png"))
print("Saved QC metrics distribution plot.")

# %% [markdown]
# ### Filtering

# %%
# Filter cells with low expression and genes that are expressed in too few cells.
print(f"Number of cells before filtering: {adata.n_obs}")
sc.pp.filter_cells(adata, min_counts=50)
print(f"Number of cells after filtering by counts: {adata.n_obs}")

print(f"Number of genes before filtering: {adata.n_vars}")
sc.pp.filter_genes(adata, min_cells=10)
print(f"Number of genes after filtering by cells: {adata.n_vars}")

# %% [markdown]
# ### Normalization and Scaling

# %%
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
print("Normalization and scaling complete.")

# %% [markdown]
# ## 5. Dimensionality Reduction and Clustering

# %%
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1.0)
print("Dimensionality reduction and clustering complete.")

# %% [markdown]
# ## 6. Visualization

# %% [markdown]
# ### UMAP

# %%
sc.pl.umap(adata, color=["leiden"], size=10, show=False, save="_leiden.png")
# Move the file to the results directory
Path("figures/umap_leiden.png").rename(f"{results_dir}/umap_leiden.png")
shutil.rmtree("figures")
print("Saved UMAP plot.")

# %% [markdown]
# ### Spatial Scatter

# %%
sq.pl.spatial_scatter(
    adata,
    library_id=library_id,
    color="leiden",
    img=False,
    figsize=(15, 15),
    size=10,
    save="spatial_leiden.png"
)
# Move the file to the results directory
if os.path.exists("figures/spatial_leiden.png"):
    Path("figures/spatial_leiden.png").rename(f"{results_dir}/spatial_leiden.png")
    shutil.rmtree("figures")
    print("Saved spatial scatter plot.")
else:
    print("Could not save spatial scatter plot.")

# %% [markdown]
# ## 7. Spatial Analysis
# Spatial organization of the clustered cells.

# %% [markdown]
# ### Neighborhood Enrichment
# Clusters that are spatially co-enriched.

# %%
sq.gr.spatial_neighbors(adata, coord_type="generic", spatial_key="spatial")
sq.gr.nhood_enrichment(adata, cluster_key="leiden")
sq.pl.nhood_enrichment(
    adata,
    cluster_key="leiden",
    method="average",
    cmap="inferno",
    vmin=-50,
    vmax=100,
    figsize=(7, 7),
    save="_enrichment.png",
)
Path("figures/_enrichment.png").rename(f"{results_dir}/neighborhood_enrichment.png")
shutil.rmtree("figures")
print("Saved neighborhood enrichment plot.")

# %% [markdown]
# ### Spatial Autocorrelation (Moran's I)
# Identify genes that show a non-random spatial distribution.

# %%
sq.gr.spatial_autocorr(adata, mode="moran")
print("Top 10 spatially autocorrelated genes:")
print(adata.uns["moranI"].head(10))

# %% [markdown]
# ### Visualize top spatially autocorrelated genes

# %%
top_autocorr = adata.uns["moranI"]["I"].sort_values(ascending=False).head(4).index.tolist()
sq.pl.spatial_scatter(
    adata,
    library_id=library_id,
    color=top_autocorr,
    size=5,
    cmap="Reds",
    img=False,
    figsize=(10, 10),
    save="spatial_top_autocorr.png"
)
Path("./figures/spatial_top_autocorr.png").rename(f"{results_dir}/spatial_top_autocorr.png")
os.rmdir("figures")
print("Saved top spatially autocorrelated genes plot.")

print("\nAnalysis complete. Results are in the 'analysis_results' directory.")

# %%
# Save the AnnData object, but first remove GeoDataFrame from spatial metadata
# as it cannot be serialized to HDF5 format
print("Preparing to save AnnData object...")

# Store the boundaries separately if they exist
boundaries_to_save = None
if library_id and 'segmentations' in adata.uns['spatial'][library_id]:
    boundaries_to_save = adata.uns['spatial'][library_id]['segmentations']
    # Remove the GeoDataFrame from spatial metadata before saving
    del adata.uns['spatial'][library_id]['segmentations']
    print("Removed GeoDataFrame from spatial metadata for HDF5 compatibility.")

# Save the AnnData object
adata.write(os.path.join(results_dir, "squidpy_adata_1.h5ad"))
print("Successfully saved AnnData object.")

# Save the boundaries separately as parquet file
if boundaries_to_save is not None:
    boundaries_path = os.path.join(results_dir, "cell_boundaries.parquet")
    boundaries_to_save.to_parquet(boundaries_path)
    print(f"Saved cell boundaries separately to: {boundaries_path}")
    
    # Restore the boundaries to the AnnData object for continued use in this session
    adata.uns['spatial'][library_id]['segmentations'] = boundaries_to_save
    print("Restored boundaries to AnnData object for current session.")
