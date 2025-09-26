import scanpy as sc
import pandas as pd
import os

# Set working directory
os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation")

# Define file paths
h5ad_path = "DATA_INTEGRATION/p30_R4_tangram_integrated.h5ad"
output_parquet_path = "DATA_INTEGRATION/p30_R4_tangram_integrated.parquet"

# Load the AnnData object
print(f"Loading AnnData from {h5ad_path}...")
adata = sc.read_h5ad(h5ad_path)
print("AnnData loaded.")

# The vizgen tool expects the index to have a name.
# Let's call it 'EntityID' to match conventions.
adata.obs.index.name = 'EntityID'

# Create a pandas DataFrame from the expression data
# The data in .X is float32, which is fine for parquet.
print("Creating DataFrame from AnnData.X...")
df = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)

# Save the DataFrame to a parquet file
print(f"Saving DataFrame to {output_parquet_path}...")
df.to_parquet(output_parquet_path)

print(f"Conversion to parquet complete. File saved at: {output_parquet_path}")