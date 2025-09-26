# %%
import os
import glob
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import numpy as np
import re

# Set project root and change directory
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
os.chdir(project_root)

# Define paths
# plots_dir = 'MERGE_SPATIAL_DATA/plots_hippo_cluster_highlight_old_rois'
# output_dir = 'MERGE_SPATIAL_DATA/sample_collages_old_rois'

plots_dir = 'MERGE_SPATIAL_DATA/plots_hippo_cluster_highlight_new_rois'
output_dir = 'MERGE_SPATIAL_DATA/sample_collages_new_rois'

os.makedirs(output_dir, exist_ok=True)

# %%
# Get all plot files
plot_files = glob.glob(os.path.join(plots_dir, 'sample_*.png'))

# Extract unique sample names
samples = set()
for plot_file in plot_files:
    filename = os.path.basename(plot_file)
    # Extract sample name from filename like "sample_p0_R3a_mut_cluster_7.png"
    match = re.match(r'sample_(.+)_cluster_\d+\.png', filename)
    if match:
        samples.add(match.group(1))

print(f"Found samples: {sorted(samples)}")

# %%
# Create collage for each sample
for sample in sorted(samples):
    print(f"Creating collage for sample: {sample}")
    
    # Get all plots for this sample
    sample_plots = glob.glob(os.path.join(plots_dir, f'sample_{sample}_cluster_*.png'))
    sample_plots.sort()
    
    if not sample_plots:
        print(f"No plots found for sample {sample}")
        continue
    
    n_plots = len(sample_plots)
    print(f"Found {n_plots} plots for sample {sample}")
    
    # Calculate grid dimensions (try to make it roughly square)
    n_cols = int(np.ceil(np.sqrt(n_plots)))
    n_rows = int(np.ceil(n_plots / n_cols))
    
    # Create figure with subplots
    fig = plt.figure(figsize=(n_cols * 6, n_rows * 3))
    gs = GridSpec(n_rows, n_cols, figure=fig, hspace=0.01, wspace=0.01)
    
    for i, plot_file in enumerate(sample_plots):
        row = i // n_cols
        col = i % n_cols
        
        # Extract cluster number from filename
        cluster_match = re.search(r'cluster_(\d+)\.png', os.path.basename(plot_file))
        cluster_num = cluster_match.group(1) if cluster_match else str(i)
        
        # Load and display image
        ax = fig.add_subplot(gs[row, col])
        img = mpimg.imread(plot_file)
        ax.imshow(img)
        ax.set_title(f'Cluster {cluster_num}', fontsize=12, fontweight='bold')
        ax.axis('off')
    
    # Hide empty subplots
    for i in range(n_plots, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        ax = fig.add_subplot(gs[row, col])
        ax.axis('off')
    
    # Add main title
    fig.suptitle(f'Sample: {sample}', fontsize=16, fontweight='bold', y=0.98)
    
    # Save collage
    collage_filename = os.path.join(output_dir, f'collage_{sample}.png')
    fig.savefig(collage_filename, dpi=150, bbox_inches='tight')
    plt.close(fig)
    
    print(f"Saved collage: {collage_filename}")

print("All collages created successfully!")
# %%