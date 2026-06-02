#!/usr/bin/env Rscript

# Script to generate marker gene expression dotplots for each sample
# Creates plots showing marker genes across cell types using Seurat
# Note: Run convert_h5ad_to_rds.py first to generate CSV files

library(Seurat)
library(ggplot2)
library(dplyr)

# Configuration
RDS_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation/MERGE_SPATIAL_DATA/RDS_FILES"
MARKER_GENES_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation/MERGE_SPATIAL_DATA/3_MARKER_GENES/marker_genes_per_sample_inclusive_new_clusters"
OUTPUT_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation/MERGE_SPATIAL_DATA/8_ADDITIONAL_PLOTS"
CLUSTER_COLUMN <- "leiden_new"
TOP_N_GENES <- 5

# Sample list
SAMPLES <- c(
  "e16_R1_mut",
  "e16_R2_ctrl",
  "p0_R5p_mut",
  "p0_R6_ctrl",
  "p7_R1_mut",
  "p7_R2_ctrl",
  "p30_R3_mut",
  "p30_R4p_ctrl"
)

cat("================================================================================\n")
cat("MARKER GENE DOTPLOT GENERATION (R/Seurat)\n")
cat("================================================================================\n\n")

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", OUTPUT_DIR, "\n\n")

# Check if RDS files exist
if (!dir.exists(RDS_DIR)) {
  cat("ERROR: RDS directory not found:", RDS_DIR, "\n")
  cat("Please run convert_h5ad_to_rds.py first to generate CSV files\n")
  quit(status = 1)
}

# Function to load marker genes
load_marker_genes <- function(sample_name, marker_genes_dir, top_n = 5) {
  folder_name <- paste0("sample_", sample_name)
  csv_path <- file.path(marker_genes_dir, folder_name, "top_30_marker_genes.csv")

  if (!file.exists(csv_path)) {
    cat("Warning: Could not find", csv_path, "\n")
    return(NULL)
  }

  # Read CSV without headers
  df <- read.csv(csv_path, header = FALSE, stringsAsFactors = FALSE)

  # First row contains cluster IDs
  cluster_ids <- as.character(df[1, ])

  # Collect top N genes for each cluster
  all_genes <- c()
  for (col_idx in 1:ncol(df)) {
    if (nrow(df) > 1) {
      # Get genes from rows 2 to (top_n + 1)
      end_row <- min(top_n + 1, nrow(df))
      genes <- df[2:end_row, col_idx]
      genes <- genes[!is.na(genes) & genes != ""]
      all_genes <- c(all_genes, genes)
    }
  }

  # Remove duplicates while preserving order
  unique_genes <- unique(all_genes)

  return(list(genes = unique_genes, clusters = cluster_ids))
}

# Process each sample
for (sample_name in SAMPLES) {
  cat("================================================================================\n")
  cat("Processing:", sample_name, "\n")
  cat("================================================================================\n")

  # Load expression and metadata
  expr_path <- file.path(RDS_DIR, paste0(sample_name, "_expression.csv"))
  meta_path <- file.path(RDS_DIR, paste0(sample_name, "_metadata.csv"))

  if (!file.exists(expr_path) || !file.exists(meta_path)) {
    cat("Warning: Files not found for", sample_name, "\n")
    cat("  Expected:", expr_path, "\n")
    cat("  Expected:", meta_path, "\n")
    next
  }

  cat("Loading expression data...\n")
  expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
  cat("Loading metadata...\n")
  meta <- read.csv(meta_path, row.names = 1, check.names = FALSE)

  # Create Seurat object
  cat("Creating Seurat object...\n")
  sample_obj <- CreateSeuratObject(counts = expr, meta.data = meta)
  cat("Sample cells:", ncol(sample_obj), "\n")

  # Load marker genes
  marker_result <- load_marker_genes(sample_name, MARKER_GENES_DIR, TOP_N_GENES)

  if (is.null(marker_result)) {
    next
  }

  marker_genes <- marker_result$genes
  cat("Loaded", length(marker_genes), "unique marker genes\n")

  # Filter genes that exist in the dataset
  available_genes <- marker_genes[marker_genes %in% rownames(sample_obj)]
  cat("Available genes:", length(available_genes), "\n")

  if (length(available_genes) == 0) {
    cat("Warning: No marker genes found in dataset\n")
    next
  }

  # Set identity to leiden_new cluster
  Idents(sample_obj) <- CLUSTER_COLUMN

  # Calculate plot dimensions
  n_genes <- length(available_genes)
  n_clusters <- length(unique(sample_obj[[CLUSTER_COLUMN]][[1]]))

  plot_width <- max(12, n_clusters * 0.5 + 3)  # Add 3 inches for legend
  plot_height <- max(8, n_genes * 0.25)

  # Create dotplot
  cat("Creating dotplot...\n")

  p <- DotPlot(
    sample_obj,
    features = available_genes,
    cols = c("blue", "red"),
    dot.scale = 12,  # Larger dots
    dot.min = 0,
    cluster.idents = TRUE  # Cluster the identities (adds dendrogram)
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 16, face = "plain"),
      axis.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18, face = "bold"),
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(
      title = paste("Marker Genes -", sample_name),
      x = "Cluster",
      y = "Gene"
    ) +
    coord_flip()  # Flip to have genes on y-axis and clusters on x-axis

  # Save plot
  output_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_marker_genes_dotplot.png"))

  ggsave(
    filename = output_path,
    plot = p,
    width = plot_width,
    height = plot_height,
    dpi = 300,
    units = "in"
  )

  cat("Saved:", output_path, "\n\n")
}

cat("================================================================================\n")
cat("COMPLETED\n")
cat("================================================================================\n")
cat("All plots saved to:", OUTPUT_DIR, "\n")
