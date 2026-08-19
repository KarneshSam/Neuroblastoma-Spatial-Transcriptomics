# 05b_plots.R
# Purpose : Generate UMAP, spatial, and SNN heatmap plots for the
#           resolution chosen after inspecting clustree in 05a.
# Output  : results/tmp/<sample>/plots (PDF files)
# Requires: sample_name  (from 00_setup.R)
# ─────────────────────────────────────────────────────────────────

library(Seurat)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)
library(leiden)
library(Banksy)
library(pheatmap)
library(SeuratWrappers)
library(clustree)
library(Matrix)

# Prompt: sample name
repeat {
  sample_name <- trimws(readline(prompt = "\nEnter sample name (e.g. P1, P1.B): "))
  if (nchar(sample_name) == 0) {
    cat("Sample name cannot be empty.\n"); next
  }
  break
}
 
# Load checkpoint
# Loads the object saved at the end of 05a_cluster_clustree.R
in_rds <- file.path("results", "tmp", sample_name,
                    paste0(sample_name, "_post_cluster.rds"))
if (!file.exists(in_rds)) {
  stop("File not found: ", in_rds,
       "\nPlease run 05a_cluster_clustree.R first.")
}
 
obj <- readRDS(in_rds)
cat("Loaded:", in_rds, "\n")

# Output directory
tmp_dir <- file.path("results", "tmp", sample_name)

# Set default assay
DefaultAssay(obj) <- "Spatial.Polygons"

# Show available resolutions - for the user
# Finds all resolution columns added by 05a in obj@meta.data
# Prints each one with its cluster count to help the user decide
res_cols <- sort(grep("^banksy_lam0.2_res",
                      colnames(obj@meta.data), value = TRUE))

cat("\nAvailable resolutions for", sample_name, ":\n")
for (i in seq_along(res_cols)) {
  n_cl <- nlevels(obj[[res_cols[i], drop = TRUE]])
  cat(sprintf("  [%d] %-32s  (%d clusters)\n", i, res_cols[i], n_cl))
}

# Prompt: choose resolution 
# Keep asking until a valid column name is entered
repeat {
  chosen_res <- trimws(readline(prompt = "\nEnter the resolution column name to use for plots: "))
  if (chosen_res %in% res_cols) break
  cat("Invalid '", chosen_res, "' — please enter one of the listed column names.\n", sep = "")
}

cat("\nUsing:", chosen_res,
    "| Clusters:", nlevels(obj[[chosen_res, drop = TRUE]]), "\n")

# Make the selected resolution to default identity
Idents(obj) <- chosen_res

#######
# UMAP 
#######
# Cells coloured by cluster, labels shown at cluster centroids
p_umap <- DimPlot(obj,
                  reduction = "umap.banksy.0.2",
                  group.by  = chosen_res,
                  label     = TRUE) +
  ggtitle(paste("UMAP —", sample_name, "|", chosen_res))
 
ggsave(file.path(tmp_dir, paste0(sample_name, "_umap_", chosen_res, ".pdf")),
       p_umap, width = 8, height = 7)
cat("Saved: UMAP plot\n")

###################
# Spatial overview 
###################
# Shows spatial distribution of clusters
# Useful to see if clusters are spatially coherent or mixed
p_spatial <- ImageDimPlot(obj, group.by = chosen_res) +
  ggtitle(paste("Spatial —", sample_name, "|", chosen_res))
 
ggsave(file.path(tmp_dir, paste0(sample_name, "_spatial_", chosen_res, ".pdf")),
       p_spatial, width = 9, height = 7)
cat("Saved: Spatial overview plot\n")

###########################
# Spatial split by cluster 
###########################
# Shows spatial distribution of each cluster separately
# Useful to see if clusters have distinct spatial patterns
pdf(file.path(tmp_dir, paste0(sample_name, "_spatial_split_", chosen_res, ".pdf")),
    width = 14, height = 10)

ImageDimPlot(obj,
             group.by        = chosen_res,
             split.by        = chosen_res,
             dark.background = TRUE,
             size            = 0.3)
dev.off()
cat("Saved: Spatial split plot\n")

# SNN connectivity heatmap
# Shows average SNN edge weight between every pair of clusters
# High values = clusters share many neighbours = potentially similar
# Helps identify clusters that may be over-split
banksy_cells    <- Cells(obj[["BANKSY.0.2l"]])                   # Extract the cells
clusters_banksy <- obj[[chosen_res]][banksy_cells, 1]            # Cluster assignments for those cells
banksy_snn_mat  <- obj@graphs[["BANKSY.0.2l_snn"]]               # Extract the SNN graph
b_levels        <- as.character(sort(unique(clusters_banksy)))   # Cluster levels
n_cl            <- length(b_levels)

# Initialise empty connectivity matrix
conn_mat <- matrix(0, nrow = n_cl, ncol = n_cl,
                   dimnames = list(b_levels, b_levels))

# Generate matrix for each cluster pairs 
# By having mean value for each cluster pairs
for (ci in b_levels) {
  for (cj in b_levels) {
    cells_i <- banksy_cells[clusters_banksy == ci]
    cells_j <- banksy_cells[clusters_banksy == cj]

    # Extract submatrix for this cluster pair
    sub_mat <- as.matrix(banksy_snn_mat[cells_i, cells_j])
    # Use only non-zeros for mean
    nz      <- sub_mat[sub_mat > 0]
    conn_mat[ci, cj] <- if (length(nz) > 0) mean(nz) else 0
  }
}

# Plot heatmap - colour: white (low) -> steelblue -> red (high) {relationship} 
pdf(file.path(tmp_dir, paste0(sample_name, "_snn_heatmap_", chosen_res, ".pdf")),
    width = 10, height = 8)

pheatmap(conn_mat,
         main            = paste("SNN connectivity —",
                                 sample_name, "|", chosen_res),
         color           = colorRampPalette(
                             c("white", "steelblue", "red"))(100),
         display_numbers = TRUE,
         number_format   = "%.3f")
dev.off()

cat("Saved: SNN connectivity heatmap\n")

cat("\nPlots done. Proceed to 06_markers.R\n")