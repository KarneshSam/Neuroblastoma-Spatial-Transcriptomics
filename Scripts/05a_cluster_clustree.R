# 05a_cluster_clustree.R
# Purpose : Run Leiden clustering at 5 user-defined resolutions,
#           then plot clustree to visualise cluster stability.
#           Use clustree to decide which resolution to use in 05b.
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
# Loads the object saved at the end of 04_banksy_pca_umap.R
# Contains BANKSY embedding, PCA, UMAP, and kNN graph
in_rds <- file.path("results", "tmp", sample_name,
                    paste0(sample_name, "_post_banksy.rds"))
if (!file.exists(in_rds)) {
  stop("File not found: ", in_rds,
       "\nPlease run 04_banksy_pca_umap.R first.")
}

obj <- readRDS(in_rds)
cat("Loaded:", in_rds, "\n")

# Output directory
tmp_dir <- file.path("results", "tmp", sample_name)

# Set default assay
DefaultAssay(obj) <- "BANKSY.0.2l"

# Prompt: enter 5 resolutions
# Higher resolution = more clusters, lower = fewer
# Keep asking until exactly 5 valid numbers are entered
repeat {
  cat("\nEnter 5 cluster resolutions between 0 and 1.5.\n")
  res_input   <- trimws(readline(prompt = "Comma-separated (e.g. 0.1,0.3,0.5,0.8,1.2): "))
  
  # Split by comma and coerce to numeric
  resolutions <- suppressWarnings(as.numeric(strsplit(res_input, ",")[[1]]))
  
  # Validate: exactly 5 values, no NAs, all within 0-1.5
  if (length(resolutions) == 5 &&      # Requires 5 resolution
      !any(is.na(resolutions))  &&     
      all(resolutions >= 0)     &&     # Greater than zero and lesser than 1.5
      all(resolutions <= 1.5)) break
  cat("Invalid — please enter exactly 5 numbers between 0 and 1.5, comma-separated.\n")
}

cat("\nRunning clustering at:",
    paste(round(resolutions, 2), collapse = ", "), "\n\n")

###################################
# FindClusters at each resolution
####################################
# Uses the kNN graph built in 04 - FindNeighbors
# algorithm = 4 is Leiden, method = "igraph" is required for Leiden
# Each resolution adds a new column to obj@meta.data
for (res in resolutions) {
  col_name <- paste0("banksy_lam0.2_res", res)
  obj <- FindClusters(obj,
                      cluster.name = col_name,
                      resolution   = res,
                      algorithm    = 4,            # Leiden clustering
                      method       = "igraph")
  # Count clusters found at this resolution                    
  n_cl <- nlevels(obj[[col_name, drop = TRUE]])
  cat(sprintf("  res %.2f  ->  %-32s  (%d clusters)\n",
              res, col_name, n_cl))
}

###########
# Clustree
###########
# Visualises how cells move between clusters as resolution increases
# Stable resolutions show clean, non-mixing branches
# Use this plot to pick your preferred resolution for 05b
p_clustree <- clustree(obj, prefix = "banksy_lam0.2_res") +
  ggtitle(paste("Clustree —", sample_name))

ggsave(file.path(tmp_dir, paste0(sample_name, "_clustree.pdf")),
       p_clustree, width = 12, height = 10)
cat("Saved: Clustree plot\n")

# Saves obj with all 5 resolution cluster columns added to meta.data
# Both 05b_plots.R and 06_markers.R load from this file
out_rds <- file.path(tmp_dir, paste0(sample_name, "_post_cluster.rds"))
saveRDS(obj, file = out_rds)
cat("\nSaved:", out_rds, "\n")
cat("Inspect the clustree, then proceed to 05b_plots.R\n")