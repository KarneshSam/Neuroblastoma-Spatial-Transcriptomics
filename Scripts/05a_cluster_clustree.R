# 05a_cluster_clustree.R
# Purpose : Run Leiden clustering at 5 user-defined resolutions,
#           then plot clustree to visualise cluster stability.
#           Use clustree to decide which resolution to use in 05b.
# Requires: sample_name  (from 00_setup.R)
# ─────────────────────────────────────────────────────────────────

# Load the Pre-processed seurat object
# This should contain the BANKSY assay, PCA, UMAP reductions and kNN graph
obj <- readRDS(paste0(sample_name, "_post_banksy.rds"))
cat("Loaded:", paste0(sample_name, "_post_banksy.rds"), "\n")

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
clustree(obj, prefix = "banksy_lam0.2_res") +
  ggtitle(paste("Clustree —", sample_name))

# Saves obj with all 5 resolution cluster columns added to meta.data
# Both 05b_plots.R and 06_markers.R load from this file
saveRDS(obj, file = paste0(sample_name, "_post_cluster.rds"))
cat("\nSaved:", paste0(sample_name, "_post_cluster.rds"), "\n")
cat("Inspect the clustree, then proceed to 05b_plots.R\n")