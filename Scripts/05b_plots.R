# 05b_plots.R
# Requires: sample_name  (from 00_setup.R)
# ─────────────────────────────────────────────────────────────────

# Load the clustered dataset
obj <- readRDS(paste0(sample_name, "_post_cluster.rds"))
cat("Loaded:", paste0(sample_name, "_post_cluster.rds"), "\n")

# Show available resolutions - for the user
res_cols <- sort(grep("^banksy_lam0.2_res",
                      colnames(obj@meta.data), value = TRUE))

cat("\nAvailable resolutions for", sample_name, ":\n")
for (i in seq_along(res_cols)) {
  n_cl <- nlevels(obj[[res_cols[i], drop = TRUE]])
  cat(sprintf("  [%d] %-32s  (%d clusters)\n", i, res_cols[i], n_cl))
}

# Prompt: choose resolution 
repeat {
  chosen_res <- trimws(readline(prompt = "\nEnter the resolution column name to use for plots: "))
  if (chosen_res %in% res_cols) break
  cat("Invalid '", chosen_res, "' — please enter one of the listed column names.\n", sep = "")
}

cat("\nUsing:", chosen_res,
    "| Clusters:", nlevels(obj[[chosen_res, drop = TRUE]]), "\n")

# Make the selected resolution to default identity
Idents(obj) <- chosen_res

# UMAP 
DimPlot(obj,
        reduction = "umap.banksy.0.2",
        group.by  = chosen_res,
        label     = TRUE) +
  ggtitle(paste("UMAP —", sample_name, "|", chosen_res))

# Spatial overview 
ImageDimPlot(obj,
             group.by = chosen_res) +
  ggtitle(paste("Spatial —", sample_name, "|", chosen_res))

# Spatial split by cluster 
ImageDimPlot(obj,
             group.by        = chosen_res,
             split.by        = chosen_res,
             dark.background = TRUE,
             size            = 0.3)

# SNN connectivity heatmap
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
