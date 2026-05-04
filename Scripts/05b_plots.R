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