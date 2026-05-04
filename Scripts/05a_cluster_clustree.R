# 05a_cluster_clustree.R
# Requires: sample_name  (from 00_setup.R)
# ─────────────────────────────────────────────────────────────────

# Load the Pre-processed dataset
obj <- readRDS(paste0(sample_name, "_post_banksy.rds"))
cat("Loaded:", paste0(sample_name, "_post_banksy.rds"), "\n")

# Prompt: enter 5 resolutions
repeat {
  cat("\nEnter 5 cluster resolutions between 0 and 1.5.\n")
  res_input   <- trimws(readline(prompt = "Comma-separated (e.g. 0.1,0.3,0.5,0.8,1.2): "))
  resolutions <- suppressWarnings(as.numeric(strsplit(res_input, ",")[[1]]))
  if (length(resolutions) == 5 &&      # Requires 5 resolution
      !any(is.na(resolutions))  &&     
      all(resolutions >= 0)     &&     # Greater than zero and lesser than 1.5
      all(resolutions <= 1.5)) break
  cat("Invalid — please enter exactly 5 numbers between 0 and 1.5, comma-separated.\n")
}

cat("\nRunning clustering at:",
    paste(round(resolutions, 2), collapse = ", "), "\n\n")

# FindClusters at each resolution
for (res in resolutions) {
  col_name <- paste0("banksy_lam0.2_res", res)
  obj <- FindClusters(obj,
                      cluster.name = col_name,
                      resolution   = res,
                      algorithm    = 4,
                      method       = "igraph")
  n_cl <- nlevels(obj[[col_name, drop = TRUE]])
  cat(sprintf("  res %.2f  ->  %-32s  (%d clusters)\n",
              res, col_name, n_cl))
}

# Clustree
clustree(obj, prefix = "banksy_lam0.2_res") +
  ggtitle(paste("Clustree —", sample_name))

# Save with cluster columns added
saveRDS(obj, file = paste0(sample_name, "_post_cluster.rds"))
cat("\nSaved:", paste0(sample_name, "_post_cluster.rds"), "\n")
cat("Inspect the clustree, then proceed to 05b_plots.R\n")