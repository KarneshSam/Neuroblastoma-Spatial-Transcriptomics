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
