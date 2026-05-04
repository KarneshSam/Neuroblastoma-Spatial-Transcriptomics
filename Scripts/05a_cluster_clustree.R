# 05a_cluster_clustree.R
# Requires: sample_name  (from 00_setup.R)
# ─────────────────────────────────────────────────────────────────

# Load the Pre-processed dataset
obj <- readRDS(paste0(sample_name, "_post_banksy.rds"))
cat("Loaded:", paste0(sample_name, "_post_banksy.rds"), "\n")

