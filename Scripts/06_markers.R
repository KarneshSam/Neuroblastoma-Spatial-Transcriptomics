# 06_markers.R
# Purpose : Find marker genes for each cluster at the chosen
#           resolution, extract top 5 per cluster, plot DotPlot,
#           and save marker tables as CSV files.
# Requires: sample_name  (from 00_setup.R)
# Note    : Can be run at a different resolution than 05b if needed
# ─────────────────────────────────────────────────────────────────

# Load the clustered dataset (contains different clusters)
obj <- readRDS(paste0(sample_name, "_post_cluster.rds"))
cat("Loaded:", paste0(sample_name, "_post_cluster.rds"), "\n")

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

# Set default assay and identtity class
# FindAllMarkers must use the raw/normalised gene expression assay
# not the BANKSY assay, so we switch back to Spatial.Polygons
DefaultAssay(obj) <- "Spatial.Polygons"
Idents(obj) <- chosen_res

################
# Find Markers
################
# Finds genes upregulated in each cluster vs all other clusters
# only.pos        = TRUE  : only report upregulated markers
# min.pct         = 0.25  : gene must be detected in >25% of cells
# logfc.threshold = 0.25  : minimum log2 fold change
obj@misc[["markers"]] <- FindAllMarkers(
  obj,
  assay           = "Spatial.Polygons",
  group.by        = chosen_res,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25)

#################
# Top 10 Markers
#################
# Filtered to significant markers only (adjusted p < 0.05)
# Then top 10 by average log2 fold change per cluster
obj@misc[["top10markers"]] <- obj@misc[["markers"]] %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  as.data.frame()

# Print top markers
print(obj@misc[["top10markers"]])


