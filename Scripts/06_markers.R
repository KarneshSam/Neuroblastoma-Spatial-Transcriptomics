# 06_markers.R
# Purpose : Find marker genes for each cluster at the chosen
#           resolution, extract top 5 per cluster, plot DotPlot,
#           and save marker tables as CSV files.
# Output  : results/tmp/<sample>/<sample>_final_<res>.rds
#           results/tmp/<sample>/markers_<sample>_<res>.csv
#           results/tmp/<sample>/top10markers_<sample>_<res>.csv
# Requires: sample_name  (from 00_setup.R)
# Note    : Can be run at a different resolution than 05b if needed
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
  chosen_res <- trimws(readline(prompt = "\nEnter the resolution column name to use for markers: "))
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

############
# DotPlot 
############
# Dot size = fraction of cells expressing the gene
# Dot colour = average expression level
# Rows = clusters, columns = top marker genes
p_dot <- DotPlot(obj,
                 features = unique(obj@misc[["top10markers"]]$gene),
                 group.by = chosen_res) +
  RotatedAxis() +
  labs(title = paste("Top 10 markers —", sample_name, "|", chosen_res)) +
  theme(plot.title = element_text(size = 13, hjust = 0.5,
                                  face = "bold"))
 
ggsave(file.path(tmp_dir,
         paste0(sample_name, "_dotplot_", chosen_res, ".pdf")),
       p_dot, width = 14, height = 7)
cat("Saved: DotPlot\n")

# Save marker tables as CSV files
# Full marker table — all significant markers for all clusters
write.csv(obj@misc[["markers"]],
          file.path(tmp_dir,
            paste0("markers_", sample_name, "_", chosen_res, ".csv")),
          row.names = TRUE)

# Top 10 marker table — condensed summary
write.csv(obj@misc[["top10markers"]],
          file.path(tmp_dir,
            paste0("top10markers_", sample_name, "_", chosen_res, ".csv")),
          row.names = TRUE)

cat("\nDone. Saved marker tables for", sample_name,
    "at resolution", chosen_res, "\n")

#################################
# Drop Unused Resolution Columns
#################################
# Removes all other resolution columns to keep the object clean
# Only the chosen resolution is kept for downstream use
cols_to_remove <- res_cols[res_cols != chosen_res]

obj@meta.data <- obj@meta.data[
  , !colnames(obj@meta.data) %in% cols_to_remove, drop = FALSE]

cat("Removed columns:", paste(cols_to_remove, collapse = ", "), "\n")
cat("Kept:", chosen_res, "\n")

# Save the final object
# Named with sample name and chosen resolution for full traceability
# 07_annotation.R loads from this file
out_rds <- file.path(tmp_dir, paste0(sample_name, "_final_", chosen_res, ".rds"))
saveRDS(obj, file = out_rds)
cat("Saved:", out_rds, "\n")
cat("Proceed to 07_annotation.R\n")
