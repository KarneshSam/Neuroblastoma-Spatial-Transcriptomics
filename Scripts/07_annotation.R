# 07_annotation.R
# Purpose : Manually annotate clusters based on marker genes from
#           06_markers.R. Prompts for a cell type label per cluster,
#           handles duplicate labels with confirmation, adds a
#           cell_type column, generates labelled plots, and saves
#           the annotated object.
# Output  : results/annotated_obj/<sample>_annotated.rds
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
                    paste0(sample_name, "_final.rds"))

if (!file.exists(in_rds)) {
  stop("File not found: ", in_rds,
       "\nPlease run 06_markers.R first.")
}

obj <- readRDS(in_rds)
cat("Loaded:", in_rds, "\n")

# Detect chosen resolution from metadata
# The final object has exactly one banksy resolution column
# left after 06_markers.R cleaned up the rest — extract it directly
chosen_res <- grep("^banksy_lam0.2_res",
                   colnames(obj@meta.data), value = TRUE)

if (length(chosen_res) == 0) {
  stop("No banksy resolution column found in meta.data.")
}
if (length(chosen_res) > 1) {
  stop("More than one banksy resolution column found: ",
       paste(chosen_res, collapse = ", "),
       "\nThis should not happen — check 06_markers.R output.")
}

cat("Resolution detected:", chosen_res, "\n")
cat("Annotating:", sample_name,
    "| Clusters:", nlevels(obj[[chosen_res, drop = TRUE]]), "\n")

# Set the default assay
DefaultAssay(obj) <- "Spatial.Polygons"

# Output directory 
# Final annotated object saved to results/annotated_obj/
out_dir <- file.path("results", "annotated_obj")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("Created directory:", out_dir, "\n")
}

# Show clusters and top markers 
# Print cluster IDs and top markers to guide annotation decisions
cluster_ids <- levels(obj[[chosen_res, drop = TRUE]])

cat("\nClusters in this object:\n")
print(cluster_ids)

cat("\nTop 10 markers per cluster (from 06_markers.R):\n")
print(obj@misc[["top10markers"]])

#######################
# Cell Type Annotation
#######################
# Prompt: enter cell type label for each cluster
# Initial setup: empty vector to store labels
cell_type_labels <- c()

cat("\nEnter a cell type label for each cluster.\n")
cat("Note: the same label can be assigned to multiple clusters\n")

# Loop through clusters and prompt user for annotation
for (cl in cluster_ids) {
  repeat {
    label <- trimws(readline(prompt = paste0("  Cluster ", cl, " -> label: ")))

    # Empty input — re-ask, a label is required for every cluster
    if (nchar(label) == 0) {
      cat("  Label cannot be empty. Please enter a cell type name.\n")
      next
    }

    # Duplicate label — warn which cluster already has it and confirm
    # Same cell type in multiple clusters is biologically valid
    # but we ask to confirm it is intentional and not a typo
    if (label %in% cell_type_labels) {
      already_in <- names(cell_type_labels)[cell_type_labels == label]
      cat("  Note: '", label, "' is already assigned to cluster(s): ",
          paste(already_in, collapse = ", "), "", sep = "")
      cat("  Same cell type in multiple clusters is valid — confirm?\n")
      confirm <- trimws(readline(
        prompt = "  Enter 'yes' to confirm, anything else to re-enter: "))
      if (tolower(confirm) == "yes") break
      next
    }

    # Unique valid label — accept and move on
    break
  }

  # Store label with cluster ID as name
  cell_type_labels <- c(cell_type_labels, label)
  names(cell_type_labels)[length(cell_type_labels)] <- cl
}

# Assign annotations to seurat object
# Maps each cluster ID to its cell type label in a new metadata column
obj[["cell_type_hr"]] <- cell_type_labels[
  as.character(obj[[chosen_res, drop = TRUE]])]

# Cell type summary
cat("\nAnnotation summary:\n")
print(table(obj[["cell_type_hr"]]))

#############################
# UMAP labelled by cell type 
#############################
# repel = TRUE prevents overlapping labels on crowded UMAPs
p_umap <- DimPlot(obj,
                  reduction = "umap.banksy.0.2",
                  group.by  = "cell_type_hr",
                  label     = TRUE,
                  repel     = TRUE) +
  ggtitle(paste("UMAP — cell types |", sample_name)) +
  theme(plot.title = element_text(size = 13, hjust = 0.5,
                                  face = "bold"))
 
ggsave(file.path(out_dir,
         paste0(sample_name, "_umap_annotated_", chosen_res, ".pdf")),
       p_umap, width = 9, height = 7)

cat("Saved: UMAP annotated plot\n")

# Spatial plot labelled by cell type 
# Shows annotated cell types in their original tissue coordinates
p_spatial <- ImageDimPlot(obj, group.by = "cell_type_hr") +
  ggtitle(paste("Spatial — cell types |", sample_name)) +
  theme(plot.title = element_text(size = 13, hjust = 0.5,
                                  face = "bold"))
 
ggsave(file.path(out_dir,
         paste0(sample_name, "_spatial_annotated_", chosen_res, ".pdf")),
       p_spatial, width = 9, height = 7)

cat("Saved: Spatial annotated plot\n")

# Spatial split by cell type 
# Each cell type shown in its own panel on the tissue section
pdf(file.path(out_dir,
      paste0(sample_name, "_spatial_split_annotated_", chosen_res, ".pdf")),
    width = 14, height = 10)

ImageDimPlot(obj,
             group.by        = "cell_type_hr",
             split.by        = "cell_type_hr",
             dark.background = TRUE,
             size            = 0.3)
dev.off()

cat("Saved: Spatial split annotated plot\n")

# Save annotated object
# Final object with cell_type column added to meta.data
out_rds <- file.path(out_dir,
                     paste0(sample_name, "_annotated.rds"))
saveRDS(obj, file = out_rds)
cat("Saved:", out_rds, "\n")
 
cat("\nAnnotation complete for", sample_name, "\n")
cat("The final annotated object is ready for 08_infercnv_pipeline.R\n")