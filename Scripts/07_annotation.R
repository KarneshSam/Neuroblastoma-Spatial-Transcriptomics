# Load final object 
obj <- readRDS(paste0(sample_name, "_final_", chosen_res, ".rds"))
cat("Loaded:", paste0(sample_name, "_final_", chosen_res, ".rds"), "\n")
cat("Annotating:", sample_name,
    "| Resolution:", chosen_res,
    "| Clusters:", nlevels(obj[[chosen_res, drop = TRUE]]), "\n")

# Show clusters and top markers 
# Print cluster IDs and top markers to guide annotation decisions
cluster_ids <- levels(obj[[chosen_res, drop = TRUE]])

cat("\nClusters in this object:\n")
print(cluster_ids)

cat("\nTop 5 markers per cluster (from 06_markers.R):\n")
print(obj@misc[["top5markers"]])

#######################
# Cell Type Annotation
#######################
# Prompt: enter cell type label for each cluster
# Initial setup: empty vector to store labels
cell_type_labels <- c()

# Loop through clusters and prompt user for annotation
for (cl in cluster_ids) {
  repeat {
    label <- trimws(readline(prompt = paste0("  Cluster ", cl, " -> label: ")))

    # Empty input — re-ask, a label is required for every cluster
    if (nchar(label) == 0) {
      cat("  Label cannot be empty. Please enter a cell type name.\n")
      next
    }
  }
  # Store label with cluster ID as name
  cell_type_labels <- c(cell_type_labels, label)
  names(cell_type_labels)[length(cell_type_labels)] <- cl
}
