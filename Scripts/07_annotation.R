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


