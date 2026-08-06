library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(pheatmap)
library(ggplot2)

##################################
# Load the list of Seurat object
##################################
rds_obj <- readRDS("seu_enact_seg_unprocessed.rds")

# Handle single object or list of any length
if (inherits(rds_obj, "Seurat")) {
  cat("Detected: single Seurat object\n")
  sample_list <- list(S1 = rds_obj)

} else if (is.list(rds_obj) &&
           all(sapply(rds_obj, inherits, "Seurat"))) {
  cat("Detected:", length(rds_obj), "Seurat objects in list\n")
  if (is.null(names(rds_obj))) {
    names(rds_obj) <- paste0("S", seq_along(rds_obj))
    cat("No names found — assigned:",
        paste(names(rds_obj), collapse = ", "), "\n")
  }
  sample_list <- rds_obj

} else {
  stop("RDS does not contain a Seurat object or a list of Seurat objects.")
}

# Print summary of all samples
cat("\nSamples available:\n")
for (nm in names(sample_list)) {
  cat(sprintf("  %-6s | Cells: %d | Genes: %d\n",
              nm, ncol(sample_list[[nm]]), nrow(sample_list[[nm]])))
}

# Prompt: choose sample 
repeat {
  cat("\nAvailable samples:", paste(names(sample_list), collapse = ", "), "\n")
  sample_name <- trimws(readline(prompt = "Which sample to analyse? "))
  if (sample_name %in% names(sample_list)) break
  cat("Invalid '", sample_name,
      "' — please enter one of the listed names.\n", sep = "")
}

