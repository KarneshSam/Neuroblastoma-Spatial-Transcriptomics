library(Seurat)
library(infercnv)
library(rtracklayer)
library(dplyr)

##################################
# Load the list of Seurat object
##################################
rds_obj <- readRDS("seu_enact_seg_unprocessed.rds")

# Setup: The file can be a single Seurat object or a list of Seurat objects of any length
# Case 1: single Seurat object — wrap in list for uniform handling
# Case 2: named or unnamed list of Seurat objects of any length
if (inherits(rds_obj, "Seurat")) {
  cat("Detected: single Seurat object\n")
  sample_list <- list(S1 = rds_obj)          # assign a name S1 to the single object for downstream use

} 
else if (is.list(rds_obj) &&
           all(sapply(rds_obj, inherits, "Seurat"))) {
  cat("Detected:", length(rds_obj), "Seurat objects in list\n")

  # Use existing names if available, otherwise assign S1, S2, S3 ...
  if (is.null(names(rds_obj))) {
    names(rds_obj) <- paste0("S", seq_along(rds_obj))
    cat("No names found — assigned:",
        paste(names(rds_obj), collapse = ", "), "\n")
  }
  sample_list <- rds_obj               # assign the list to sample_list for downstream use

} 
else {
  stop("The file does not contain a Seurat object or a list of Seurat objects.")
}

# Print summary of all samples
cat("\nSamples available:\n")
for (nm in names(sample_list)) {
  cat(sprintf("  %-6s | Cells: %d | Genes: %d\n",
              nm, ncol(sample_list[[nm]]), nrow(sample_list[[nm]])))
}

#########################
# Prompt: choose sample 
#########################
# Provide the available sample names to the user and ask which one to run InferCNV on
repeat {
  cat("\nAvailable samples:", paste(names(sample_list), collapse = ", "), "\n")
  sample_name <- trimws(readline(prompt = "Which sample to run InferCNV on? "))
  if (sample_name %in% names(sample_list)) break
  cat("Invalid '", sample_name,
      "' — please enter one of the listed names.\n", sep = "")
}

# Load that Seurat object into the obj variable
obj <- sample_list[[sample_name]]
cat("\nLoaded:", sample_name,
    "| Cells:", ncol(obj),
    "| Genes:", nrow(obj), "\n")

# Show available cell types
all_cell_types <- sort(unique(obj$cell_type_hr))
cat("\nCell types available in", sample_name, ":\n")
for (i in seq_along(all_cell_types)) {
  cat(sprintf("  [%d] %s\n", i, all_cell_types[i]))
}

# Prompt: tumour cell types 
# These are the populations InferCNV will test for CNV signal
repeat {
  cat("\nEnter tumour cell type names, comma-separated.\n")
  cat("These will be tested for CNV signal.\n")
  tumour_input <- trimws(readline(prompt = "Tumour types: "))
  tumour_types <- trimws(strsplit(tumour_input, ",")[[1]])
  tumour_types <- tumour_types[nzchar(tumour_types)]

  # If any of the entered names are not found in the object, warn and re-ask 
  invalid <- tumour_types[!tumour_types %in% all_cell_types]
  if (length(tumour_types) == 0) {
    cat("Please enter at least one cell type.\n"); next
  }
  if (length(invalid) > 0) {
    cat("Not found in object:", paste(invalid, collapse = ", "), "\n")
    cat("Please enter names exactly as listed above.\n"); next
  }
  break
}
cat("Tumour types:", paste(tumour_types, collapse = ", "), "\n")