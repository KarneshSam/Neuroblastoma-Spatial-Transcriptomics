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

obj <- sample_list[[sample_name]]
cat("\nLoaded:", sample_name,
    "| Cells:", ncol(obj),
    "| Genes:", nrow(obj), "\n")

#####################################
# Prompt: InferCNV output directory
#####################################
# Directory containing the InferCNV run outputs
# Should contain run.final.infercnv_obj and HMM output files
repeat {
  infercnv_dir <- trimws(readline(
    prompt = "\nEnter InferCNV output directory path: "))
  if (nchar(infercnv_dir) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!dir.exists(infercnv_dir)) {
    cat("Directory not found: '", infercnv_dir,
        "' — please check the path.\n", sep = ""); next
  }
  break
}

# Prompt: plot output directory
# Where all plots from this script will be saved
repeat {
  plot_dir <- trimws(readline(
    prompt = "\nEnter directory to save plots: "))
  if (nchar(plot_dir) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  break
}

# Create plot directory if it does not exist
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
  cat("Created plot directory:", plot_dir, "\n")
}

##############################################
# Prompt: Tumour & Normal cell type prefixes
##############################################
# Used to filter CNV regions to tumour subclones only
# e.g. "NE" to keep NE1, NE2, NE3 etc.
repeat {
  tumour_prefix <- trimws(readline(
    prompt = "\nEnter tumour subclone prefix to filter CNV regions (e.g. NE): "))
  if (nchar(tumour_prefix) == 0) {
    cat("Prefix cannot be empty.\n"); next
  }
  break
}

# The Infercnv produce subclone names with prefixes for normal cell types too
# These prefixes will be simplified to just the cell type name 
# e.g. "Plasma,Endo,VSMCs,TAM1,TAM2" -> kept the first part alone "."
repeat {
  cat("\nEnter normal cell type prefixes to simplify subclone names.\n")
  cat("Comma-separated (e.g. Plasma,Endo,VSMCs,TAM1,TAM2)\n")
  normal_input   <- trimws(readline(prompt = "Normal prefixes: "))
  normal_prefixes <- trimws(strsplit(normal_input, ",")[[1]])
  normal_prefixes <- normal_prefixes[nzchar(normal_prefixes)]
  if (length(normal_prefixes) == 0) {
    cat("Please enter at least one prefix.\n"); next
  }
  break
}

# Build regex pattern from the normal cell type prefixes
# e.g. c("Plasma","Endo") -> "^(Plasma|Endo)"
normal_pattern <- paste0("^(", paste(normal_prefixes, collapse = "|"), ")")
cat("Normal subclone pattern:", normal_pattern, "\n")

#####################
# InferCNV analysis
#####################
# Load final InferCNV object from the specified directory
# Contains the CNV expression matrix (genes x cells)
infercnv_rds <- file.path(infercnv_dir, "run.final.infercnv_obj")
if (!file.exists(infercnv_rds)) {
  stop("run.final.infercnv_obj not found in: ", infercnv_dir)
}

cat("\nLoading InferCNV object...\n")
infercnv_obj <- readRDS(infercnv_rds)

# CNV score per cell
# Variance of CNV score across genes per cell
# Higher variance = more CNV events = more aneuploid
# Neutral (diploid) cells have CNV values close to 1 = low variance
cnv_mat   <- infercnv_obj@expr.data

cat("CNV matrix dimensions:", nrow(cnv_mat), "genes x",
    ncol(cnv_mat), "cells\n")
# Compute CNV score for each cell based on the variance
cnv_score <- apply(cnv_mat, 2, var)

# Assign CNV scores back to Seurat object
# Only cells in the Infercnv object will have a score
# Cells not in InferCNV (e.g. filtered out) get NA
cells_in_infercnv     <- names(cnv_score)
obj$cnv_score         <- NA
obj$cnv_score[match(cells_in_infercnv, colnames(obj))] <- cnv_score

cat("CNV scores assigned to", sum(!is.na(obj$cnv_score)), "cells\n")


