library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)

####################
# Load Sample name
####################
# All file paths are constructed from sample_name
# Keep asking until a non-empty string is entered
repeat {
  sample_name <- trimws(readline(prompt = "\nEnter sample name (e.g. P1, P2): "))
  if (nchar(sample_name) == 0) {
    cat("Sample name cannot be empty.\n"); next
  }
  break
}

# Prompt: post_cnv RDS file
# Seurat object saved at the end of 09_cnv_downstream_analysis.R
# Contains cnv_score and infercnv_subclone in meta.data
repeat {
  rds_path <- trimws(readline(
    prompt = paste0("\nEnter path to ", sample_name, "_post_cnv.rds: ")))
  if (nchar(rds_path) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!file.exists(rds_path)) {
    cat("File not found: '", rds_path, "'\n", sep = ""); next
  }
  break
}

obj <- readRDS(rds_path)
cat("Loaded:", rds_path,
    "| Cells:", ncol(obj), "| Genes:", nrow(obj), "\n")

# Prompt: chromosome-level csv file
# chr_wide CSV saved by 09_cnv_downstream_analysis.R
# Rows = subclones, columns = chromosomes, values = dominant CNV state
repeat {
  chr_wide_path <- trimws(readline(
    prompt = paste0("\nEnter path to ", sample_name, "_cnv_chr_wide.csv: ")))
  if (nchar(chr_wide_path) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!file.exists(chr_wide_path)) {
    cat("File not found: '", chr_wide_path, "'\n", sep = ""); next
  }
  break
}

chr_wide <- read.csv(chr_wide_path, check.names = FALSE)
cat("Loaded chr_wide:", nrow(chr_wide), "subclones x",
    ncol(chr_wide) - 1, "chromosomes\n")

# Prompt: arm-level csv file 
# arm_wide CSV saved by 09_cnv_downstream_analysis.R
# Rows = subclones, columns = chromosome arms, values = dominant CNV state
repeat {
  arm_wide_path <- trimws(readline(
    prompt = paste0("\nEnter path to ", sample_name, "_cnv_arm_wide.csv: ")))
  if (nchar(arm_wide_path) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!file.exists(arm_wide_path)) {
    cat("File not found: '", arm_wide_path, "'\n", sep = ""); next
  }
  break
}

arm_wide <- read.csv(arm_wide_path, check.names = FALSE)
cat("Loaded arm_wide:", nrow(arm_wide), "subclones x",
    ncol(arm_wide) - 1, "arms\n")

# Prompt: InferCNV output directory
# Needed to reload cell_groupings (cell barcode -> subclone mapping)
repeat {
  infercnv_dir <- trimws(readline(
    prompt = "\nEnter InferCNV output directory path: "))
  if (nchar(infercnv_dir) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!dir.exists(infercnv_dir)) {
    cat("Directory not found: '", infercnv_dir, "'\n", sep = ""); next
  }
  break
}

# Prompt: plot output directory
repeat {
  plot_dir <- trimws(readline(
    prompt = "\nEnter directory to save plots: "))
  if (nchar(plot_dir) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  break
}

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
  cat("Created plot directory:", plot_dir, "\n")
}

# Prompt: tumour cell type prefix
# Used to filter meta.data to tumour cells for proportion plots
# e.g. "NE" to keep NE1, NE2 ... regardless of how many exist
repeat {
  tumour_prefix <- trimws(readline(
    prompt = "\nEnter tumour cell type prefix (e.g. NE): "))
  if (nchar(tumour_prefix) == 0) {
    cat("Prefix cannot be empty.\n"); next
  }
  break
}

# Detect all tumour subtypes present in the object
# Uses cell_type_hr column — works with any number of subtypes
tumour_types <- sort(unique(obj$cell_type_hr[
  grepl(paste0("^", tumour_prefix), obj$cell_type_hr)]))
cat("Tumour subtypes detected:", paste(tumour_types, collapse = ", "), "\n")

if (length(tumour_types) == 0) {
  stop("No cell types found with prefix '", tumour_prefix, "'")
}

# Load cell groupings
# Maps each cell barcode to its InferCNV subclone assignment
hmm_prefix     <- "17_HMM_predHMMi6.leiden.hmm_mode-subclusters"
groupings_file <- file.path(infercnv_dir,
  paste0(hmm_prefix, ".cell_groupings"))

if (!file.exists(groupings_file)) {
  stop("Cell groupings file not found: ", groupings_file)
}

cell_groupings <- read.table(groupings_file,
                             header    = TRUE,
                             sep       = "\t",
                             col.names = c("subcluster", "cell_barcode"))

cat("Cell groupings loaded:", nrow(cell_groupings), "cells\n")

