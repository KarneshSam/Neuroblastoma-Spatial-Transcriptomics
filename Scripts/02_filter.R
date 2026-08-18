# 02_filter.R
# Purpose : Load post-QC object, filter cells and genes iteratively
#           based on QC metrics, and save the filtered object.
# Output  : Filtered Seurat object saved as _post_filter.rds
# Requires: sample_name, plot_dir  (from 01_qc.R)
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
# Required to construct file paths for loading and saving
repeat {
  sample_name <- trimws(readline(prompt = "\nEnter sample name (e.g. P1, P1.B): "))
  if (nchar(sample_name) == 0) {
    cat("Sample name cannot be empty.\n"); next
  }
  break
}

# ── Prompt: input RDS file ────────────────────────────────────────
# Should be the _post_qc.rds saved by 01_qc.R
repeat {
  rds_path <- trimws(readline(
    prompt = paste0("\nEnter path to ", sample_name, "_post_qc.rds: ")))
  if (nchar(rds_path) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!file.exists(rds_path)) {
    cat("File not found: '", rds_path, "'\n", sep = ""); next
  }
  break
}

obj <- readRDS(rds_path)
cat("Loaded:", rds_path, "\n")
cat("Cells:", ncol(obj), "| Genes:", nrow(obj), "\n")

# Prompt: plot output directory
# Must match the directory used in 01_qc.R
repeat {
  plot_dir <- trimws(readline(prompt = "\nEnter plot output directory: "))
  if (nchar(plot_dir) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!dir.exists(plot_dir)) {
    cat("Directory not found: '", plot_dir, "' — please check.\n", sep = ""); next
  }
  break
}

# Prompt: filter thresholds
# Press Enter to accept the default value shown in brackets
# Defaults are based on typical spatial transcriptomics thresholds

# Helper function — prompts for a numeric value with a default
ask_numeric <- function(prompt_text, default_val, min_val = 0) {
  repeat {
    input <- trimws(readline(prompt = prompt_text))
    # Empty input — use default
    if (nchar(input) == 0) {
      cat("  Using default:", default_val, "\n")
      return(default_val)
    }
    val <- suppressWarnings(as.numeric(input))
    if (is.na(val) || val < min_val) {
      cat("  Invalid — please enter a number >=", min_val, "\n"); next
    }
    return(val)
  }
}

cat("\n── Filter thresholds ────────────────────────────────────────\n")
cat("Press Enter to use the default value shown in [brackets].\n\n")

# Minimum unique genes per cell
min_features <- ask_numeric(
  prompt_text = "Min unique genes per cell [default: 50]: ",
  default_val = 50,
  min_val     = 0
)

# Minimum unique spatial bins per cell
min_bins <- ask_numeric(
  prompt_text = "Min unique spatial bins per cell [default: 3]: ",
  default_val = 3,
  min_val     = 0
)

# Maximum mitochondrial % per cell
max_mt <- ask_numeric(
  prompt_text = "Max mitochondrial % per cell [default: 25]: ",
  default_val = 25,
  min_val     = 0
)

# Minimum total counts per gene across all cells
min_gene_count <- ask_numeric(
  prompt_text = "Min total counts per gene [default: 50]: ",
  default_val = 50,
  min_val     = 0
)

# Safety cap on number of iterations
max_iter <- as.integer(ask_numeric(
  prompt_text = "Max filter iterations [default: 10]: ",
  default_val = 10,
  min_val     = 1
))

cat("\nRunning with thresholds:\n")
cat(sprintf("  Min features   : %d\n", as.integer(min_features)))
cat(sprintf("  Min bins       : %d\n", as.integer(min_bins)))
cat(sprintf("  Max MT %%       : %.1f\n", max_mt))
cat(sprintf("  Min gene count : %d\n", as.integer(min_gene_count)))
cat(sprintf("  Max iterations : %d\n", max_iter))

#####################
# Filter thresholds
#####################
# Adjust these based on QC plots from 01_qc.R
min_features   <- 50    # Minimum unique genes per cell
min_bins       <- 3     # Minimum unique spatial bins per cell
max_mt         <- 25    # Maximum mitochondrial % per cell
min_gene_count <- 50    # Minimum total counts per gene across all cells
max_iter       <- 10    # Safety cap on iterations

cat("\nIterative filtering —", sample_name, "\n")

for (iter in seq_len(max_iter)) {
  
  # Record counts before round of filtering
  n_cells_before <- ncol(obj)
  n_genes_before <- nrow(obj)

  # Recompute feature count from current matrix
  obj[["nFeature_fresh"]] <- colSums(
    obj[["Spatial.Polygons"]]$counts > 0)

  # Cell filter
  obj <- subset(obj,
                subset = nFeature_fresh  > min_features &
                         num_unique_bins > min_bins     &
                         percent.mt      < max_mt)

  # Gene filter
  keep_genes <- rownames(obj)[
    Matrix::rowSums(obj[["Spatial.Polygons"]]$counts) > min_gene_count]
  obj <- subset(obj, features = keep_genes)
  
  # Record counts after filtering
  n_cells_after <- ncol(obj)
  n_genes_after <- nrow(obj)

  # Report for each round
  cat("Round", iter, "→",
      "Cells:", n_cells_before, "→", n_cells_after,
      "| Genes:", n_genes_before, "→", n_genes_after, "\n")
  
  # Convergence check: 
  # if no more cells or genes are being filtered out, stop
  if (n_cells_before == n_cells_after &
      n_genes_before == n_genes_after) {
    cat("Converged at round", iter, "\n")
    break
  }
  
  # Max iteration check
  # If still changing after max_iter rounds, warn and stop
  if (iter == max_iter) {
    cat("Warning: max iterations reached — check thresholds\n")
    break
  }

  # Warning if losing too many cells per round
  pct_lost <- (n_cells_before - n_cells_after) / n_cells_before * 100
  if (pct_lost > 10) {
    cat("Warning: losing", round(pct_lost, 1),
        "% of cells in round", iter,
        "— consider relaxing thresholds\n")
  }
}

# Recompute new UMI count per cells after filtering 
obj[["nCount_fresh"]] <- colSums(obj[["Spatial.Polygons"]]$counts)

cat("\nFiltering done:", ncol(obj), "cells |", nrow(obj), "genes\n")
cat("Proceed to 03_normalise.R\n")

