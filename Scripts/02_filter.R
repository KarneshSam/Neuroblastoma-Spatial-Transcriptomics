# 02_filter.R
# Purpose : Load post-QC object, filter cells and genes iteratively
#           based on QC metrics, and save the filtered object.
# Output  : results/tmp/<sample>/<sample>_post_filter.rds
#           results/tmp/<sample>/<sample>_post_filter_qc.pdf
# Requires: sample_name, tmp_dir  (from 01_qc.R)
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

# Input RDS file 
# Should be the _post_qc.rds saved by 01_qc.R
in_rds <- file.path("results", "tmp", sample_name,
                    paste0(sample_name, "_post_qc.rds"))
if (!file.exists(in_rds)) {
  stop("File not found: ", in_rds,
       "\nPlease run 01_qc.R first.")
}

obj <- readRDS(in_rds)
cat("Loaded:", in_rds, "\n")
cat("Before filtering:", ncol(obj), "cells |", nrow(obj), "genes\n")

# Set default assay
DefaultAssay(obj) <- "Spatial.Polygons"

# plot output directory
# Must match the directory used in 01_qc.R
tmp_dir <- file.path("results", "tmp", sample_name)

if (!dir.exists(tmp_dir)) {
  dir.create(tmp_dir, recursive = TRUE)
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

# Iterative filtering loop
cat("\nIterative filtering —", sample_name, "\n")

for (iter in seq_len(max_iter)) {

  # Record counts before this round of filtering
  n_cells_before <- ncol(obj)
  n_genes_before <- nrow(obj)

  # Recompute feature count from current matrix
  # Necessary because subsetting changes the matrix each round
  obj[["nFeature_fresh"]] <- colSums(
    obj[["Spatial.Polygons"]]$counts > 0)

  # Cell filter 
  # Keep cells passing all three thresholds simultaneously
  obj <- subset(obj,
                subset = nFeature_fresh  > min_features &
                         num_unique_bins > min_bins     &
                         percent.mt      < max_mt)

  # Gene filter
  # Keep genes with total counts above threshold across remaining cells
  keep_genes <- rownames(obj)[
    Matrix::rowSums(obj[["Spatial.Polygons"]]$counts) > min_gene_count]
  obj <- subset(obj, features = keep_genes)

  # Record counts after this round
  n_cells_after <- ncol(obj)
  n_genes_after <- nrow(obj)

  # Report this round 
  cat("Round", iter, "->",
      "Cells:", n_cells_before, "->", n_cells_after,
      "| Genes:", n_genes_before, "->", n_genes_after, "\n")

  # Convergence check 
  # If nothing changed this round, filtering has stabilised
  if (n_cells_before == n_cells_after &
      n_genes_before  == n_genes_after) {
    cat("  Converged at round", iter, "\n")
    break
  }

  # Max iterations check 
  if (iter == max_iter) {
    cat("  Warning: max iterations reached — check thresholds\n")
    break
  }

  # Warning if losing too many cells per round
  pct_lost <- (n_cells_before - n_cells_after) / n_cells_before * 100
  if (pct_lost > 10) {
    cat("  Warning: losing", round(pct_lost, 1),
        "% of cells in round", iter,
        "— consider relaxing thresholds\n")
  }
}

# Recompute new UMI count per cells after filtering 
obj[["nCount_fresh"]] <- colSums(obj[["Spatial.Polygons"]]$counts)

cat("\nFiltering done:", ncol(obj), "cells |", nrow(obj), "genes\n")

# Plot: post-filter QC summary
# Re-plot QC metrics after filtering to confirm thresholds worked
p1 <- VlnPlot(obj, features = "nFeature_Spatial.Polygons",
              pt.size = 0, group.by = "orig.ident") +
  ggtitle("Features (post-filter)") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

p2 <- VlnPlot(obj, features = "nCount_Spatial.Polygons",
              pt.size = 0, group.by = "orig.ident") +
  ggtitle("UMI counts (post-filter)") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

p3 <- VlnPlot(obj, features = "percent.mt",
              pt.size = 0, group.by = "orig.ident") +
  ggtitle("Mitochondrial % (post-filter)") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

p_filter <- (p1 | p2 | p3) +
  plot_annotation(
    title = paste("Post-filter QC metrics —", sample_name),
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5,
                                            face = "bold")))

filter_plot_path <- file.path(tmp_dir,
  paste0(sample_name, "_post_filter_qc.pdf"))
ggsave(filter_plot_path, p_filter, width = 14, height = 5)
cat("Saved:", filter_plot_path, "\n")

# Save filtered object
# 03_normalise.R loads from this file
out_rds <- file.path(tmp_dir, paste0(sample_name, "_post_filter.rds"))
saveRDS(obj, file = out_rds)
cat("Saved:", out_rds, "\n")

cat("Proceed to 03_normalise.R\n")