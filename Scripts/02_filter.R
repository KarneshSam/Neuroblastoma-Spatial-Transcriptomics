# 02_filter.R
# Purpose : Filter the Seurat object based on QC metrics.
# Output  : Filtered Seurat object with improved data quality.
# Requires: obj, sample_name  (from 01_qc.R)
# ─────────────────────────────────────────────────────────────────

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

