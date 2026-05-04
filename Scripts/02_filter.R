# 02_filter.R
# Requires: obj, sample_name  (from 01_qc.R)
# ─────────────────────────────────────────────────────────────────

# Parameters
min_features   <- 50
min_bins       <- 3
max_mt         <- 25
min_gene_count <- 50
max_iter       <- 10

# Filtering based on
# No of Features
# No of bins
# Mitochondrial percentage
# Counts per gene
for (iter in seq_len(max_iter)) {

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
  
  # After filtering
  n_cells_after <- ncol(obj)
  n_genes_after <- nrow(obj)

  # Report each round
  cat("Round", iter, "→",
      "Cells:", n_cells_before, "→", n_cells_after,
      "| Genes:", n_genes_before, "→", n_genes_after, "\n")
  
  # Stop conditions
  if (n_cells_before == n_cells_after &
      n_genes_before == n_genes_after) {
    cat("Converged at round", iter, "\n")
    break
  }
  
  # Max iteration check
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

# Fresh UMI count per cells after filtering 
obj[["nCount_fresh"]] <- colSums(obj[["Spatial.Polygons"]]$counts)

cat("\nFiltering done:", ncol(obj), "cells |", nrow(obj), "genes\n")
cat("Proceed to 03_normalise.R\n")

