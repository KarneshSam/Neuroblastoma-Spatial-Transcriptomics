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
                         num_unique_bins > min_bins      &
                         percent.mt      < max_mt)

  # Gene filter
  keep_genes <- rownames(obj)[
    Matrix::rowSums(obj[["Spatial.Polygons"]]$counts) > min_gene_count]
  obj <- subset(obj, features = keep_genes)

  n_cells_after <- ncol(obj)
  n_genes_after <- nrow(obj)

  # Report each round
  cat("Round", iter, "→",
      "Cells:", n_cells_before, "→", n_cells_after,
      "| Genes:", n_genes_before, "→", n_genes_after, "\n")
}
