# 04_banksy_pca_umap.R
# Purpose: Run BANKSY spatial embedding, PCA, kNN graph, and UMAP.
#           Saves the complete object to disk — this is the main
#           checkpoint. 05a onwards loads from this RDS.
# Output  : results/tmp/<sample>/<sample>_post_banksy.rds
# Requires: obj, sample_name  (from 03_normalise.R)
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
repeat {
  sample_name <- trimws(readline(prompt = "\nEnter sample name (e.g. P1, P1.B): "))
  if (nchar(sample_name) == 0) {
    cat("Sample name cannot be empty.\n"); next
  }
  break
}

# Load checkpoint
# Loads the object saved at the end of 03_normalise.R
in_rds <- file.path("results", "tmp", sample_name,
                    paste0(sample_name, "_post_normalise.rds"))
if (!file.exists(in_rds)) {
  stop("File not found: ", in_rds,
       "\nPlease run 03_normalise.R first.")
}
 
obj <- readRDS(in_rds)
cat("Loaded:", in_rds, "\n")
cat("Cells:", ncol(obj), "| Genes:", nrow(obj), "\n")

# Output directory 
tmp_dir <- file.path("results", "tmp", sample_name)

if (!dir.exists(tmp_dir)) {
  dir.create(tmp_dir, recursive = TRUE)
}

# Set the default assay
DefaultAssay(obj) <- "Spatial.Polygons"

###########################################
# BANKSY - Spatial neighbourhood embedding
###########################################
# BANKSY augments each cell's gene expression with the average
# expression of its spatial neighbours (k_geom = 50 neighbours)
# lambda = 0.2 controls the neighbour contribution weight:
#   0 = pure expression, 1 = pure neighbourhood
# Result stored in a new assay: BANKSY.0.2l
obj <- RunBanksy(obj,
                 lambda     = 0.2,
                 assay_name = "BANKSY.0.2l",
                 k_geom     = 50,
                 assay      = "Spatial.Polygons")

# Set the defaualt to BANKSY assay for further analysis
DefaultAssay(obj) <- "BANKSY.0.2l"

#######
# PCA
#######
# Runs PCA on all genes in the BANKSY assay
# Uses 30 principal components
# Result stored under reduction name "pca.banksy.0.2"
obj <- RunPCA(obj,
              assay          = "BANKSY.0.2l",
              features       = rownames(obj),
              npcs           = 30,
              reduction.name = "pca.banksy.0.2")

# Visualise top genes driving the first 4 PCs
p_pca <- VizDimLoadings(obj,
                        reduction = "pca.banksy.0.2",
                        nfeatures = 5,
                        dims      = 1:4) +
  plot_annotation(title = paste("PCA loadings —", sample_name))

ggsave(file.path(tmp_dir, paste0(sample_name, "_pca_loadings.pdf")),
       p_pca, width = 12, height = 8)
cat("Saved: PCA loadings plot\n")

#############
# kNN graph
#############
# Constructs a k-nearest neighbour graph based on the PCA embedding
obj <- FindNeighbors(obj,
                     reduction = "pca.banksy.0.2",
                     dims      = 1:30)

########
# UMAP
########
# Non-linear dimensionality reduction for visualisation
obj <- RunUMAP(obj,
               reduction      = "pca.banksy.0.2",
               reduction.name = "umap.banksy.0.2",
               return.model   = TRUE,
               dims           = 1:30)

# Save the object
# At this point obj contains:
#   - Spatial.Polygons assay (raw, normalised, scaled, HVGs)
#   - BANKSY.0.2l assay (spatial neighbourhood embedding)
#   - pca.banksy.0.2 reduction
#   - umap.banksy.0.2 reduction
#   - BANKSY.0.2l_snn and _nn graphs
#   - All meta.data columns (percent.mt, nFeature_fresh, etc.)
out_rds <- file.path(tmp_dir, paste0(sample_name, "_post_banksy.rds"))
saveRDS(obj, file = out_rds)
cat("Saved:", out_rds, "\n")
cat("Proceed to 05a_cluster_clustree.R\n")