# 03_normalise.R
# Purpose : Normalise counts, identify highly variable genes (HVGs),
#           and scale the data while regressing out mitochondrial %.
# Output  : results/tmp/<sample>/<sample>_post_normalise.rds
#           results/tmp/<sample>/<sample>_hvg_plot.pdf
# Requires: obj, sample_name  (from 02_filter.R)
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
# Loads the object saved at the end of 02_filter.R
in_rds <- file.path("results", "tmp", sample_name,
                    paste0(sample_name, "_post_filter.rds"))

if (!file.exists(in_rds)) {
  stop("File not found: ", in_rds,
       "\nPlease run 02_filter.R first.")
}

obj <- readRDS(in_rds)
cat("Loaded:", in_rds, "\n")
cat("Cells:", ncol(obj), "| Genes:", nrow(obj), "\n")

# Output directory
tmp_dir <- file.path("results", "tmp", sample_name)

if (!dir.exists(tmp_dir)) {
  dir.create(tmp_dir, recursive = TRUE)
}

# Set default assay
DefaultAssay(obj) <- "Spatial.Polygons"

# Normalise
# LogNormalize: divides each cell's counts by total counts,
# multiplies by scale factor, then log1p transforms
obj <- NormalizeData(obj,
                     normalization.method = "LogNormalize",
                     scale.factor         = 10000)

# Highly variable genes
# VST method models mean-variance relationship and selects
# Top 2000 variable genes
obj <- FindVariableFeatures(obj,
                            selection.method = "vst",
                            nfeatures        = 2000)

##################
# Visualise HVGs
##################
# Visualise mean expression vs standardised variance

# Extract HVG metadata from the assay
hvf          <- obj[["Spatial.Polygons"]]@meta.data
# Add gene names and variable gene status for plotting
hvf$gene     <- rownames(obj[["Spatial.Polygons"]])
hvf$variable <- hvf$gene %in% VariableFeatures(obj)
# Get top 10 HVGs for labelling
top10        <- head(VariableFeatures(obj), 10)

# Red = selected variable genes, Blue = non-variable
# Top 10 HVGs are labelled by name
p_hvg <- ggplot(hvf, aes(x = vf_vst_counts_mean,
                          y = vf_vst_counts_variance.expected,
                          color = variable)) +
  geom_point(size = 1.2, alpha = 0.5) +
  scale_color_manual(
    values = c("FALSE" = "blue", "TRUE" = "red"),
    labels = c("FALSE" = "Non-variable", "TRUE" = "Variable")) +
  scale_x_log10() +
  # Highlight top 10 with black points
  geom_point(data  = hvf[hvf$gene %in% top10, ],
             color = "black", size = 1.2) +
  # Label top 10 gene names, repelled to avoid overlap
  geom_text_repel(data  = hvf[hvf$gene %in% top10, ],
                  aes(label = gene), color = "black", size = 2.5) +
  labs(title = paste("Highly variable genes (2,000) —", sample_name),
       x     = "Mean expression (log scale)",
       y     = "Standardised variance",
       color = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 13, hjust = 0.5,
                                  face = "bold"))
 
 # Save HVG plot
hvg_path <- file.path(tmp_dir, paste0(sample_name, "_hvg_plot.pdf"))
ggsave(hvg_path, p_hvg, width = 8, height = 6)
cat("Saved:", hvg_path, "\n")

# Scale — regress out mitochondrial %
obj <- ScaleData(obj, vars.to.regress = "percent.mt")

# Save checkpoint
# 04_banksy_pca_umap.R loads from this file
out_rds <- file.path(tmp_dir, paste0(sample_name, "_post_normalise.rds"))
saveRDS(obj, file = out_rds)
cat("Saved:", out_rds, "\n")

cat("Normalisation complete. Proceed to 04_banksy_pca_umap.R\n")