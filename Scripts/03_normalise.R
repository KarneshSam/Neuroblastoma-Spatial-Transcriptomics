# 03_normalise.R
# Requires: obj, sample_name  (from 02_filter.R)
# ─────────────────────────────────────────────────────────────────

DefaultAssay(obj) <- "Spatial.Polygons"

# Normalise
obj <- NormalizeData(obj,
                     normalization.method = "LogNormalize",
                     scale.factor         = 10000)

# Highly variable genes
# Compute top 2000 variable genes
obj <- FindVariableFeatures(obj,
                            selection.method = "vst",
                            nfeatures        = 2000)

# HVG scatter plot
hvf          <- obj[["Spatial.Polygons"]]@meta.data
hvf$gene     <- rownames(obj[["Spatial.Polygons"]])
hvf$variable <- hvf$gene %in% VariableFeatures(obj)
top10        <- head(VariableFeatures(obj), 10)

ggplot(hvf, aes(x = vf_vst_counts_mean,
                y = vf_vst_counts_variance.expected,
                color = variable)) +
  geom_point(size = 1.2, alpha = 0.5) +
  scale_color_manual(
    values = c("FALSE" = "blue", "TRUE" = "red"),
    labels = c("FALSE" = "Non-variable", "TRUE" = "Variable")) +
  scale_x_log10() +
  geom_point(data  = hvf[hvf$gene %in% top10, ],
             color = "black", size = 1.2) +
  geom_text_repel(data  = hvf[hvf$gene %in% top10, ],
                  aes(label = gene), color = "black", size = 2.5) +
  labs(title = paste("Highly variable genes (2,000) —", sample_name),
       x     = "Mean expression (log scale)",
       y     = "Standardised variance",
       color = "") +
  theme_classic()

# Scale — regress out mitochondrial %
obj <- ScaleData(obj, vars.to.regress = "percent.mt")

cat("Normalisation complete. Proceed to 04_banksy_pca_umap.R\n")