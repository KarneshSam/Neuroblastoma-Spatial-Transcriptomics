# 04_banksy_pca_umap.R
# Requires: obj, sample_name  (from 03_normalise.R)
# ─────────────────────────────────────────────────────────────────

DefaultAssay(obj) <- "Spatial.Polygons"

# BANKSY spatial neighbourhood embedding
obj <- RunBanksy(obj,
                 lambda     = 0.2,
                 assay_name = "BANKSY.0.2l",
                 k_geom     = 50,
                 assay      = "Spatial.Polygons")

# Set the defaualt to BANKSY assay for further analysis
DefaultAssay(obj) <- "BANKSY.0.2l"

# PCA
obj <- RunPCA(obj,
              assay          = "BANKSY.0.2l",
              features       = rownames(obj),
              npcs           = 30,
              reduction.name = "pca.banksy.0.2")

VizDimLoadings(obj,
               reduction = "pca.banksy.0.2",
               nfeatures = 5,
               dims      = 1:4) +
  plot_annotation(title = paste("PCA loadings —", sample_name))


