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