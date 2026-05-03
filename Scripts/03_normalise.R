# 03_normalise.R
# Requires: obj, sample_name  (from 02_filter.R)
# ─────────────────────────────────────────────────────────────────

DefaultAssay(obj) <- "Spatial.Polygons"

# Normalise
obj <- NormalizeData(obj,
                     normalization.method = "LogNormalize",
                     scale.factor         = 10000)
