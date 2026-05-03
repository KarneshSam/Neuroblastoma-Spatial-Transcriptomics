# 01_qc.R
# Requires: obj, sample_name  (from 00_setup.R)
# ─────────────────────────────────────────────────────────────────

# Mitochondrial percentage per cell
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Per-cell violin plots
p1 <- VlnPlot(obj, features = "nFeature_Spatial.Polygons",
              pt.size = 0, group.by = "orig.ident") +
  ggtitle("Number of features") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

p2 <- VlnPlot(obj, features = "nCount_Spatial.Polygons",
              pt.size = 0, group.by = "orig.ident") +
  ggtitle("Total UMI counts") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

p3 <- VlnPlot(obj, features = "percent.mt",
              pt.size = 0, group.by = "orig.ident") +
  ggtitle("Mitochondrial %") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

(p1 | p2 | p3) +
  plot_annotation(
    title = paste("QC metrics —", sample_name),
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5,
                                            face = "bold"))
  )

