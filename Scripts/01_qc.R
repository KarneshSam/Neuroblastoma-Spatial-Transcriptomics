# 01_qc.R
# Purpose : Load libraries, read RDS, select sample, compute QC
#           metrics, and save QC plots.
# Output  : obj, sample_name, plot_dir — used by all downstream scripts
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

# Prompt: RDS file path
repeat {
  rds_path <- trimws(readline(prompt = "\nEnter path to RDS file: "))
  if (nchar(rds_path) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!file.exists(rds_path)) {
    cat("File not found: '", rds_path, "' — please check the path.\n", sep = ""); next
  }
  break
}

###########
# Load RDS
############
cat("Loading:", rds_path, "\n")
rds_obj <- readRDS(rds_path)

# Mitochondrial percentage per cell
# High mitochondrial % indicates damaged or dying cells
# Genes starting with MT- are mitochondrial in human data
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

#########################
# Per-cell violin plots
#########################
# nFeature: number of genes detected per cell
p1 <- VlnPlot(obj, features = "nFeature_Spatial.Polygons",
              pt.size = 0, group.by = "orig.ident") +
  ggtitle("Number of features") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

# nCount: total UMI counts per cell
p2 <- VlnPlot(obj, features = "nCount_Spatial.Polygons",
              pt.size = 0, group.by = "orig.ident") +
  ggtitle("Total UMI counts") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

# percent.mt: mitochondrial gene proportion per cell
p3 <- VlnPlot(obj, features = "percent.mt",
              pt.size = 0, group.by = "orig.ident") +
  ggtitle("Mitochondrial %") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

# Combine the three plots side by side with a shared title
(p1 | p2 | p3) +
  plot_annotation(
    title = paste("QC metrics —", sample_name),
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5,
                                            face = "bold")))

##################################
# Per-gene UMI count distribution
###################################
# Shows how counts are distributed across genes
gene_counts <- rowSums(obj[["Spatial.Polygons"]]$counts)

# Violin Plot
ggplot(data.frame(gene_counts), aes(x = "", y = gene_counts)) +
  geom_violin() +
  ylim(0, 10000) +
  labs(title = paste("Counts per gene —", sample_name),
       x = NULL, y = "Counts") +
  theme_classic()

cat("QC plots done. Inspect, then proceed to 02_filter.R\n")