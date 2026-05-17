# 00_setup.R
# Purpose : Load all libraries, read the RDS.
# Output  : obj, sample_name — used by all downstream scripts
# ─────────────────────────────────────────────────────────────────

# Libraries
library(Seurat)            # Core single-cell analysis
library(ggplot2)           # Plotting
library(ggrepel)           # Non-overlapping labels on plots
library(patchwork)         # Combine multiple plots
library(dplyr)             # Data wrangling
library(leiden)            # Clustering algorithm
library(Banksy)            # Spatial analysis tools
library(pheatmap)          # Heatmaps
library(SeuratWrappers)    # Additional Seurat functionalities
library(clustree)          # Visualize cluster relationships
library(Matrix)            # Sparse matrix handling

# Load the list of Seurat object
seu <- readRDS("seu_enact_seg_unprocessed.rds")

# Assign each sample to a variable for easy access
P1 <- seu[[1]]
P1.B <- seu[[2]]
P2 <- seu[[3]]
P2.B <- seu[[4]]

########################
# Prompt: choose sample
########################
# Get from the user --> select the dataset (Seurat object)
repeat {
  cat("\nAvailable samples: P1, P1.B, P2, P2.B\n")
  sample_name <- trimws(readline(
    prompt = "Which sample to run? (P1 / P1.B / P2 / P2.B): "))
  if (sample_name %in% c("P1", "P1.B", "P2", "P2.B")) break
  cat("Invalid input '", sample_name, "' — please enter P1, P1.B, P2, or P2.B.\n", sep = "")
}

# Set it to the obj variable for downstream use
obj <- get(sample_name)   

# Print the matrix (genes x cells)
cat("\nLoaded:", sample_name,
    "| Cells:", ncol(obj),
    "| Genes:", nrow(obj), "\n")
cat("Proceed to 01_qc.R\n")
