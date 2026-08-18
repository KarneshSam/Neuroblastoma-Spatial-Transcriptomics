suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(tibble)
})

########################
# Prompt: sample name 
#########################
repeat {
  sample_name <- trimws(readline(prompt = "\nEnter sample name (e.g. P1, P2): "))
  if (nchar(sample_name) == 0) {
    cat("Sample name cannot be empty.\n"); next
  }
  break
}

# Prompt: NE cells RDS path
# Loads the NE-only subset saved at the end of 09_cnv_metadata.R
# File is named: <sample_name>_ne_cells.rds
repeat {
  rds_path <- trimws(readline(
    prompt = paste0("\nEnter path to ", sample_name, "_ne_cells.rds: ")))
  if (nchar(rds_path) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!file.exists(rds_path)) {
    cat("File not found: '", rds_path, "'\n", sep = ""); next
  }
  break
}

############################
# Prompt: output directory
############################
# All GSEA results and CSVs will be saved here
repeat {
  output_dir <- trimws(readline(
    prompt = "\nEnter output directory for GSEA results: "))
  if (nchar(output_dir) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  break
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}


# Prompt: minimum cells per subclone
# Subclones with fewer cells than this threshold are excluded
# Avoids running GSEA on unreliable average expression estimates
repeat {
  min_cells_input <- trimws(readline(
    prompt = "\nMinimum cells per subclone to include (default 10): "))

  # Accept empty input — use default of 10
  if (nchar(min_cells_input) == 0) {
    min_cells <- 10
    cat("Using default: 10\n")
    break
  }

  min_cells <- suppressWarnings(as.integer(min_cells_input))
  if (is.na(min_cells) || min_cells < 1) {
    cat("Please enter a positive integer.\n"); next
  }
  break
}

# Load NE cells object 
cat("\nLoading:", rds_path, "\n")
ne_cells <- readRDS(rds_path)
cat("Loaded:", sample_name,
    "| Cells:", ncol(ne_cells),
    "| Genes:", nrow(ne_cells), "\n")

# Filter valid subclones (cell counts >= min_cells)
# Remove subclones with too few cells for reliable average expression
cat("Setting identities to infercnv_subclone...\n")
Idents(ne_cells) <- "infercnv_subclone"

subclone_sizes  <- table(ne_cells@meta.data$infercnv_subclone)
valid_subclones <- names(subclone_sizes[subclone_sizes >= min_cells])

cat("Valid subclones (>=", min_cells, "cells):", length(valid_subclones), "\n")
cat("Removed subclones (<", min_cells, "cells):",
    sum(subclone_sizes < min_cells), "\n")

# Subset NE cells to only include valid subclones
ne_cells_filtered <- subset(ne_cells,
                             infercnv_subclone %in% valid_subclones)

###################################
# Average expression per subclone
####################################
# Uses normalised 'data' layer — log-normalised counts
# Reflects each subclone's own absolute expression profile
# NOT one-vs-rest DE — avoids penalising genes expressed everywhere
cat("\nComputing average expression per subclone...\n")

# Run NormalizeData first if data layer is missing
if (!"data" %in% Layers(ne_cells_filtered, assay = "Spatial.Polygons")) {
  cat("No 'data' layer found — running NormalizeData() first...\n")
  ne_cells_filtered <- NormalizeData(ne_cells_filtered,
                                     assay = "Spatial.Polygons")
}

avg_expr_list <- AverageExpression(
  ne_cells_filtered,
  assays   = "Spatial.Polygons",
  layer    = "data",
  group.by = "infercnv_subclone"
)

