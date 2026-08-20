# 13_neighbourhood.R
# Purpose : Compute spatial neighbourhood composition for each NE
#           subtype using a fixed-radius nearest-neighbour search.
#           Generates a stacked bar plot per NE subtype coloured by
#           broad cell class (Tumour, Stromal, Endothelial, Immune).
# Output  : results/neighbour/<sample>_neighbourhood_composition.pdf
# ─────────────────────────────────────────────────────────────────

library(Seurat)
library(RANN)
library(tidyverse)
library(ggh4x)

# Prompt: sample name
repeat {
  sample_name <- trimws(readline(prompt = "\nEnter sample name (e.g. P1, P1.B): "))
  if (nchar(sample_name) == 0) {
    cat("Sample name cannot be empty.\n"); next
  }
  break
}

# Detect available annotated objects
# Lists all annotated RDS files saved by 07_annotation.R
final_files <- list.files(
  path       = file.path("results", "annotated_obj"),
  pattern    = paste0("^", sample_name, "_annotated_.*\\.rds$"),
  full.names = TRUE
)

if (length(final_files) == 0) {
  stop("No annotated RDS files found in results/annotated_obj/", sample_name,
       "/\nPlease run 07_annotation.R first.")
}
 
cat("\nAvailable annotated objects for", sample_name, ":\n")
for (i in seq_along(final_files)) {
  cat(sprintf("  [%d] %s\n", i, basename(final_files[i])))
}

# Prompt: choose file
repeat {
  chosen_file <- trimws(readline(
    prompt = "\nEnter the filename to load (e.g. P1_annotated_banksy_lam0.2_res0.5.rds): "))
  if (chosen_file %in% basename(final_files)) break
  cat("Invalid — please enter one of the listed filenames.\n")
}
 
in_rds <- file.path("results", "annotated_obj", chosen_file)
obj    <- readRDS(in_rds)
cat("Loaded:", in_rds, "\n")
cat("Cells:", ncol(obj), "| Genes:", nrow(obj), "\n")

# Set default assay
DefaultAssay(obj) <- "Spatial.Polygons"

# Prompt: search radius
# Radius in µm for the fixed-radius nearest-neighbour search
# Default is 200 µm based on the study methodology
repeat {
  radius_input <- trimws(readline(
    prompt = "\nEnter spatial search radius in µm [default: 200]: "))
  if (nchar(radius_input) == 0) {
    radius <- 200
    cat("Using default: 200 µm\n")
    break
  }
  radius <- suppressWarnings(as.numeric(radius_input))
  if (is.na(radius) || radius <= 0) {
    cat("Invalid — please enter a positive number.\n"); next
  }
  break
}
 
# Prompt: tumour subtype prefix
# Used to detect which cell types to run neighbourhood analysis for
# e.g. "NE" will find NE1, NE2, NE3 etc. from cell_type_hr
repeat {
  tumour_prefix <- trimws(readline(
    prompt = "\nEnter tumour subtype prefix (e.g. NE): "))
  if (nchar(tumour_prefix) == 0) {
    cat("Prefix cannot be empty.\n"); next
  }
  break
}
 
# Detect all tumour subtypes present in cell_type_hr
ne_types <- sort(unique(obj$cell_type_hr[
  grepl(paste0("^", tumour_prefix), obj$cell_type_hr)]))
cat("Tumour subtypes detected:", paste(ne_types, collapse = ", "), "\n")
 
if (length(ne_types) == 0) {
  stop("No cell types found with prefix '", tumour_prefix,
       "' in cell_type_hr column.")
}

# Output directory
out_dir <- file.path("results", "neighbour")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("Created directory:", out_dir, "\n")
}

# Step 1: Extract coordinates and cell type annotations
# GetTissueCoordinates returns spatial x/y coordinates per cell
# Match cell_type_hr from metadata to the coordinate table
coords <- GetTissueCoordinates(obj)
coords$cell_type <- obj$cell_type_hr[
  match(coords$cell, colnames(obj))
]
 
# Remove cells with no cell type annotation
coords_all     <- coords[!is.na(coords$cell_type), ]
coords_mat_all <- as.matrix(coords_all[, c("x", "y")])
 
cat("\nTotal cells with coordinates:", nrow(coords_all), "\n")
 
# ── Step 2: Neighbourhood function ───────────────────────────────
# For each cell of the target type, find all cells within radius µm
# Exclude the query cell itself from its own neighbour list
# Pool all neighbours across all query cells and compute type proportions
get_neighbourhood <- function(target_type, radius = 200) {
 
  target_idx    <- which(coords_all$cell_type == target_type)
  target_coords <- coords_mat_all[target_idx, ]
 
  # Fixed-radius nearest-neighbour search
  # k = 100 sets an upper limit on neighbours returned per cell
  nn <- RANN::nn2(
    data       = coords_mat_all,
    query      = target_coords,
    k          = 100,
    searchtype = "radius",
    radius     = radius
  )
 
  # For each query cell, collect neighbour cell types
  # Remove zero-index entries (unfilled slots) and self
  neighbour_types <- lapply(seq_len(nrow(nn$nn.idx)), function(i) {
    idx <- nn$nn.idx[i, ]
    idx <- idx[idx > 0]           # remove empty slots
    idx <- idx[idx != target_idx[i]]  # exclude self
    coords_all$cell_type[idx]
  })
 
  # Pool all neighbour types and compute relative frequencies
  all_neighbours <- unlist(neighbour_types)
  type_counts    <- table(all_neighbours)
  type_freq      <- prop.table(type_counts)
 
  data.frame(
    neighbour_type = names(type_freq),
    frequency      = as.numeric(type_freq),
    target         = target_type
  )
}