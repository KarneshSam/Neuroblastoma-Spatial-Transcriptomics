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
 