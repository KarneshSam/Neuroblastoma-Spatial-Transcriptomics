# 11_gsea_results.R
# Purpose : Load GSEA results from 10_gsea_subclone.R, filter by
#           p.adjust and user-defined NES threshold, generate
#           heatmap and bubble plot, and save all outputs.
# Requires: sample_name  (from 00_setup.R)
#           gsea_all_subclones.csv from 10_gsea_subclone.R output
# ─────────────────────────────────────────────────────────────────

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(tibble)
library(scales)

# Prompt: sample name
repeat {
  sample_name <- trimws(readline(prompt = "\nEnter sample name (e.g. P1, P2): "))
  if (nchar(sample_name) == 0) {
    cat("Sample name cannot be empty.\n"); next
  }
  break
}

# Prompt: GSEA results CSV
# gsea_all_subclones.csv saved by 10_gsea_subclone.R
repeat {
  gsea_path <- trimws(readline(
    prompt = "\nEnter path to gsea_all_subclones.csv: "))
  if (nchar(gsea_path) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!file.exists(gsea_path)) {
    cat("File not found: '", gsea_path, "'\n", sep = ""); next
  }
  break
}

# Prompt: output directory
repeat {
  plot_dir <- trimws(readline(
    prompt = "\nEnter directory to save plots: "))
  if (nchar(plot_dir) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  break
}

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
  cat("Created plot directory:", plot_dir, "\n")
}

# Load GSEA results
df_g <- read.csv(gsea_path)
cat("\nLoaded:", nrow(df_g), "pathways across",
    n_distinct(df_g$subclone), "subclones\n")

# Quick overview
cat("\nPathways per subclone:\n")
print(table(df_g$subclone))

cat("\nNES summary:\n")
print(summary(df_g$NES))
