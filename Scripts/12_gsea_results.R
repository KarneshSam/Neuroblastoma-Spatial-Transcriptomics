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

# Filter by p.adjust
# Keep only statistically significant pathways
sig <- df_g %>% filter(p.adjust < 0.05)

cat("\nTotal significant pathways (p.adjust < 0.05):", nrow(sig), "\n")
cat("Per subclone:\n")
print(table(sig$subclone))
cat("Per NE type:\n")
print(table(sig$NE_type))

# Top 10 per subclone
# Ranked by p.adjust then absolute NES
sig_top10 <- sig %>%
  group_by(subclone) %>%
  arrange(p.adjust, desc(abs(NES))) %>%
  slice_head(n = 10) %>%
  ungroup()

cat("\nTop 10 per subclone — total rows:", nrow(sig_top10), "\n")
cat("NES summary in top 10:\n")
print(summary(sig_top10$NES))

# Prompt: NES threshold
# Inspect the histogram and NES summary above to choose a threshold
# Only pathways with |NES| > threshold are included in the plots
# Higher threshold = fewer, stronger pathways
# Typical values: 1.5, 2.0, 2.5
cat("\nNES distribution in top 10 pathways:\n")
hist(sig_top10$NES,
     breaks = 30,
     main   = paste("NES distribution — top 10 per subclone |", sample_name),
     xlab   = "NES",
     col    = "steelblue",
     border = "white")

cat("\nPathways with |NES| above common thresholds:\n")
for (thresh in c(1.0, 1.5, 2.0, 2.5, 3.0)) {
  n <- sum(abs(sig_top10$NES) > thresh, na.rm = TRUE)
  cat(sprintf("  |NES| > %.1f  ->  %d pathways\n", thresh, n))
}

# Keep asking until a valid numeric threshold is entered
repeat {
  nes_input <- trimws(readline(
    prompt = "\nEnter NES absolute value threshold (e.g. 2.0): "))
  nes_thresh <- suppressWarnings(as.numeric(nes_input))
  if (is.na(nes_thresh) || nes_thresh < 0) {
    cat("Please enter a positive number (e.g. 2.0).\n"); next
  }
  n_kept <- sum(abs(sig_top10$NES) > nes_thresh, na.rm = TRUE)
  cat(sprintf("  |NES| > %.1f keeps %d pathways\n", nes_thresh, n_kept))
  if (n_kept == 0) {
    cat("  No pathways kept — threshold too high. Please try a lower value.\n"); next
  }
  # Confirm if very few pathways remain
  if (n_kept < 5) {
    cat("  Warning: only", n_kept, "pathway(s) kept.\n")
    confirm <- trimws(readline(
      prompt = "  Enter 'yes' to proceed, anything else to re-enter: "))
    if (tolower(confirm) != "yes") next
  }
  break
}

cat("Using NES threshold: |NES| >", nes_thresh, "\n")

# Apply NES filter
sig_1_g <- sig_top10 %>% filter(abs(NES) > nes_thresh)

cat("Pathways after NES filter:", nrow(sig_1_g), "\n")

