# 11_gsea_results.R
# Purpose : Load GSEA results from 10_gsea_subclone.R, filter by
#           p.adjust and user-defined NES threshold, generate
#           heatmap and bubble plot, GSEA-CNV plotand save all outputs.
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

# Parse leading edge
# Extracts tags, list, signal percentages from leading_edge string
# tags   : % of pathway genes in the leading edge
# list   : % of ranked list at the point of max enrichment
# signal : enrichment signal strength
sig_1_g <- sig_1_g %>%
  mutate(
    tags   = as.numeric(str_extract(leading_edge, "(?<=tags=)\\d+")),
    list   = as.numeric(str_extract(leading_edge, "(?<=list=)\\d+")),
    signal = as.numeric(str_extract(leading_edge, "(?<=signal=)\\d+"))
  ) %>%
  # Leading edge gene count = (tags% / 100) * pathway gene set size
  mutate(leading_edge_count = round((tags / 100) * setSize))

# Add direction column
# Derived from NES sign — positive = Activated, negative = Suppressed
sig_1_g <- sig_1_g %>%
  mutate(
    direction = case_when(
      NES > 0 ~ "Activated",
      NES < 0 ~ "Suppressed",
      TRUE    ~ NA_character_
    )
  ) %>%
  filter(!is.na(direction))

###########################################
# Heatmap - Pathway enrichment per subclone
############################################

# Build direction matrix
# Rows = pathways, columns = subclones
# Values: +2 = activated, -2 = suppressed, 0 = absent
heatmap_base_gsea <- sig_1_g %>%
  dplyr::select(Description, subclone, NES, direction) %>%
  distinct() %>%
  mutate(value = ifelse(NES > 0, nes_thresh, -nes_thresh))

mat_direction <- heatmap_base_gsea %>%
  dplyr::select(Description, subclone, value) %>%
  pivot_wider(names_from  = subclone,
              values_from = value,
              values_fill = 0) %>%
  column_to_rownames("Description") %>%
  as.matrix()

# Column annotation
# Annotates each subclone column with its NE type
col_ann <- sig_1_g %>%
  dplyr::select(subclone, NE_type) %>%
  distinct() %>%
  filter(subclone %in% colnames(mat_direction)) %>%
  arrange(NE_type) %>%
  column_to_rownames("subclone")

# Order columns by NE type
col_order     <- rownames(col_ann)
mat_direction <- mat_direction[, col_order]

# ── NE type colours — dynamic ─────────────────────────────────────
n_ne <- n_distinct(col_ann_gsea$NE_type)
ne_colors_gsea <- setNames(
  brewer.pal(max(n_ne, 3), "Set1")[seq_len(n_ne)],
  unique(sort(col_ann_gsea$NE_type))
)
ann_colors_gsea <- list(NE_type = ne_colors_gsea)

# Plot: heatmap
# breaks and legend_breaks both derived from nes_thresh
heatmap_path <- file.path(plot_dir,
  paste0(sample_name, "_gsea_pathway_heatmap.pdf"))

pdf(heatmap_path, width = 16, height = 12)
pheatmap(mat_direction,
         annotation_col     = col_ann_gsea,
         annotation_colors  = ann_colors_gsea,
         color              = c("#2166AC", "#D9D9D9", "#B2182B"),
         breaks             = c(-(nes_thresh + 0.5),
                                      -0.5,
                                       0.5,
                                       nes_thresh + 0.5),
         cluster_rows             = TRUE,
         clustering_distance_cols = "binary",
         clustering_method        = "ward.D2",
         show_colnames            = TRUE,
         show_rownames            = TRUE,
         fontsize_row             = 6,
         fontsize_col             = 6,
         angle_col                = 90,
         main                     = paste("GSEA Reactome Pathways x Subclones |",
                                          sample_name,
                                          "\n(Activated | Absent | Suppressed)",
                                          "| |NES| >", nes_thresh),
         gaps_col                 = cumsum(table(col_ann_gsea$NE_type)),
         legend_breaks            = c(-nes_thresh, 0, nes_thresh),
         legend_labels            = c("Suppressed", "Absent", "Activated"),
         border_color             = "grey90")
dev.off()
cat("Saved:", heatmap_path, "\n")

#################################################
# Bubble plot - Pathway enrichment per subclone
#################################################
# Prepare bubble data 
# Truncate long pathway names for readability
bubble_overall <- sig_1_g %>%
  mutate(
    NE_type           = str_extract(subclone, "^[^.]+"),
    Description_short = case_when(
      str_length(Description) > 40 ~ paste0(str_sub(Description, 1, 40), "..."),
      TRUE                         ~ Description
    ),
    direction = case_when(
      NES > 0 ~ "Activated",
      NES < 0 ~ "Suppressed",
      TRUE    ~ NA_character_
    )
  ) %>%
  filter(!is.na(direction))

# Order subclones by NE type then name for consistent x axis ordering
# avoids reorder() warning when column is non-numeric
bubble_overall <- bubble_overall %>%
  mutate(subclone = factor(subclone,
                           levels = unique(subclone[order(NE_type, subclone)])))

# ── Build bubble plot ─────────────────────────────────────────────
# Size  = leading edge gene count (signal strength)
# Fill  = p.adjust (red = more significant)
# Shape = direction (triangle up = activated, down = suppressed)
p_bubble <- ggplot(bubble_overall,
                   aes(x     = subclone,
                       y     = Description_short,
                       size  = leading_edge_count,
                       fill  = p.adjust,
                       shape = direction)) +
  geom_point(alpha = 0.8, stroke = 0.3, color = "grey30") +
  scale_size_continuous(
    range = c(2, 8),
    name  = "Leading edge\ngene count",
    guide = guide_legend(override.aes = list(shape = 24, fill = "grey50"))
  ) +
  scale_fill_gradientn(
    colours = c("#B2182B", "#F4A582", "#D1E5F0"),
    values  = scales::rescale(c(0, 0.01, 0.05)),
    name    = "p.adjust",
    limits  = c(0, 0.05)
  ) +
  scale_shape_manual(
    values       = c("Activated" = 24, "Suppressed" = 25),
    na.translate = FALSE,
    name         = "Direction"
  ) +
  facet_grid(. ~ NE_type,
             scales = "free_x",
             space  = "free_x") +
  labs(
    title = paste("Reactome Pathway Enrichment per Subclone |",
                  sample_name, "| |NES| >", nes_thresh),
    x     = "Subclone",
    y     = "Reactome Pathway"
  ) +
  theme_bw() +
  theme(
    axis.text.x   = element_text(angle = 90, hjust = 1,
                                 vjust = 0.5, size = 7),
    axis.text.y   = element_text(size = 7),
    strip.text    = element_text(face = "bold", size = 10),
    plot.title    = element_text(face = "bold", hjust = 0.5),
    panel.spacing = unit(0.3, "lines"),
    legend.title  = element_text(face = "bold")
  )

# ── Save bubble plot ──────────────────────────────────────────────
# Width scales with number of subclones, height with number of pathways
# limitsize = FALSE allows very wide plots when many subclones exist
bubble_path <- file.path(plot_dir,
  paste0(sample_name, "_gsea_bubble_plot.pdf"))

ggsave(
  filename  = bubble_path,
  plot      = p_bubble,
  width     = max(10, n_distinct(bubble_overall$subclone) * 0.4),
  height    = max(8,  n_distinct(bubble_overall$Description_short) * 0.2),
  units     = "in",
  dpi       = 300,
  limitsize = FALSE
)
cat("Saved:", bubble_path, "\n")

# Summary
cat("\n========== SUMMARY ==========\n")
cat("Sample:                  ", sample_name, "\n")
cat("Total pathways loaded:   ", nrow(df_g), "\n")
cat("After p.adjust < 0.05:  ", nrow(sig), "\n")
cat("After top 10 per clone: ", nrow(sig_top10), "\n")
cat("After |NES| >", nes_thresh, ":  ", nrow(sig_1_g), "\n")
cat("Activated:               ", sum(sig_1_g$direction == "Activated"), "\n")
cat("Suppressed:              ", sum(sig_1_g$direction == "Suppressed"), "\n")
cat("Plots saved to:          ", plot_dir, "\n")
cat("==============================\n")

#######################################
# GSEA - CNV genes relationship analysis
########################################

# Prompt: InferCNV CNV genes file
# pred_cnv_genes.dat from the InferCNV HMM output directory
repeat {
  cnv_genes_path <- trimws(readline(
    prompt = "\nEnter path to pred_cnv_genes.dat: "))
  if (nchar(cnv_genes_path) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!file.exists(cnv_genes_path)) {
    cat("File not found: '", cnv_genes_path, "'\n", sep = ""); next
  }
  break
}

cnv_genes <- read.table(cnv_genes_path, header = TRUE, sep = "\t")
cat("CNV genes loaded:", nrow(cnv_genes), "rows\n")
head(cnv_genes)

