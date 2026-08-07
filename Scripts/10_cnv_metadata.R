library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)

####################
# Load Sample name
####################
# All file paths are constructed from sample_name
# Keep asking until a non-empty string is entered
repeat {
  sample_name <- trimws(readline(prompt = "\nEnter sample name (e.g. P1, P2): "))
  if (nchar(sample_name) == 0) {
    cat("Sample name cannot be empty.\n"); next
  }
  break
}

# Prompt: post_cnv RDS file
# Seurat object saved at the end of 09_cnv_downstream_analysis.R
# Contains cnv_score and infercnv_subclone in meta.data
repeat {
  rds_path <- trimws(readline(
    prompt = paste0("\nEnter path to ", sample_name, "_post_cnv.rds: ")))
  if (nchar(rds_path) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!file.exists(rds_path)) {
    cat("File not found: '", rds_path, "'\n", sep = ""); next
  }
  break
}

obj <- readRDS(rds_path)
cat("Loaded:", rds_path,
    "| Cells:", ncol(obj), "| Genes:", nrow(obj), "\n")

# Prompt: chromosome-level csv file
# chr_wide CSV saved by 09_cnv_downstream_analysis.R
# Rows = subclones, columns = chromosomes, values = dominant CNV state
repeat {
  chr_wide_path <- trimws(readline(
    prompt = paste0("\nEnter path to ", sample_name, "_cnv_chr_wide.csv: ")))
  if (nchar(chr_wide_path) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!file.exists(chr_wide_path)) {
    cat("File not found: '", chr_wide_path, "'\n", sep = ""); next
  }
  break
}

chr_wide <- read.csv(chr_wide_path, check.names = FALSE)
cat("Loaded chr_wide:", nrow(chr_wide), "subclones x",
    ncol(chr_wide) - 1, "chromosomes\n")

# Prompt: arm-level csv file 
# arm_wide CSV saved by 09_cnv_downstream_analysis.R
# Rows = subclones, columns = chromosome arms, values = dominant CNV state
repeat {
  arm_wide_path <- trimws(readline(
    prompt = paste0("\nEnter path to ", sample_name, "_cnv_arm_wide.csv: ")))
  if (nchar(arm_wide_path) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!file.exists(arm_wide_path)) {
    cat("File not found: '", arm_wide_path, "'\n", sep = ""); next
  }
  break
}

arm_wide <- read.csv(arm_wide_path, check.names = FALSE)
cat("Loaded arm_wide:", nrow(arm_wide), "subclones x",
    ncol(arm_wide) - 1, "arms\n")

# Prompt: InferCNV output directory
# Needed to reload cell_groupings (cell barcode -> subclone mapping)
repeat {
  infercnv_dir <- trimws(readline(
    prompt = "\nEnter InferCNV output directory path: "))
  if (nchar(infercnv_dir) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!dir.exists(infercnv_dir)) {
    cat("Directory not found: '", infercnv_dir, "'\n", sep = ""); next
  }
  break
}

# Prompt: plot output directory
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

# Prompt: tumour cell type prefix
# Used to filter meta.data to tumour cells for proportion plots
# e.g. "NE" to keep NE1, NE2 ... regardless of how many exist
repeat {
  tumour_prefix <- trimws(readline(
    prompt = "\nEnter tumour cell type prefix (e.g. NE): "))
  if (nchar(tumour_prefix) == 0) {
    cat("Prefix cannot be empty.\n"); next
  }
  break
}

# Detect all tumour subtypes present in the object
# Uses cell_type_hr column — works with any number of subtypes
tumour_types <- sort(unique(obj$cell_type_hr[
  grepl(paste0("^", tumour_prefix), obj$cell_type_hr)]))
cat("Tumour subtypes detected:", paste(tumour_types, collapse = ", "), "\n")

if (length(tumour_types) == 0) {
  stop("No cell types found with prefix '", tumour_prefix, "'")
}

# Load cell groupings
# Maps each cell barcode to its InferCNV subclone assignment
hmm_prefix     <- "17_HMM_predHMMi6.leiden.hmm_mode-subclusters"
groupings_file <- file.path(infercnv_dir,
  paste0(hmm_prefix, ".cell_groupings"))

if (!file.exists(groupings_file)) {
  stop("Cell groupings file not found: ", groupings_file)
}

cell_groupings <- read.table(groupings_file,
                             header    = TRUE,
                             sep       = "\t",
                             col.names = c("subcluster", "cell_barcode"))

cat("Cell groupings loaded:", nrow(cell_groupings), "cells\n")

###############################################
# Add chr-level cnv metadata to Seurat object
###############################################
# Convert chr_wide to dataframe
chr_wide_df <- as.data.frame(chr_wide)

# Get chromosome columns — everything except the subclone name column
chr_cols <- colnames(chr_wide_df)[colnames(chr_wide_df) != "cell_group_name"]

# Build cell-level metadata from subclone assignments 
# Join chr_wide onto cell_groupings by subclone name
# Each cell inherits the CNV state of its subclone 
# NE cells with no subclone assignment get state "3" (neutral)
cell_cnv_meta <- cell_groupings %>%
  left_join(chr_wide_df, by = c("subcluster" = "cell_group_name")) %>%
  mutate(across(
    all_of(chr_cols),
    ~ ifelse(grepl(paste0("^", tumour_prefix), subcluster) &
               is.na(.x), "3", .x)
  ))

# Create cell barcode indexed dataframe 
chr_meta     <- cell_cnv_meta[, c("cell_barcode", chr_cols)]
rownames(chr_meta)    <- chr_meta$cell_barcode
chr_meta$cell_barcode <- NULL

cat("\nChr metadata dimensions:", nrow(chr_meta), "cells x",
    ncol(chr_meta), "chromosomes\n")

# Create full matrix aligned to all Seurat cells
# Initialise with NA for every cell — cells not in InferCNV stay NA
full_meta <- matrix(
  NA,
  nrow     = ncol(obj),
  ncol     = length(chr_cols),
  dimnames = list(colnames(obj), chr_cols)
)

# Fill in only matched cells
matched_cells               <- intersect(rownames(chr_meta), colnames(obj))
full_meta[matched_cells, ]  <- as.matrix(chr_meta[matched_cells, ])
full_meta_df                <- as.data.frame(full_meta)

cat("Matched cells for chr metadata:", length(matched_cells), "\n")

# Add chromosome CNV columns to Seurat object
obj <- AddMetaData(obj, metadata = full_meta_df)
cat("Chromosome CNV state columns added to metadata\n")

# Validate chr columns added
# Keep only chr columns that are actually in meta.data
chr_cols <- paste0("chr", 1:22)
chr_cols_valid <- chr_cols[chr_cols %in% colnames(obj@meta.data)]
cat("Valid chromosome columns:", length(chr_cols_valid), "\n")

# Build proportion table for all chromosomes
# For each chromosome: proportion of cells in each CNV state 
# Filtered to tumour cells only
state_prop_all <- lapply(chr_cols_valid, function(chr) {
  obj@meta.data %>%
    filter(cell_type %in% tumour_types) %>%
    filter(!is.na(.data[[chr]])) %>%
    group_by(cell_type, state = .data[[chr]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(cell_type) %>%
    mutate(prop       = round(n / sum(n), 3),
           chromosome = chr) %>%
    ungroup()
}) %>% bind_rows()

# Print per chromosome summary
for (chr in chr_cols_valid) {
  cat("\n===", chr, "===\n")
  state_prop_all %>%
    filter(chromosome == chr) %>%
    arrange(cell_type, state) %>%
    print(n = Inf)
}

# Plot: stacked bar — CNV state proportion per chromosome
# Each facet = one chromosome
# x axis = tumour subtype, y axis = proportion, fill = CNV state
p_chr_bar <- state_prop_all %>%
  mutate(
    state      = factor(state, levels = c("1","2","3","4","5","6")),
    chromosome = factor(chromosome,
                        levels = paste0("chr", c(1:22, "X"))),
    cell_type  = factor(cell_type, levels = tumour_types)
  ) %>%
  ggplot(aes(x = cell_type, y = prop, fill = state)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ chromosome, ncol = 6) +
  scale_fill_manual(
    values = c("1" = "#2166AC", "2" = "#92C5DE", "3" = "grey90",
               "4" = "#F4A582", "5" = "#D6604D", "6" = "#B2182B"),
    name   = "CNV state",
    labels = c("1" = "Complete loss", "2" = "Loss", "3" = "Neutral",
               "4" = "Gain",          "5" = "Amplification",
               "6" = "High amp")
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    strip.text      = element_text(face = "bold", size = 8),
    legend.position = "bottom"
  ) +
  labs(x     = NULL,
       y     = "Proportion of cells",
       title = paste("CNV state distribution per tumour subtype",
                     "per chromosome —", sample_name)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

chr_bar_path <- file.path(plot_dir,
  paste0(sample_name, "_cnv_chr_state_barplot.pdf"))
ggsave(chr_bar_path, p_chr_bar, width = 16, height = 10)
cat("Saved:", chr_bar_path, "\n")

###############################################
# Add arm-level cnv metadata to Seurat object
###############################################
# Convert arm_wide to dataframe
# arm_wide was saved as a CSV with cell_group_name as a column
# Convert back to matrix format for joining
arm_wide_df <- as.data.frame(arm_wide)

# Get arm columns
arm_cols <- colnames(arm_wide_df)

# Build cell-level arm metadata
# Same approach as chromosome-level — join by subclone name
# NE cells with no assignment get state "3" (neutral)
cell_cnv_arm_meta <- cell_groupings %>%
  left_join(
    arm_wide_df %>% rownames_to_column("cell_group_name"),
    by = c("subcluster" = "cell_group_name")) %>%
  mutate(across(
    all_of(arm_cols),
    ~ ifelse(grepl(paste0("^", tumour_prefix), subcluster) &
               is.na(.x), "3", .x)
  ))

# Create cell barcode indexed dataframe
arm_meta              <- cell_cnv_arm_meta[, c("cell_barcode", arm_cols)]
rownames(arm_meta)    <- arm_meta$cell_barcode
arm_meta$cell_barcode <- NULL

cat("\nArm metadata dimensions:", nrow(arm_meta), "cells x",
    ncol(arm_meta), "arms\n")

# Create full matrix aligned to all Seurat cells
full_arm_meta <- matrix(
  NA,
  nrow     = ncol(obj),
  ncol     = length(arm_cols),
  dimnames = list(colnames(obj), arm_cols)
)

matched_cells                  <- intersect(rownames(arm_meta), colnames(obj))
full_arm_meta[matched_cells, ] <- as.matrix(arm_meta[matched_cells, ])
full_arm_meta_df               <- as.data.frame(full_arm_meta)

cat("Matched cells for arm metadata:", length(matched_cells), "\n")

# Add arm CNV columns to Seurat object
obj <- AddMetaData(obj, metadata = full_arm_meta_df)
cat("Arm CNV state columns added to metadata\n")

# Validate arm columns added
arm_cols <- paste0("chr", rep(1:22, each = 2), c("p", "q"))
arm_cols_valid <- arm_cols[arm_cols %in% colnames(obj@meta.data)]
cat("Valid arm columns:", length(arm_cols_valid), "\n")

# Build proportion table for all arms 
arm_prop_all <- lapply(arm_cols_valid, function(arm) {
  obj@meta.data %>%
    filter(cell_type %in% tumour_types) %>%
    filter(!is.na(.data[[arm]])) %>%
    group_by(cell_type, state = .data[[arm]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(cell_type) %>%
    mutate(prop = round(n / sum(n), 3),
           arm  = arm) %>%
    ungroup()
}) %>% bind_rows()

# Print per arm summary
for (arm in arm_cols_valid) {
  cat("\n===", arm, "===\n")
  arm_prop_all %>%
    filter(arm == !!arm) %>%
    arrange(cell_type, state) %>%
    print(n = Inf)
}

# Plot: stacked bar — CNV state proportion per arm
p_arm_bar <- arm_prop_all %>%
  mutate(
    state     = factor(state, levels = c("1","2","3","4","5","6")),
    arm       = factor(arm, levels = arm_cols_valid),
    cell_type = factor(cell_type, levels = tumour_types)
  ) %>%
  ggplot(aes(x = cell_type, y = prop, fill = state)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ arm, ncol = 6) +
  scale_fill_manual(
    values = c("1" = "#2166AC", "2" = "#92C5DE", "3" = "grey90",
               "4" = "#F4A582", "5" = "#D6604D", "6" = "#B2182B"),
    name   = "CNV state",
    labels = c("1" = "Complete loss", "2" = "Loss", "3" = "Neutral",
               "4" = "Gain",          "5" = "Amplification",
               "6" = "High amp")
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y     = element_text(size = 7),
    strip.text      = element_text(face = "bold", size = 7),
    legend.position = "bottom"
  ) +
  labs(x     = NULL,
       y     = "Proportion of cells",
       title = paste("CNV state distribution per tumour subtype",
                     "per chromosome arm —", sample_name)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

arm_bar_path <- file.path(plot_dir,
  paste0(sample_name, "_cnv_arm_state_barplot.pdf"))
ggsave(arm_bar_path, p_arm_bar, width = 18, height = 12)
cat("Saved:", arm_bar_path, "\n")



