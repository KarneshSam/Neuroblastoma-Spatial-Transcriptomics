library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(pheatmap)
library(ggplot2)

##################################
# Load the list of Seurat object
##################################
rds_obj <- readRDS("seu_enact_seg_unprocessed.rds")

# Handle single object or list of any length
if (inherits(rds_obj, "Seurat")) {
  cat("Detected: single Seurat object\n")
  sample_list <- list(S1 = rds_obj)

} else if (is.list(rds_obj) &&
           all(sapply(rds_obj, inherits, "Seurat"))) {
  cat("Detected:", length(rds_obj), "Seurat objects in list\n")
  if (is.null(names(rds_obj))) {
    names(rds_obj) <- paste0("S", seq_along(rds_obj))
    cat("No names found — assigned:",
        paste(names(rds_obj), collapse = ", "), "\n")
  }
  sample_list <- rds_obj

} else {
  stop("RDS does not contain a Seurat object or a list of Seurat objects.")
}

# Print summary of all samples
cat("\nSamples available:\n")
for (nm in names(sample_list)) {
  cat(sprintf("  %-6s | Cells: %d | Genes: %d\n",
              nm, ncol(sample_list[[nm]]), nrow(sample_list[[nm]])))
}

# Prompt: choose sample 
repeat {
  cat("\nAvailable samples:", paste(names(sample_list), collapse = ", "), "\n")
  sample_name <- trimws(readline(prompt = "Which sample to analyse? "))
  if (sample_name %in% names(sample_list)) break
  cat("Invalid '", sample_name,
      "' — please enter one of the listed names.\n", sep = "")
}

obj <- sample_list[[sample_name]]
cat("\nLoaded:", sample_name,
    "| Cells:", ncol(obj),
    "| Genes:", nrow(obj), "\n")

#####################################
# Prompt: InferCNV output directory
#####################################
# Directory containing the InferCNV run outputs
# Should contain run.final.infercnv_obj and HMM output files
repeat {
  infercnv_dir <- trimws(readline(
    prompt = "\nEnter InferCNV output directory path: "))
  if (nchar(infercnv_dir) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  if (!dir.exists(infercnv_dir)) {
    cat("Directory not found: '", infercnv_dir,
        "' — please check the path.\n", sep = ""); next
  }
  break
}

# Prompt: plot output directory
# Where all plots from this script will be saved
repeat {
  plot_dir <- trimws(readline(
    prompt = "\nEnter directory to save plots: "))
  if (nchar(plot_dir) == 0) {
    cat("Path cannot be empty.\n"); next
  }
  break
}

# Create plot directory if it does not exist
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
  cat("Created plot directory:", plot_dir, "\n")
}

##############################################
# Prompt: Tumour & Normal cell type prefixes
##############################################
# Used to filter CNV regions to tumour subclones only
# e.g. "NE" to keep NE1, NE2, NE3 etc.
repeat {
  tumour_prefix <- trimws(readline(
    prompt = "\nEnter tumour subclone prefix to filter CNV regions (e.g. NE): "))
  if (nchar(tumour_prefix) == 0) {
    cat("Prefix cannot be empty.\n"); next
  }
  break
}

# The Infercnv produce subclone names with prefixes for normal cell types too
# These prefixes will be simplified to just the cell type name 
# e.g. "Plasma,Endo,VSMCs,TAM1,TAM2" -> kept the first part alone "."
repeat {
  cat("\nEnter normal cell type prefixes to simplify subclone names.\n")
  cat("Comma-separated (e.g. Plasma,Endo,VSMCs,TAM1,TAM2)\n")
  normal_input   <- trimws(readline(prompt = "Normal prefixes: "))
  normal_prefixes <- trimws(strsplit(normal_input, ",")[[1]])
  normal_prefixes <- normal_prefixes[nzchar(normal_prefixes)]
  if (length(normal_prefixes) == 0) {
    cat("Please enter at least one prefix.\n"); next
  }
  break
}

# Build regex pattern from the normal cell type prefixes
# e.g. c("Plasma","Endo") -> "^(Plasma|Endo)"
normal_pattern <- paste0("^(", paste(normal_prefixes, collapse = "|"), ")")
cat("Normal subclone pattern:", normal_pattern, "\n")

#####################
# InferCNV analysis
#####################
# Load final InferCNV object from the specified directory
# Contains the CNV expression matrix (genes x cells)
infercnv_rds <- file.path(infercnv_dir, "run.final.infercnv_obj")
if (!file.exists(infercnv_rds)) {
  stop("run.final.infercnv_obj not found in: ", infercnv_dir)
}

cat("\nLoading InferCNV object...\n")
infercnv_obj <- readRDS(infercnv_rds)

# CNV score per cell
# Variance of CNV score across genes per cell
# Higher variance = more CNV events = more aneuploid
# Neutral (diploid) cells have CNV values close to 1 = low variance
cnv_mat   <- infercnv_obj@expr.data

cat("CNV matrix dimensions:", nrow(cnv_mat), "genes x",
    ncol(cnv_mat), "cells\n")
# Compute CNV score for each cell based on the variance
cnv_score <- apply(cnv_mat, 2, var)

# Assign CNV scores back to Seurat object
# Only cells in the Infercnv object will have a score
# Cells not in InferCNV (e.g. filtered out) get NA
cells_in_infercnv     <- names(cnv_score)
obj$cnv_score         <- NA
obj$cnv_score[match(cells_in_infercnv, colnames(obj))] <- cnv_score

cat("CNV scores assigned to", sum(!is.na(obj$cnv_score)), "cells\n")

#############################
# HMM subcluster assignments
##############################
# Maps each cell barcode to its InferCNV subclone assignment
hmm_prefix    <- "17_HMM_predHMMi6.leiden.hmm_mode-subclusters"
groupings_file <- file.path(infercnv_dir,
  paste0(hmm_prefix, ".cell_groupings"))

if (!file.exists(groupings_file)) {
  stop("Cell groupings file not found: ", groupings_file)
}

cell_groupings <- read.table(groupings_file,
                             header    = TRUE,
                             sep       = "\t",
                             col.names = c("subcluster", "cell_barcode"))
                    
cat("\nCell groupings loaded:", nrow(cell_groupings), "cells\n")
cat("Unique subclones:", length(unique(cell_groupings$subcluster)), "\n")

# Assign subclone labels to Seurat object 
# Initialise as NA then fill from cell_groupings
obj$infercnv_subclone <- NA
obj$infercnv_subclone[
  match(cell_groupings$cell_barcode, colnames(obj))
] <- as.character(cell_groupings$subcluster)

# Simplify normal subclone names
# Normal cells get long subclone names like "Endo.s1.c2"
# Simplify to just the cell type prefix e.g. "Endo"
obj$infercnv_subclone <- ifelse(
  grepl(normal_pattern, obj$infercnv_subclone),
  str_extract(obj$infercnv_subclone, "^[^.]+"),
  obj$infercnv_subclone)

# Fill NA with the cell type label
# Cells not in InferCNV keep their annotated cell type names
obj$infercnv_subclone <- ifelse(
  is.na(obj$infercnv_subclone),
  as.character(obj$cell_type),
  obj$infercnv_subclone)

cat("\nSubclone distribution:\n")
print(table(obj$infercnv_subclone, useNA = "ifany"))

################
# HMM CNV genes
################
# Contains per-gene CNV predictions per subclone
cnv_genes_file <- file.path(infercnv_dir,
  paste0(hmm_prefix, ".pred_cnv_genes.dat"))

if (!file.exists(cnv_genes_file)) {
  stop("CNV genes file not found: ", cnv_genes_file)
}

cnv_genes <- read.table(cnv_genes_file, header = TRUE, sep = "\t")
cat("CNV genes table loaded:", nrow(cnv_genes), "rows\n")

#############################
# CNV regions - Gain or Loss
#############################
# Contains per-region CNV state predictions per subclone
# State 3 = neutral (diploid), <3 = loss, >3 = gain
cnv_regions_file <- file.path(infercnv_dir,
  paste0(hmm_prefix, ".pred_cnv_regions.dat"))

if (!file.exists(cnv_regions_file)) {
  stop("CNV regions file not found: ", cnv_regions_file)
}

cnv_regions <- read.table(cnv_regions_file, header = TRUE, sep = "\t")
cat("\nCNV regions loaded:", nrow(cnv_regions), "rows\n")

# Filter to tumour subclones only
# Keep only rows where subclone name starts with tumour prefix
cnv_regions_tumour <- cnv_regions[
  grepl(paste0("^", tumour_prefix), cnv_regions$cell_group_name), ]

n_subclones <- length(unique(cnv_regions_tumour$cell_group_name))
cat("Tumour subclones found:", n_subclones, "\n")

if (n_subclones == 0) {
  stop("No subclones found with prefix '", tumour_prefix,
       "' — check tumour prefix.")
}

########################################
# CNV state per subclone per chromosome
########################################
# Summarise dominant CNV state per subclone per chromosome
# Excludes neutral state (3) before finding dominant state
# Dominant state = most frequent non-neutral state in that chromosome
chr_summary <- cnv_regions_tumour %>%
  filter(state != 3) %>%
  group_by(cell_group_name, chr) %>%
  summarise(
    dominant_state = as.character(names(which.max(table(state)))),
    .groups = "drop"
  )

cat("\nChromosome summary dimensions:",
    nrow(chr_summary), "rows\n")

# Verify one row per subclone + chromosome combination
dups <- chr_summary %>%
  count(cell_group_name, chr) %>%
  filter(n > 1)
if (nrow(dups) > 0) {
  warning("Duplicate subclone + chromosome combinations found: ",
          nrow(dups))
} else {
  cat("Verified: no duplicate subclone + chromosome combinations\n")
}

#################################
# Chromosome matrix for heatmap
#################################
# Pivot to wide matrix
# Rows = subclones, columns = chromosomes
# Missing = neutral (state 3, filled as "3")
chr_wide <- chr_summary %>%
  pivot_wider(
    names_from  = chr,
    values_from = dominant_state,
    values_fill = "3")

# Save chr_wide table
# Wide matrix of dominant CNV state per subclone per chromosome
# Rows = subclones, columns = chromosomes, values = CNV state
chr_wide_path <- file.path(plot_dir,
  paste0(sample_name, "_cnv_chr_wide.csv"))

write.csv(chr_wide,
          file      = chr_wide_path,
          row.names = FALSE)

cat("Saved:", chr_wide_path, "\n")

# Order columns by genomic position
# chr1 -> chr22 -> chrX — standard genomic order
chr_order <- paste0("chr", c(1:22, "X"))
chr_order <- chr_order[chr_order %in% colnames(chr_wide)]

# Convert to numeric matrix
# pheatmap requires a numeric matrix
chr_mat <- chr_wide %>%
  column_to_rownames("cell_group_name") %>%
  dplyr::select(all_of(chr_order)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

cat("CNV matrix:", nrow(chr_mat), "subclones x",
    ncol(chr_mat), "chromosomes\n")

# Row annotation: NE subtype
# Extracts the cell type prefix from subclone name
# e.g. NE1.NE1_s1       NE1
#      NE1.NE1_s10      NE1
#      ....             ..
row_ann_chr <- data.frame(
  NE_type   = gsub("\\..*", "", rownames(chr_mat)),
  row.names = rownames(chr_mat))

# NE subtype colours — generated dynamically
# Extract unique NE subtypes present in this sample's subclones
# Handles any number of NE types — not hardcoded to 8
ne_types <- sort(unique(row_ann_chr$NE_type))
cat("NE subtypes found:", paste(ne_types, collapse = ", "), "\n")

# Generate a colour palette scaled to however many NE types exist
# RColorBrewer "Set1" up to 9, beyond that use colorRampPalette
if (length(ne_types) <= 9) {
  ne_palette <- RColorBrewer::brewer.pal(
    n    = max(length(ne_types), 3),  # brewer.pal needs at least 3
    name = "Set1"
  )[seq_along(ne_types)]
} else {
  # More than 9 types — interpolate a larger palette
  ne_palette <- colorRampPalette(
    RColorBrewer::brewer.pal(9, "Set1")
  )(length(ne_types))
}

# Name the palette by NE type
names(ne_palette) <- ne_types

ne_colors <- list(NE_type = ne_palette)

cat("Colours assigned:\n")
for (nm in names(ne_palette)) {
  cat(sprintf("  %-6s -> %s\n", nm, ne_palette[nm]))
}

# CNV state colours: blue (loss) -> white (neutral) -> red (gain)
# States 1-6: 1=deep loss, 2=loss, 3=neutral, 4=gain, 5=amp, 6=high amp
cnv_colors <- colorRampPalette(c(
  "#2166AC",   # state 1 — deep loss
  "#92C5DE",   # state 2 — loss
  "white",     # state 3 — neutral
  "#F4A582",   # state 4 — gain
  "#D6604D",   # state 5 — amplification
  "#B2182B"    # state 6 — high amplification
))(6)

###############################
# Plot: CNV landscape heatmap
###############################
# Rows = tumour subclones, columns = chromosomes
# Colour = dominant CNV state per chromosome per subclone
cnv_heatmap_path <- file.path(plot_dir,
  paste0(sample_name, "_cnv_landscape_heatmap.pdf"))

pdf(cnv_heatmap_path, width = 14, height = 8)

pheatmap(chr_mat,
         color             = cnv_colors,
         breaks            = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
         annotation_row    = row_ann_chr,
         annotation_colors = ne_colors,
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         fontsize_col      = 8,
         main              = paste("CNV landscape per subclone —",
                                   sample_name),
         border_color      = NA)
dev.off()
cat("Saved:", cnv_heatmap_path, "\n")

#############################################
# CNV state per subclone per chromosome arm
#############################################
# Centromere positions (hg38)
# Used to classify each CNV region as p arm, q arm, or spanning both
# Positions sourced from UCSC hg38 centromere annotations
centromeres <- data.frame(
  chr        = paste0("chr", c(1:22, "X", "Y")),
  centromere = c(
    123400000, 93900000, 90900000, 50000000, 48800000, 59800000,
    60100000,  45200000, 43000000, 39800000, 53400000, 35500000,
    17700000,  17200000, 19000000, 36800000, 25100000, 18500000,
    26200000,  28100000, 12000000, 15000000, 61000000, 10400000
  )
)

# Add arm annotation to tumour CNV subclones
# Joins centromere positions then classifies each region:
#   end   <= centromere -> p arm (short arm)
#   start >= centromere -> q arm (long arm)
#   spanning both       -> p/q  (dropped below — ambiguous)
cnv_regions_ne <- cnv_regions_tumour %>%
  left_join(centromeres, by = "chr") %>%
  mutate(
    arm     = case_when(
      end   <= centromere ~ "p",
      start >= centromere ~ "q",
      TRUE                ~ "p/q"
    ),
    chr_arm = paste0(chr, arm)
  )

# Quick check — inspect arm assignments
cat("\nArm assignment preview:\n")
print(head(cnv_regions_ne %>%
             select(cell_group_name, chr, start, end,
                    state, arm, chr_arm)))

# Summarise dominant CNV state per subclone per arm
# Same logic as chromosome-level:
# Excludes neutral (3) then finds most frequent remaining state
arm_summary <- cnv_regions_ne %>%
  filter(state != 3) %>%
  group_by(cell_group_name, chr_arm) %>%
  summarise(
    dominant_state = as.character(names(which.max(table(state)))),
    .groups = "drop"
  )

cat("Arm summary:", nrow(arm_summary), "rows\n")

# Pivot to wide matrix
# Rows = subclones, columns = chromosome arms
# Arms with no non-neutral state filled as "3" (neutral)
arm_wide <- arm_summary %>%
  pivot_wider(
    names_from  = chr_arm,
    values_from = dominant_state,
    values_fill = "3")

# Order columns by genomic position
# chr1p, chr1q, chr2p, chr2q ... chr22p, chr22q, chrXp, chrXq
arm_order <- paste0("chr", rep(c(1:22, "X"), each = 2), c("p", "q"))
arm_order <- arm_order[arm_order %in% colnames(arm_wide)]

# Drop any p/q spanning columns — ambiguous, cannot assign to one arm
other_cols <- colnames(arm_wide)[
  !colnames(arm_wide) %in% c("cell_group_name", arm_order)]
if (length(other_cols) > 0) {
  message("Dropping centromere-spanning columns: ",
          paste(other_cols, collapse = ", "))
  arm_wide <- arm_wide %>% select(-all_of(other_cols))
}

# Convert to numeric matrix 
arm_mat <- arm_wide %>%
  column_to_rownames("cell_group_name") %>%
  select(all_of(arm_order)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# Fill any remaining NAs with 3 (neutral)
# Can occur if an arm has no CNV regions at all for a subclone
arm_mat[is.na(arm_mat)] <- 3

cat("Arm matrix:", nrow(arm_mat), "subclones x",
    ncol(arm_mat), "arms\n")

# Save arm_wide table
# Wide matrix of dominant CNV state per subclone per chromosome arm
arm_wide_path <- file.path(plot_dir,
  paste0(sample_name, "_cnv_arm_wide.csv"))
write.csv(arm_mat, file = arm_wide_path, row.names = FALSE)
cat("Saved:", arm_wide_path, "\n")

# Row annotation: tumour subtype
# Extracts subtype prefix from subclone name e.g. "NE1.s1.c2" -> "NE1"
row_ann_arm <- data.frame(
  NE_type   = gsub("\\..*", "", rownames(arm_mat)),
  row.names = rownames(arm_mat)
)

# NE subtype colours — generated dynamically
# Reads actual subtypes in arm matrix — may differ from chr matrix
ne_types_arm <- sort(unique(row_ann_arm$NE_type))
cat("Tumour subtypes in arm matrix:",
    paste(ne_types_arm, collapse = ", "), "\n")

if (length(ne_types_arm) <= 9) {
  ne_palette_arm <- RColorBrewer::brewer.pal(
    n    = max(length(ne_types_arm), 3),
    name = "Set1"
  )[seq_along(ne_types_arm)]
} else {
  ne_palette_arm <- colorRampPalette(
    RColorBrewer::brewer.pal(9, "Set1")
  )(length(ne_types_arm))
}

names(ne_palette_arm) <- ne_types_arm
ne_colors_arm <- list(NE_type = ne_palette_arm)

# Plot: arm-level CNV landscape heatmap
# Finer resolution than chromosome-level
# Each column is one chromosome arm (p or q)
# Helps identify focal events e.g. 17p loss vs 17q gain
arm_heatmap_path <- file.path(plot_dir,
  paste0(sample_name, "_cnv_arm_landscape_heatmap.pdf"))

pdf(arm_heatmap_path, width = 16, height = 8)
pheatmap(arm_mat,
         color             = cnv_colors,
         breaks            = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
         annotation_row    = row_ann_arm,
         annotation_colors = ne_colors_arm,
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         fontsize_col      = 7,
         main              = paste("CNV landscape per subclone — arm level |",
                                   sample_name),
         border_color      = NA)
dev.off()
cat("Saved:", arm_heatmap_path, "\n")

# Plot: subclone labels in tissue space 
# Shows spatial distribution of InferCNV subclones on tissue section
p_spatial_subclone <- ImageDimPlot(obj,
                                    group.by = "infercnv_subclone") +
  labs(title = paste("InferCNV subclones — spatial |", sample_name)) +
  theme(plot.title = element_text(size = 13, hjust = 0.5,
                                  face = "bold"))

spatial_subclone_path <- file.path(plot_dir,
  paste0(sample_name, "_infercnv_subclones_spatial.pdf"))
ggsave(spatial_subclone_path, p_spatial_subclone, width = 9, height = 7)
cat("Saved:", spatial_subclone_path, "\n")

# Save updated Seurat object
# Saves obj with two new metadata columns added by this script:
#   - cnv_score          : per-cell CNV variance score
#   - infercnv_subclone  : subclone assignment from InferCNV HMM
saveRDS(obj,
        file = paste0(sample_name, "_post_cnv.rds"))
cat("Saved:", paste0(sample_name, "_post_cnv.rds"), "\n")

cat("\nAll plots saved to:", plot_dir, "\n")
cat("CNV analysis complete for", sample_name, "\n")