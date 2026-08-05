library(Seurat)
library(infercnv)
library(rtracklayer)
library(dplyr)

##################################
# Load the list of Seurat object
##################################
rds_obj <- readRDS("seu_enact_seg_unprocessed.rds")

# Setup: The file can be a single Seurat object or a list of Seurat objects of any length
# Case 1: single Seurat object — wrap in list for uniform handling
# Case 2: named or unnamed list of Seurat objects of any length
if (inherits(rds_obj, "Seurat")) {
  cat("Detected: single Seurat object\n")
  sample_list <- list(S1 = rds_obj)          # assign a name S1 to the single object for downstream use

} 
else if (is.list(rds_obj) &&
           all(sapply(rds_obj, inherits, "Seurat"))) {
  cat("Detected:", length(rds_obj), "Seurat objects in list\n")

  # Use existing names if available, otherwise assign S1, S2, S3 ...
  if (is.null(names(rds_obj))) {
    names(rds_obj) <- paste0("S", seq_along(rds_obj))
    cat("No names found — assigned:",
        paste(names(rds_obj), collapse = ", "), "\n")
  }
  sample_list <- rds_obj               # assign the list to sample_list for downstream use

} 
else {
  stop("The file does not contain a Seurat object or a list of Seurat objects.")
}

# Print summary of all samples
cat("\nSamples available:\n")
for (nm in names(sample_list)) {
  cat(sprintf("  %-6s | Cells: %d | Genes: %d\n",
              nm, ncol(sample_list[[nm]]), nrow(sample_list[[nm]])))
}

# Prompt: output directory 
# All InferCNV outputs including gene order file go here
repeat {
  out_dir <- trimws(readline(prompt = "\nEnter output directory name: "))
  if (nchar(out_dir) == 0) {
    cat("Output directory cannot be empty.\n"); next
  }
  break
}

# Create output directory if it does not exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("Created output directory:", out_dir, "\n")
}

#########################
# Prompt: choose sample 
#########################
# Provide the available sample names to the user and ask which one to run InferCNV on
repeat {
  cat("\nAvailable samples:", paste(names(sample_list), collapse = ", "), "\n")
  sample_name <- trimws(readline(prompt = "Which sample to run InferCNV on? "))
  if (sample_name %in% names(sample_list)) break
  cat("Invalid '", sample_name,
      "' — please enter one of the listed names.\n", sep = "")
}

# Load that Seurat object into the obj variable
obj <- sample_list[[sample_name]]
cat("\nLoaded:", sample_name,
    "| Cells:", ncol(obj),
    "| Genes:", nrow(obj), "\n")

# Show available cell types
all_cell_types <- sort(unique(obj$cell_type_hr))
cat("\nCell types available in", sample_name, ":\n")
for (i in seq_along(all_cell_types)) {
  cat(sprintf("  [%d] %s\n", i, all_cell_types[i]))
}

#############################
# Prompt: tumour cell types 
##############################
# These are the populations InferCNV will test for CNV signal
repeat {
  cat("\nEnter tumour cell type names, comma-separated.\n")
  cat("These will be tested for CNV signal.\n")
  tumour_input <- trimws(readline(prompt = "Tumour types: "))
  tumour_types <- trimws(strsplit(tumour_input, ",")[[1]])
  tumour_types <- tumour_types[nzchar(tumour_types)]

  # If any of the entered names are not found in the object, warn and re-ask 
  invalid <- tumour_types[!tumour_types %in% all_cell_types]
  if (length(tumour_types) == 0) {
    cat("Please enter at least one cell type.\n"); next
  }
  if (length(invalid) > 0) {
    cat("Not found in object:", paste(invalid, collapse = ", "), "\n")
    cat("Please enter names exactly as listed above.\n"); next
  }
  break
}
cat("Tumour types:", paste(tumour_types, collapse = ", "), "\n")

######################################
# Prompt: normal reference cell types
######################################
# Used as the reference baseline for CNV inference
repeat {
  cat("\nEnter normal reference cell type names, comma-separated.\n")
  cat("These will be used as the diploid reference baseline.\n")
  normal_input <- trimws(readline(prompt = "Normal types: "))
  normal_types <- trimws(strsplit(normal_input, ",")[[1]])
  normal_types <- normal_types[nzchar(normal_types)]

  invalid <- normal_types[!normal_types %in% all_cell_types]
  if (length(normal_types) == 0) {
    cat("Please enter at least one cell type.\n"); next
  }
  if (length(invalid) > 0) {
    cat("Not found in object:", paste(invalid, collapse = ", "), "\n")
    cat("Please enter names exactly as listed above.\n"); next
  }

  # Warn if any normal type overlaps with tumour types
  overlap <- intersect(tumour_types, normal_types)
  if (length(overlap) > 0) {
    cat("Warning: these types appear in both tumour and normal:",
        paste(overlap, collapse = ", "), "\n")
    cat("This will cause issues. Please re-enter normal types.\n"); next
  }
  break
}
cat("Normal types:", paste(normal_types, collapse = ", "), "\n")

#####################################
# Generate gene order file from GTF
#####################################
# Gene order must match the genes actually present in the object
# so it is generated fresh per sample rather than using a static file
# Only standard chromosomes (chr1-22, X, Y) are kept
# Duplicate gene names are removed keeping the first genomic occurrence
cat("\nImporting GTF file — this may take a moment...\n")
gtf <- import("gencode.v44.basic.annotation.gtf.gz")

cat("Building gene order file for", sample_name, "...\n")
genes_df <- as.data.frame(gtf) |>
  # Keep gene-level entries only (not transcript, exon etc.)
  dplyr::filter(type == "gene") |>
  # Select only the columns InferCNV needs: name, chr, start, end
  dplyr::select(gene_name, seqnames, start, end) |>
  # Keep only standard chromosomes
  dplyr::filter(seqnames %in% paste0("chr", c(1:22, "X", "Y"))) |>
  # Keep only genes present in the Seurat object
  dplyr::filter(gene_name %in% rownames(obj)) |>
  # Sort by chromosome then position for correct CNV ordering
  dplyr::arrange(seqnames, start) |>
  # Remove duplicate gene names — keep first genomic occurrence
  dplyr::distinct(gene_name, .keep_all = TRUE)

# Summary check
cat("Genes in object:", nrow(obj), "\n")
cat("Genes matched in GTF:", nrow(genes_df), "\n")
cat("Genes not found in GTF:",
    length(setdiff(rownames(obj), genes_df$gene_name)), "\n")

# Write gene order file to output directory
gene_order_path <- file.path(out_dir, "gene_order_file.txt")
write.table(genes_df,
            file      = gene_order_path,
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Verify no duplicates — should always be 0 after distinct()
gene_order_check <- read.table(gene_order_path, sep = "\t")
n_dups <- sum(duplicated(gene_order_check$V1))
if (n_dups > 0) {
  warning("Duplicate gene names found in gene order file: ", n_dups)
} else {
  cat("Gene order file verified — no duplicates\n")
}
cat("Gene order file saved:", gene_order_path, "\n")

###############################
# Prepare InferCNV input files
###############################
# Get raw counts
# InferCNV requires raw integer counts, not normalised values
counts_mat <- GetAssayData(obj,
                           assay = "Spatial.Polygons",
                           layer = "counts")

# Build cell annotation data frame
# Two columns: cell barcode and cell type label
# Only tumour and normal cells are kept 
cell_ann <- data.frame(
  cell_barcode = colnames(counts_mat),
  cell_type    = obj$cell_type_hr[
    match(colnames(counts_mat), colnames(obj))])

# Filter to tumour and normal cells only
cell_ann_filt <- cell_ann[
  cell_ann$cell_type_hr %in% c(tumour_types, normal_types), ]

# Write annotation file 
# InferCNV reads annotations from a tab-separated text file
ann_file <- file.path(out_dir, "cell_annotations.txt")
write.table(cell_ann_filt,
            file      = ann_file,
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE,
            col.names = FALSE)
cat("Annotation file written:", ann_file, "\n")

# Subset counts to annotated cells only 
# Removes any cells not in the filtered annotation table
counts_filt <- counts_mat[,
                          cell_ann_filt$cell_barcode]

cat("Counts matrix dimensions:",
    nrow(counts_filt), "genes x",
    ncol(counts_filt), "cells\n")

#########################
# Create InferCNV object
#########################
# ref_group_names: normal cell types 
# min_max_counts_per_cell = c(0, Inf): include all cells 
# since filtering was already done upstream in the Seurat pipeline
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix       = counts_filt,
  annotations_file        = ann_file,
  delim                   = "\t",
  gene_order_file         = gene_order_path,
  ref_group_names         = normal_types,
  min_max_counts_per_cell = c(0, Inf))

# Check mean expression per gene
# Helps choose a sensible cutoff value for infercnv::run
# cutoff removes lowly expressed genes before CNV inference
expr_means <- rowMeans(infercnv_obj@expr.data)

cat("\nMean expression per gene summary:\n")
print(summary(expr_means))

cat("\nGenes above cutoff thresholds:\n")
# Select a range of cutoff values to show how many genes would be included at each threshold
for (thresh in c(0.01, 0.05, 0.1, 0.2, 0.5)) {
  n_genes <- sum(expr_means >= thresh)
  cat(sprintf("  cutoff >= %.2f  ->  %d genes (%.1f%%)\n",
              thresh,
              n_genes,
              n_genes / length(expr_means) * 100))
}

