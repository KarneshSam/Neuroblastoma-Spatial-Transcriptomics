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

# Step 3: Run for all detected NE subtypes
neighbours_all <- map_dfr(ne_types, function(nt) {
  cat("Computing", nt, "neighbourhood...\n")
  get_neighbourhood(nt, radius = radius)
})
 
cat("\nNeighbourhood computation done.\n")

# Step 4: Show available cell types
# Lists all unique neighbour cell types found in the data
# User assigns each to a broad class and defines display order
all_types <- sort(unique(neighbours_all$neighbour_type))

cat("\nNeighbour cell types found in", sample_name, ":\n")
for (i in seq_along(all_types)) {
  cat(sprintf("  [%d] %s\n", i, all_types[i]))
}

# Prompt: assign each cell type to a broad class
# Available classes: Tumour, Stromal, Endothelial, Immune, Ambiguous
valid_classes <- c("Tumour", "Stromal", "Endothelial", "Immune", "Ambiguous")

cat("\nAssign each cell type to a broad class.\n")
cat("Options:", paste(valid_classes, collapse = ", "), "\n\n")

class_map <- c()

for (ct in all_types) {
  repeat {
    cls <- trimws(readline(
      prompt = paste0("  ", ct, " -> class: ")))
    if (nchar(cls) == 0) {
      cat("  Cannot be empty.\n"); next
    }
    # Case-insensitive match
    matched <- valid_classes[tolower(valid_classes) == tolower(cls)]
    if (length(matched) == 0) {
      cat("  Invalid — choose from:",
          paste(valid_classes, collapse = ", "), "\n"); next
    }
    class_map[ct] <- matched
    break
  }
}

cat("\nClass assignments:\n")
for (ct in names(class_map)) {
  cat(sprintf("  %-30s -> %s\n", ct, class_map[ct]))
}

# Prompt: display order
# User defines the order of cell types on the x-axis
# Enter as comma-separated list — types not listed are dropped
cat("\nDefine the display order for the x-axis.\n")
cat("Enter cell type names comma-separated in the order you want.\n")
cat("Available:\n")
for (ct in all_types) cat("  ", ct, "\n")

repeat {
  order_input <- trimws(readline(prompt = "\nX-axis order: "))
  type_order  <- trimws(strsplit(order_input, ",")[[1]])
  type_order  <- type_order[nzchar(type_order)]

  # Warn about any names not in the data but still accept
  unknown <- type_order[!type_order %in% all_types]
  if (length(unknown) > 0) {
    cat("Warning: these names not found in data and will be ignored:",
        paste(unknown, collapse = ", "), "\n")
  }
  type_order <- type_order[type_order %in% all_types]
  if (length(type_order) == 0) {
    cat("No valid cell types entered. Please try again.\n"); next
  }
  break
}

# Colour palette by broad cell class
class_colours <- c(
  "Tumour"      = "#D62728",
  "Immune"      = "#1F77B4",
  "Stromal"     = "#2CA02C",
  "Endothelial" = "#FF7F0E",
  "Ambiguous"   = "#BDBDBD"
)

# Step 5: Build plot data
# Assign each neighbour type to its user-defined broad cell class
# Factor levels for x-axis follow user-defined type_order
plot_data <- neighbours_all %>%
  mutate(
    cell_class     = class_map[neighbour_type],
    neighbour_type = factor(neighbour_type, levels = type_order),
    cell_class     = factor(cell_class,
                            levels = c("Tumour", "Stromal",
                                       "Endothelial", "Immune",
                                       "Ambiguous")),
    target       = factor(target, levels = ne_types)
  ) %>%
  filter(!is.na(cell_class))

# Step 6: Plot
p <- ggplot(plot_data,
            aes(x    = neighbour_type,
                y    = frequency,
                fill = cell_class)) +
 
  # Stacked bar per neighbour cell type
  geom_bar(stat = "identity", width = 0.75) +
 
  # Label bars where frequency >= 0.03
  geom_text(
    data     = plot_data %>% filter(frequency >= 0.03),
    aes(label = sprintf("%.2f", frequency)),
    vjust    = -0.4,
    size     = 2.6,
    colour   = "grey20",
    fontface = "bold"
  ) +
 
  # One panel per NE subtype
  facet_wrap(~ target, ncol = 3, scales = "free_y",
             axes = "all_x", axis.labels = "all_x") +
 
  scale_fill_manual(values = class_colours, name = "Cell types") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.22)),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  scale_x_discrete(drop = FALSE) +
 
  # Colour x-axis tick labels by cell class
  guides(
    x = guide_axis_color(
      color = class_colours[
        as.character(
          plot_data %>%
            distinct(neighbour_type, cell_class) %>%
            arrange(neighbour_type) %>%
            pull(cell_class)
        )
      ]
    ),
    fill = guide_legend(ncol = 2)
  ) +
 
  theme_classic(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold", size = 12,
                                    colour = "white"),
    strip.background = element_rect(fill = "grey30", colour = NA),
 
    axis.text.x  = element_text(angle = 90, hjust = 1,
                                vjust = 0.5, size = 8),
    axis.text.y  = element_text(size = 9),
    axis.title   = element_text(size = 11),
    axis.line    = element_line(colour = "grey60"),
    axis.ticks   = element_line(colour = "grey60"),
 
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
 
    legend.position        = "inside",
    legend.position.inside = c(0.45, 0.10),
    legend.direction       = "vertical",
    legend.title           = element_text(face = "bold", size = 10),
    legend.text            = element_text(size = 9),
    legend.key.size        = unit(0.45, "cm"),
 
    panel.spacing = unit(0.8, "lines"),
    plot.title    = element_text(face = "bold", size = 13,
                                 hjust = 0.5),
    plot.caption  = element_text(size = 9, colour = "grey40",
                                 hjust = 0),
    plot.margin   = margin(10, 15, 10, 10)
  ) +
 
  labs(
    title   = paste("Neighbourhood composition —", sample_name),
    x       = NULL,
    y       = "Proportion of neighbours",
    caption = paste0("*Spatial radius r = ", radius,
                     " µm  |  Bars labelled when >= 0.03")
  ) 

# Save plot
out_path <- file.path(out_dir,
  paste0(sample_name, "_neighbourhood_composition.pdf"))
ggsave(out_path, p,
       width  = max(14, length(ne_types) * 5),
       height = 10,
       units  = "in")
cat("Saved:", out_path, "\n")

# Save neighbourhood frequency table
# Raw frequency table for further analysis if needed
csv_path <- file.path(out_dir,
  paste0(sample_name, "_neighbourhood_frequencies.csv"))
write.csv(neighbours_all, csv_path, row.names = FALSE)
cat("Saved:", csv_path, "\n")
 
cat("\nNeighbourhood analysis complete for", sample_name, "\n")