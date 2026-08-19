# Exploring Spatial Intratumour Heterogeneity in Neuroblastoma

## Overview

A modular R pipeline for processing, clustering, annotating, and analysing
spatial transcriptomics data from neuroblastoma tumour tissue sections,
integrating CNV inference, GSEA, spatial neighbourhood analysis, and
cell-cell communication.

---

## Purpose of the Analysis

- Quality control and filtering of spatially resolved single-cell data
- Spatially-informed dimensionality reduction and clustering using BANKSY
- Cell type annotation based on marker genes
- CNV inference to identify tumour subclones and chromosomal alterations
- GSEA to characterise biological programmes enriched per CNV-defined subclone
- Spatial neighbourhood analysis to quantify the tumour microenvironment composition
- Cell-cell communication analysis using LIANA+ (Python)

---

## Background

Neuroblastoma (NB) is the most common extracranial solid malignancy in children,
characterised by chromosomal copy number variations (CNVs), transcriptional
plasticity between adrenergic (ADRN) and mesenchymal (MES) cell states, and a
complex tumour microenvironment (TME) including tumour-associated macrophages
(TAMs) and cancer-associated fibroblasts (CAFs). Standard spatial transcriptomics
platforms aggregate signal across multiple cells; Visium HD provides 2 µm bin
resolution enabling near single-nucleus profiling of intact tissue sections. This
pipeline integrates CNV analysis with spatial context to examine how the
surrounding microenvironment shapes clonal dynamics and tumour heterogeneity.

---

## Data Description

### Input Data

A Seurat object containing the cell × gene count matrix output from SpaceRanger
and ENACT cell segmentation, stored under the assay name `Spatial.Polygons`.
The file may contain a single Seurat object or a named list of Seurat objects —
the pipeline handles both automatically.

---

## Directory Structure

```
.
├── data/
│   ├── seu_enact_seg_unprocessed.rds         # Input Seurat RDS (not available in the github)
│   └── gencode.v44.basic.annotation.gtf.gz   # Gencode GTF for InferCNV
│
├── scripts/
│   ├── 00_setup.R                             # Install all dependencies
│   ├── 01_qc.R
│   ├── 02_filter.R
│   ├── 03_normalise.R
│   ├── 04_banksy_pca_umap.R
│   ├── 05a_cluster_clustree.R
│   ├── 05b_plots.R
│   ├── 06_markers.R
│   ├── 07_annotation.R
│   ├── 08_infercnv_pipeline.R
│   ├── 09_cnv_downstream_analysis.R
│   ├── 10_cnv_metadata.R
│   ├── 11_gsea_subclone.R
│   ├── 12_gsea_results.R
│   ├── 13_neighbourhood.R
│   └── 14_liana.py
│
└── results/
    ├── tmp/
    │   └── <sample_name>/                     # Intermediate checkpoints (01–06)
    │       ├── <sample>_post_qc.rds
    │       ├── <sample>_post_filter.rds
    │       ├── <sample>_post_normalise.rds
    │       ├── <sample>_post_banksy.rds
    │       ├── <sample>_post_cluster.rds
    │       └── <sample>_final.rds
    │
    ├── <sample_name>/
    │     └──<sample>_annotated.rds   # Final annotated object (07)
    │
    ├── <infercnv_output_dir>/                 # InferCNV outputs (08)
    │    └── infercnv_<sample_name>/ 
    │        ├── run.final.infercnv_obj
    │        ├── gene_order_file.txt
    │        ├── cell_annotations.txt
    │        └── 17_HMM_predHMMi6.leiden.*
    │
    ├── <cnv_downstream_output_dir>/           # CNV downstream plots (09)
    │   ├── <sample>_cnv_chr_wide.csv
    │   ├── <sample>_cnv_arm_wide.csv
    │   ├── <sample>_cnv_landscape_heatmap.pdf
    │   ├── <sample>_cnv_arm_landscape_heatmap.pdf
    │   ├── <sample>_infercnv_subclones_spatial.pdf
    │   └── <sample>_post_cnv.rds
    │
    ├── subclone/
    │   └── <sample_name>/                     # CNV metadata outputs (10)
    │       ├── <sample>_post_cnv_meta.rds
    │       ├── <sample>_ne_cells.rds
    │       ├── <sample>_cnv_chr_state_barplot.pdf
    │       ├── <sample>_cnv_arm_state_barplot.pdf
    │       ├── <sample>_cnv_directionality_all_arms.pdf
    │       └── <sample>_cnv_directionality_classic_arms.pdf
    │
    ├── gsea/
    │   └── <sample_name>/                     # GSEA outputs (11 + 12)
    │       ├── avg_expression_per_subclone.csv
    │       ├── gsea_all_subclones.csv
    │       ├── gsea_activated_subclones.csv
    │       ├── gsea_suppressed_subclones.csv
    │       ├── <sample>_gsea_pathway_heatmap.pdf
    │       ├── <sample>_gsea_bubble_plot.pdf
    │       ├── <sample>_cnv_gsea_pathway_bubble.pdf
    │       └── <sample>_<NE>_gsea_core_enrichment.png
    │
    ├── neighbour/
    │   └── <sample>_neighbourhood_composition.png   # Neighbourhood plots (13)
    │
    └── liana/
        └── <sample_name>/                     # LIANA+ outputs (14)
            ├── chord_diagram_overall.png
            ├── chord_diagram_<cell_type>.png
            └── lr_inflow_scores.pdf
```

---

## Required Software

### R Pipeline (scripts 00–13)

| Package | Source | Purpose |
|---------|--------|---------|
| Seurat | CRAN | Core single-cell / spatial analysis |
| SeuratWrappers | Bioconductor | RunBanksy extension |
| Banksy | Bioconductor | Spatial neighbourhood embedding |
| clustree | CRAN | Cluster resolution visualisation |
| leiden | CRAN | Leiden clustering algorithm |
| Matrix | CRAN | Sparse matrix operations |
| infercnv | Bioconductor | Copy number variation inference |
| rtracklayer | Bioconductor | GTF file import |
| clusterProfiler | Bioconductor | GSEA framework |
| ReactomePA | Bioconductor | Reactome pathway analysis |
| org.Hs.eg.db | Bioconductor | Human gene ID conversion |
| ggplot2 / ggrepel / patchwork | CRAN | Plotting |
| pheatmap / RColorBrewer / ggh4x | CRAN | Heatmaps and colour palettes |
| RANN | CRAN | Spatial nearest neighbour search |
| dplyr / tidyr / tibble / stringr / tidyverse / scales | CRAN | Data wrangling |

### Python Pipeline (script 14)

| Package | Purpose |
|---------|---------|
| liana 1.2 | Cell-cell communication analysis |
| squidpy 1.8.3 | Spatially variable gene detection |
| anndata | AnnData object handling |
| pyCirclize | Chord diagram visualisation |

---

## Quick Start

### 1. Install all R dependencies (run once)

```bash
Rscript scripts/00_setup.R
```

Installs all system-level (Ubuntu/Debian) and R package dependencies and
verifies they load correctly.

### 2. Run the pipeline

Run all scripts from the **project root directory**:

```bash
Rscript scripts/01_qc.R
Rscript scripts/02_filter.R
Rscript scripts/03_normalise.R
Rscript scripts/04_banksy_pca_umap.R
Rscript scripts/05a_cluster_clustree.R
Rscript scripts/05b_plots.R
Rscript scripts/06_markers.R
Rscript scripts/07_annotation.R
Rscript scripts/08_infercnv_pipeline.R
Rscript scripts/09_cnv_downstream_analysis.R
Rscript scripts/10_cnv_metadata.R
Rscript scripts/11_gsea_subclone.R
Rscript scripts/12_gsea_results.R
Rscript scripts/13_neighbourhood.R
python scripts/14_liana.py
```

Each script prompts you interactively for inputs — input file, sample
selection, thresholds, and output directories. No hardcoded paths. All
threshold prompts show a default value; press **Enter** to accept it.

---

## Workflow

### Stage 1 — Preprocessing (scripts 01–07)

| Script | What it does |
|--------|-------------|
| `01_qc.R` | Loads RDS, selects sample, computes mitochondrial %, generates QC violin plots and gene count distribution |
| `02_filter.R` | Iteratively filters low-quality cells (gene count, spatial bins, MT%) and lowly expressed genes until convergence |
| `03_normalise.R` | LogNormalizes counts (scale factor 10,000), identifies top 2,000 HVGs using VST, scales data regressing out MT% |
| `04_banksy_pca_umap.R` | Runs BANKSY spatial embedding (λ=0.2, k=50 neighbours), PCA (30 PCs), kNN graph, and UMAP |
| `05a_cluster_clustree.R` | Runs Leiden clustering at 5 user-defined resolutions; plots clustree to assess stability |
| `05b_plots.R` | Generates UMAP, spatial ImageDimPlot, and SNN connectivity heatmap for the chosen resolution |
| `06_markers.R` | Finds marker genes per cluster (Wilcoxon, only.pos, min.pct=0.25, logFC≥0.25, adj p<0.05), plots DotPlot |
| `07_annotation.R` | Interactively labels each cluster with a cell type name; saves final annotated object to `results/annotated_obj/` |

### Stage 2 — CNV Inference (scripts 08–10)

| Script | What it does |
|--------|-------------|
| `08_infercnv_pipeline.R` | Generates gene order file from Gencode GTF; runs InferCNV with HMM state calling (6 states), denoising, and user-defined cutoff |
| `09_cnv_downstream_analysis.R` | Computes per-cell CNV variance scores, assigns HMM subclone labels, summarises dominant CNV state per chromosome and arm, generates heatmaps |
| `10_cnv_metadata.R` | Adds chromosome and arm-level CNV states to metadata, computes directionality scores (proportion gain − proportion loss), subsets NE cells |

### Stage 3 — GSEA (scripts 11–12)

| Script | What it does |
|--------|-------------|
| `11_gsea_subclone.R` | Computes average expression per subclone, converts to Entrez IDs, runs gsePathway (Reactome) ranked by expression per subclone |
| `12_gsea_results.R` | Filters by p.adjust < 0.05, applies user-defined NES threshold, generates pathway heatmap and bubble plot, overlaps core enrichment genes with InferCNV CNV genes, produces per-NE combined plots |

### Stage 4 — Spatial Analysis (scripts 13–14)

| Script | What it does |
|--------|-------------|
| `13_neighbourhood.R` | Fixed-radius nearest-neighbour search (RANN, r=200 µm, k=100) per NE subtype; computes neighbour type proportions; plots stacked bar charts coloured by cell class |
| `14_liana.py` | Builds spatial neighbourhood graphs, detects spatially variable genes (Moran's I), computes ligand-receptor inflow scores, generates chord diagrams  |

---

## Output Files

### `results/tmp/<sample_name>/` — Intermediate checkpoints

| File | Description |
|------|-------------|
| `<sample>_post_qc.rds` | Object with percent.mt added |
| `<sample>_post_filter.rds` | Filtered object |
| `<sample>_post_normalise.rds` | Normalised and scaled object |
| `<sample>_post_banksy.rds` | BANKSY + PCA + UMAP + kNN |
| `<sample>_post_cluster.rds` | All cluster resolutions added |
| `<sample>_final.rds` | Single chosen resolution, clean |

### `results/annotated_obj/` — Final object

| File | Description |
|------|-------------|
| `<sample>_annotated.rds` | Final annotated Seurat object with `cell_type_hr` column |

