# 00_setup.R
# Purpose : Install ALL system-level and R package dependencies
#           required by the full pipeline.
#           Run this ONCE from the terminal before starting:
#
#             Rscript 00_setup.R
#
#           After this completes, all scripts 01_qc.R through
#           12_gsea_results.R will load without errors.
# ─────────────────────────────────────────────────────────────────

cat("==========================================================\n")
cat(" Pipeline Dependency Installer\n")
cat("==========================================================\n\n")

# STEP 1 — System-level dependencies (Ubuntu / Debian)
# These must be installed before any R packages can compile.
# Requires sudo — if you do not have sudo, ask your sysadmin to
# run this block, then proceed from Step 2.

cat("Step 1: System libraries (requires sudo)\n")
cat("Checking OS...\n")

os_info <- tryCatch(readLines("/etc/os-release"), error = function(e) "")

if (any(grepl("Ubuntu|Debian", os_info))) {

  cat("Ubuntu/Debian detected — installing system libraries...\n\n")

  sys_packages <- c(
    # curl / SSL / XML (Seurat, rtracklayer, BiocManager)
    "libcurl4-openssl-dev",
    "libssl-dev",
    "libxml2-dev",

    # HDF5 (hdf5r, used by Seurat for .h5 files) 
    "libhdf5-dev",

    # Graphics / font rendering (ggplot2, patchwork)
    "libfontconfig1-dev",
    "libharfbuzz-dev",
    "libfribidi-dev",
    "libfreetype6-dev",
    "libpng-dev",
    "libtiff5-dev",
    "libjpeg-dev",
    "libcairo2-dev",
    "libxt-dev",

    # Compression (rtracklayer, Bioconductor I/O)
    "libbz2-dev",
    "liblzma-dev",
    "zlib1g-dev",

    # Graph / linear algebra (igraph, leiden, Matrix) 
    "libglpk-dev",
    "libgmp-dev",
    "liblapack-dev",
    "libblas-dev",

    # GSL (infercnv dependencies) 
    "libgsl-dev",

    # CMake (needed to compile some C++ R packages) 
    "cmake",

    # Git (gert, credentials, devtools) 
    "libgit2-dev",

    # Units / spatial (pulled by tidyverse dependencies) 
    "libudunits2-dev",
    "libgdal-dev",
    "libproj-dev",

    # Java (AnnotationDbi, org.Hs.eg.db) 
    "default-jdk",
    "default-jre",

    # Python (leiden uses leidenalg Python package)
    "python3",
    "python3-pip",
    "python3-dev"
  )

  for (pkg in sys_packages) {
    ret <- system(paste("dpkg -s", pkg, "> /dev/null 2>&1"))
    if (ret != 0) {
      cat("  Installing system package:", pkg, "...\n")
      system(paste("sudo apt-get install -y", pkg,
                   "2>&1 | tail -1"),
             ignore.stdout = FALSE)
    } else {
      cat(sprintf("  %-35s already installed\n", pkg))
    }
  }

  # Reconfigure Java for R (needed by rJava / AnnotationDbi)
  cat("\nReconfiguring Java for R (sudo R CMD javareconf)...\n")
  system("sudo R CMD javareconf 2>&1 | tail -3")

  # Install leidenalg Python package (required by leiden R package)
  cat("\nInstalling Python leidenalg (required by leiden R package)...\n")
  system("pip3 install leidenalg igraph --quiet")

} else if (any(grepl("CentOS|Red Hat|Fedora|Rocky", os_info))) {

  cat("CentOS/RHEL/Fedora detected.\n")
  cat("Please install system dependencies manually:\n\n")
  cat("  sudo yum install -y \\\n")
  cat("    openssl-devel libcurl-devel libxml2-devel hdf5-devel \\\n")
  cat("    fontconfig-devel freetype-devel libpng-devel libtiff-devel \\\n")
  cat("    libjpeg-devel cairo-devel bzip2-devel xz-devel zlib-devel \\\n")
  cat("    glpk-devel gmp-devel lapack-devel blas-devel gsl-devel \\\n")
  cat("    cmake git java-11-openjdk-devel python3 python3-pip\n\n")
  cat("Then re-run this script.\n")

} else {
  cat("OS not recognised as Ubuntu/Debian or RHEL.\n")
  cat("Please install system dependencies manually before proceeding.\n")
  cat("See the README for the full list.\n")
}

# STEP 2 — Bootstrap R installers

cat("\nStep 2: Bootstrap R installers\n")

install_if_missing <- function(pkg, source = "CRAN", repo = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing:", pkg, "from", source, "...\n")
    tryCatch({
      if (source == "CRAN") {
        install.packages(pkg,
                         repos        = "https://cloud.r-project.org",
                         dependencies = TRUE,
                         quiet        = TRUE)
      } else if (source == "Bioconductor") {
        BiocManager::install(pkg,
                             update = FALSE,
                             ask    = FALSE,
                             quiet  = TRUE)
      }
      if (requireNamespace(pkg, quietly = TRUE)) {
        cat("    OK\n")
      } else {
        cat("    FAILED — package did not install correctly\n")
      }
    }, error = function(e) {
      cat("    FAILED:", conditionMessage(e), "\n")
    })
  } else {
    cat(sprintf("  %-30s already installed\n", pkg))
  }
}

install_if_missing("BiocManager", "CRAN")
install_if_missing("remotes",     "CRAN")

# STEP 3 — CRAN packages

cat("\nStep 3: CRAN packages\n")

cran_packages <- c(
  # Core data wrangling
  "dplyr",        # Data manipulation
  "tidyr",        # Reshaping (pivot_wider, pivot_longer)
  "tibble",       # Modern data frames
  "stringr",      # String manipulation
  "tidyverse",    # Meta-package (ggplot2, dplyr, tidyr etc.)
  "scales",       # Scale functions for ggplot2

  # Plotting
  "ggplot2",      # Core plotting
  "ggrepel",      # Non-overlapping text labels
  "patchwork",    # Combining multiple ggplots
  "pheatmap",     # Heatmaps
  "RColorBrewer", # Colour palettes
  "ggh4x",        # Extended ggplot2 axes (guide_axis_color)

  # Bioinformatics utilities
  "Matrix",       # Sparse matrix operations
  "leiden",       # Leiden clustering algorithm
  "RANN",         # Fast nearest neighbour search

  # Seurat ecosystem
  "Seurat",       # Core single-cell / spatial analysis
  "clustree"      # Cluster resolution visualisation
)

for (pkg in cran_packages) {
  install_if_missing(pkg, "CRAN")
}

# STEP 4 — Bioconductor packages
cat("\nStep 4: Bioconductor packages\n")

bioc_packages <- c(
  "Banksy",          # Spatial neighbourhood embedding
  "SeuratWrappers",  # RunBanksy and other Seurat extensions
  "infercnv",        # Copy number variation inference
  "rtracklayer",     # Import GTF / GFF files
  "clusterProfiler", # Gene set enrichment analysis framework
  "ReactomePA",      # Reactome pathway analysis
  "org.Hs.eg.db",   # Human gene annotation (Entrez ID conversion)
  "AnnotationDbi",   # Database interface (dependency of org.Hs.eg.db)
  "BiocParallel"     # Parallel processing (used by clusterProfiler)
)

for (pkg in bioc_packages) {
  install_if_missing(pkg, "Bioconductor")
}

