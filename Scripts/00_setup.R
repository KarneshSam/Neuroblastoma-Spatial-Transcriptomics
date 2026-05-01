# 00_setup.R
# ─────────────────────────────────────────────────────────────────

library(Seurat)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)
library(leiden)
library(Banksy)
library(pheatmap)
library(SeuratWrappers)
library(clustree)
library(Matrix)

seu <- readRDS("seu_enact_seg_unprocessed.rds")

P1 <- seu[[1]]
P1.B <- seu[[2]]
P2 <- seu[[3]]
P2.B <- seu[[4]]

########################
# Prompt: choose sample
########################
repeat {
  cat("\nAvailable samples: P1, P1.B, P2, P2.B\n")
  sample_name <- trimws(readline(
    prompt = "Which sample to run? (P1 / P1.B / P2 / P2.B): "))
  if (sample_name %in% c("P1", "P1.B", "P2", "P2.B")) break
  cat("Invalid input '", sample_name, "' — please enter P1, P1.B, P2, or P2.B.\n", sep = "")
}

