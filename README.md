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
