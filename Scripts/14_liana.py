# 14_liana.py
# Purpose : Run LIANA+ cell-cell communication analysis on the
#           annotated spatial transcriptomics object. Computes
#           spatial neighbours, spatially variable genes (Moran's I),
#           ligand-receptor inflow scores, global specificity, and
#           generates chord diagrams and dotplots per NE subtype.
# Output  : results/liana/<sample>/
# Requires: sample name 
# ─────────────────────────────────────────────────────────────────

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import squidpy as sq
import liana as li
import matplotlib.pyplot as plt
import plotnine as p9
from plotnine import theme, element_text, element_blank, labs
from pycirclize import Circos

# Prompt: sample name
sample_name = input("\nEnter sample name: ").strip()
while not sample_name:
    print("Sample name cannot be empty.")
    sample_name = input("Enter sample name: ").strip()

# Prompt: h5ad file path
# AnnData object converted from the annotated Seurat object
while True:
    h5ad_path = input("\nEnter path to .h5ad file: ").strip()
    if not h5ad_path:
        print("Path cannot be empty.")
        continue
    if not os.path.exists(h5ad_path):
        print(f"File not found: '{h5ad_path}' — please check the path.")
        continue
    break

# Output directory
out_dir = os.path.join("results", "liana", sample_name)
os.makedirs(out_dir, exist_ok=True)
print(f"\nOutputs will be saved to: {out_dir}")

