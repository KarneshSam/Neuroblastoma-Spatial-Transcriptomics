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

# Prompt: custom LR pairs CSV
# Optional CSV with additional ligand-receptor pairs
# Must have columns: ligand, receptor
# Press Enter to skip and use only the consensus resource
while True:
    lr_csv = input(
        "\nEnter path to custom LR pairs CSV (or press Enter to skip): "
    ).strip()
    if lr_csv == "":
        lr_csv = None
        print("Skipping custom LR pairs — using consensus resource only.")
        break
    if not os.path.exists(lr_csv):
        print(f"File not found: '{lr_csv}' — please check the path.")
        continue
    break

# Prompt: tumour subtype prefix
tumour_prefix = input(
    "\nEnter tumour subtype prefix (e.g. NE): ").strip()

while not tumour_prefix:
    print("Prefix cannot be empty.")
    tumour_prefix = input("Enter tumour subtype prefix: ").strip()

# Prompt: spatial bandwidth
while True:
    bw_input = input(
        "\nEnter spatial bandwidth in µm [default: 200]: "
    ).strip()
    if bw_input == "":
        bandwidth = 200
        print("Using default: 200 µm")
        break
    try:
        bandwidth = float(bw_input)
        if bandwidth <= 0:
            raise ValueError
        break
    except ValueError:
        print("Invalid — please enter a positive number.")

# Prompt: lr_mean filter threshold
while True:
    lr_input = input(
        "\nEnter lr_mean filter threshold for dotplots [default: 0.01]: "
    ).strip()
    if lr_input == "":
        lr_threshold = 0.01
        print("Using default: 0.01")
        break
    try:
        lr_threshold = float(lr_input)
        if lr_threshold < 0:
            raise ValueError
        break
    except ValueError:
        print("Invalid — please enter a non-negative number.")

# Output directory
out_dir = os.path.join("results", "liana", sample_name)
os.makedirs(out_dir, exist_ok=True)
print(f"\nOutputs will be saved to: {out_dir}")

##############################################
# STEP 1 — Load AnnData
##############################################
print("\nLoading AnnData...")
adata = ad.read_h5ad(h5ad_path)
 
# Use log-normalised counts as the expression layer
adata.X = adata.layers["logcounts"].copy()
print(adata)
 
##############################################
# STEP 2 — Check spatial coordinates range
##############################################
print("\nSpatial coordinates range:")
print("X:", adata.obsm["spatial"][:, 0].min(),
      "to", adata.obsm["spatial"][:, 0].max())
print("Y:", adata.obsm["spatial"][:, 1].min(),
      "to", adata.obsm["spatial"][:, 1].max())
 
##############################################
# STEP 3 — Query bandwidth and visualise
##############################################
print("\nQuerying bandwidth options...")
plot, df = li.ut.query_bandwidth(
    coordinates=adata.obsm["spatial"],
    start=10,
    end=500,
    interval_n=40
)
 
# Save bandwidth plot with chosen bandwidth marked
bw_plot = plot + \
    p9.geom_vline(xintercept=bandwidth, linetype="dashed", color="red") + \
    p9.annotate("text", x=bandwidth, y=df.neighbours.max(),
                label=f"chosen={bandwidth}µm",
                angle=90, va="bottom", ha="right", color="red")
bw_plot.save(os.path.join(out_dir, f"{sample_name}_bandwidth_selection.pdf"),
             dpi=300)
print(f"Saved: bandwidth selection plot")

###############################################
# STEP 4 — Compute spatial neighbours
###############################################
print(f"\nComputing spatial neighbours (bandwidth={bandwidth} µm)...")
li.ut.spatial_neighbors(
    adata=adata,
    bandwidth=bandwidth,
    spatial_key="spatial"
)
print("Spatial neighbours computed.")

################################################
# STEP 5 — Spatially variable genes (Moran's I)
################################################
print("\nRunning Moran's I spatial autocorrelation...")
sq.gr.spatial_autocorr(adata, mode='moran', use_raw=False,
                        show_progress_bar=True)
 
# Filter to significant SVGs: BH-adjusted p < 0.05 and I > 0.01
svgs = adata.uns['moranI'].index[
    (adata.uns['moranI']['pval_norm_fdr_bh'] < 0.05) &
    (adata.uns['moranI']['I'] > 0.01)
]
print(f"Spatially variable genes (SVGs): {len(svgs)}")

# Save SVG table
adata.uns['moranI'].loc[svgs].sort_values("I", ascending=False).to_csv(
    os.path.join(out_dir, f"{sample_name}_svgs_moranI.csv"))
print("Saved: SVG Moran's I table")

 