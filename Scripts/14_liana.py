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

#################################################
# STEP 6 — Load LR resource
#################################################
print("\nLoading LR resource (consensus)...")
resource = li.rs.select_resource('consensus')
print(f"Consensus resource: {len(resource)} LR pairs")
 
# Optional: load and merge custom LR pairs
if lr_csv is not None:
    print(f"\nLoading custom LR pairs from: {lr_csv}")
    your_pasted_df = pd.read_csv(lr_csv)
    print(your_pasted_df.head())
 
    # Extract unique ligand-receptor pairs
    custom_pairs = your_pasted_df[['ligand', 'receptor']].drop_duplicates().reset_index(drop=True)
    print(f"Unique custom LR pairs: {len(custom_pairs)}")
 
    # Check which pairs have genes present in the panel
    def gene_in_panel(pair_str, panel_genes):
        return all(g in panel_genes for g in pair_str.split('_'))
 
    panel_genes = set(adata.var_names)
    custom_pairs['ligand_in_panel']   = custom_pairs['ligand'].apply(
        lambda x: gene_in_panel(x, panel_genes))
    custom_pairs['receptor_in_panel'] = custom_pairs['receptor'].apply(
        lambda x: gene_in_panel(x, panel_genes))
 
    recoverable  = custom_pairs[
        custom_pairs['ligand_in_panel'] & custom_pairs['receptor_in_panel']]
    truly_missing = custom_pairs[
        ~(custom_pairs['ligand_in_panel'] & custom_pairs['receptor_in_panel'])]
 
    print(f"\nCustom pairs usable with panel: {len(recoverable)}")
    print(f"Custom pairs missing from panel (dropped): {len(truly_missing)}")
 
    # Extend resource with recoverable custom pairs
    resource_extended = pd.concat(
        [resource, recoverable[['ligand', 'receptor']]],
        ignore_index=True
    ).drop_duplicates()
 
    print(f"Extended resource: {len(resource_extended)} pairs")
 
    # Collect all genes needed for custom pairs
    custom_genes_needed = set(recoverable['ligand']) | \
        {g for pair in recoverable['receptor'] for g in pair.split('_')}
else:
    resource_extended  = resource
    custom_genes_needed = set()

##################################################
# STEP 7 — Subset adata to SVGs + custom LR genes
##################################################
genes_to_keep = list(set(list(svgs)) | custom_genes_needed)
genes_to_keep = [g for g in genes_to_keep if g in adata.var_names]
 
adata_targeted = adata[:, genes_to_keep].copy()
print(f"\nTargeted adata shape: {adata_targeted.shape}")

# Recompute spatial neighbours on the gene subset
li.ut.spatial_neighbors(
    adata=adata_targeted,
    bandwidth=bandwidth,
    spatial_key="spatial"
)
print("Spatial neighbours recomputed on targeted adata.")

##################################################
# STEP 8 — Run inflow scores
##################################################
print("\nRunning LIANA inflow scores...")
lrdata = li.mt.inflow(
    adata_targeted,
    groupby='cell_type_hr',
    resource=resource_extended,
    use_raw=False
)
print("Inflow scores computed.")
print(lrdata.shape)
 
# Save global interactions CSV
lrdata.uns['global_interactions'].to_csv(
    os.path.join(out_dir, f"{sample_name}_lr_interactions.csv"),
    index=False)
print("Saved: LR interactions CSV")

##################################################
# STEP 9 — Moran's I on LR interactions
##################################################
svis = lrdata.uns['moranI'].index[
    (lrdata.uns['moranI']['pval_norm_fdr_bh'] <= 0.05) &
    (lrdata.uns['moranI']['I'] > 0.01)
]
print(f"\nSpatially variable LR interactions: {len(svis)}")
 
lrdata.uns['moranI'].loc[lrdata.var_names].sort_values(
    "I", ascending=False).to_csv(
    os.path.join(out_dir, f"{sample_name}_lr_moranI.csv"))
print("Saved: LR Moran's I table")