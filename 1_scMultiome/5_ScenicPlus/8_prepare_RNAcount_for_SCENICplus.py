#!/usr/bin/env python

import os
import scanpy as sc
import pandas as pd
import scipy.io

# config -----------------------------------------------------------------------
OUT_DIR = "outs"
REGION = "BLA"
CELL_TYPE = "BA" # choose cell type of interest.

# Prepare RNA adata
X = scipy.io.mmread(f"path/to/TRAPed_{REGION}_{CELL_TYPE}_rna_counts.mtx")
X = X.T.tocsr()

genes = pd.read_csv(f"path/to/TRAPed_{REGION}_{CELL_TYPE}_genes.tsv", header=None, sep="\t")[0].values

cells = pd.read_csv(f"path/to/TRAPed_{REGION}_{CELL_TYPE}_barcodes.tsv", header=None, sep="\t")[0].values

adata = sc.AnnData(X=X)
adata.var_names = genes
adata.obs_names = cells

meta = pd.read_csv(f"path/to/TRAPed_{CELL_TYPE}_metadata.tsv", sep="\t")
meta = meta.set_index("barcode")
adata.obs = meta.loc[adata.obs_names]

adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.obs_names = [
    bc if bc.endswith("-10x_multiome_brain") else f"{bc}-10x_multiome_brain"
    for bc in adata.obs_names
]

adata.write_h5ad("rna_raw_for_SCENICplus.h5ad")