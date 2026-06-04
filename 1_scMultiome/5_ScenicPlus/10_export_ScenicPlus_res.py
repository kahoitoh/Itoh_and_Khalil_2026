#!/usr/bin/env python

import mudata
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
from scenicplus.networks import create_nx_tables
from scenicplus.scenicplus_class import mudata_to_scenicplus

# import output
scplus_mdata = mudata.read("scplusmdata.h5mu")

scplus_mdata.uns["direct_e_regulon_metadata"]
scplus_mdata.uns["extended_e_regulon_metadata"]

# export to csv
scplus_mdata.uns["direct_e_regulon_metadata"].to_csv("df_direct_e_regulon_uns.csv", index=True)
scplus_mdata.uns["extended_e_regulon_metadata"].to_csv("df_extended_e_regulon_uns.csv", index=True)

scplus_mdata['direct_region_based_AUC'].to_df().to_csv("df_direct_region_based_AUC.csv", index=True)
scplus_mdata['extended_region_based_AUC'].to_df().to_csv("df_extended_region_based_AUC.csv", index=True)

# perform UMAP
eRegulon_region_AUC = anndata.concat(
    [scplus_mdata["direct_region_based_AUC"], scplus_mdata["extended_region_based_AUC"]],
    axis = 1,
)

eRegulon_region_AUC.obs = scplus_mdata.obs.loc[eRegulon_region_AUC.obs_names]
sc.pp.neighbors(eRegulon_region_AUC, use_rep = "X")
sc.tl.umap(eRegulon_region_AUC)

sc.pl.umap(eRegulon_region_AUC, color = "scATAC_counts:Condition", palette={"High": "#EE82EE", "Low": "#4169E1"})
fig = plt.gcf() 
fig.savefig("umap_eRegulonregion_cond.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# export TF to Gene table
scplus_obj = mudata_to_scenicplus(
    mdata = scplus_mdata,
    path_to_cistarget_h5 = "ctx_results.hdf5",
    path_to_dem_h5 = "dem_results.hdf5"
)

nx_tables_full = create_nx_tables(
    scplus_obj,                           
    eRegulon_metadata_key="eRegulon_metadata",
    subset_eRegulons=None, 
    subset_regions=None,
    subset_genes = None
)

nx_tables_full["Edge"]["TF2G"].to_csv("edges_TF2G_full.csv", index=False)