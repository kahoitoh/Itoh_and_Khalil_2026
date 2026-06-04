#!/usr/bin/env python

# ==============================================================================
# cisTopic analysis and DAR detection for AP-1-high and AP-1-low cells
# ------------------------------------------------------------------------------
# This script creates a cisTopic object from 10x Multiome ATAC fragments,
# selects an LDA topic model, binarizes topics, imputes accessibility scores,
# and identifies differentially accessible regions between AP-1-high and
# AP-1-low cells.
#
# Main steps:
#   1. Load cell annotation metadata and filtered barcodes.
#   2. Create a cisTopic object using consensus peak regions.
#   3. Run MALLET LDA models and select the final topic model.
#   4. Visualize topic-cell associations.
#   5. Binarize topics using ntop and Otsu methods.
#   6. Impute accessibility and detect DARs between AP-1-high and AP-1-low cells.
#   7. Export topic regions and DARs as BED files.
#
# Notes:
#   - Consensus peak regions are used as the feature set.
#   - The selected model should be adjusted based on model evaluation results.
#   - Cell-type-specific contrasts should be updated for each analysis.
# ==============================================================================

import os
import pickle
import numpy as np
import pandas as pd
import polars as pl
import matplotlib.pyplot as plt

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
from pycisTopic.lda_models import run_cgs_models_mallet, evaluate_models
from pycisTopic.clust_vis import cell_topic_heatmap
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
from pycisTopic.utils import region_names_to_coordinates

# config
out_dir = "outs"
path_to_regions = os.path.join(
    out_dir,
    "consensus_peak_calling",
    "consensus_regions.bed"
)
path_to_blacklist = "mm10.blacklist.bed"
pycistopic_qc_output_dir = "outs/qc"
fragments_dict = {
    "10x_multiome_brain": "path_to_your_CellRanger_Aggr_dir/outs/atac_fragments.tsv.gz"
}
mallet_path="Mallet-202108/bin/mallet"

REGION = "BLA"
CELL_TYPE = "BA" # choose cell type of interest.
SELECT_MODEL = 15 # 10~15 depending on cell numbers. 
N_CPU_MODEL = 40
N_CPU_DAR = 5

# import data
cell_data = pd.read_table(f"path/to/TRAPed_{REGION}_{CELL_TYPE}_anno.tsv", index_col = 0) # TRAPed annotation file of cell type of interest
    
sample_id_to_barcodes_passing_filters = {
    sid: grp['barcode'].values
    for sid, grp in cell_data.groupby('Sample')
}

# create cisTopic object
cistopic_obj_list = []
for sample_id in fragments_dict:
    sample_metrics = pl.read_parquet(
        os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
    ).to_pandas().set_index("CB").loc[ sample_id_to_barcodes_passing_filters[sample_id] ]
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,
        path_to_blacklist = path_to_blacklist,
        metrics = sample_metrics,
        valid_bc = sample_id_to_barcodes_passing_filters[sample_id],
        n_cpu = 1,
        project = sample_id,
        split_pattern = '-'
    )
    cistopic_obj_list.append(cistopic_obj)

cistopic_obj = cistopic_obj_list[0]
print(cistopic_obj)

# add metadata
cistopic_obj.add_cell_data(cell_data, split_pattern='-')

cistopic_obj.cell_data["Condition"].value_counts(dropna=False)
cistopic_obj.cell_data["Cell"].value_counts(dropna=False)

# Run models to decide n_topics
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    n_cpu=N_CPU_MODEL,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path="mallet",
    save_path="mallet",
    mallet_path=mallet_path,
)

model = evaluate_models(
    models,
    select_model = SELECT_MODEL, 
    return_model = True
)

cistopic_obj.add_LDA_model(model)

fig = plt.gcf()
fig.savefig("outs/model.png", dpi=300, bbox_inches="tight") 
plt.close(fig)

# visual check of topics
cistopic_obj.cell_data["Cell_cond"] = (
    cistopic_obj.cell_data["Cell"].astype(str)
    + "-"
    + cistopic_obj.cell_data["High_low"].astype(str)
)

cell_topic_heatmap(
    cistopic_obj,
    variables = ['Cell_cond'],
    scale = False,
    legend_loc_x = 1.0,
    legend_loc_y = -1.2,
    legend_dist_y = -1,
    figsize = (15, 15)
)

fig = plt.gcf()   
fig.savefig("outs/topic_heatmap_Cell_cond.png", dpi=300, bbox_inches="tight")
plt.close(fig)


# Topic binarization
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000,
    plot=True, num_columns=5
)

fig = plt.gcf()   
fig.savefig("outs/topic_binarization_ntop.png", dpi=300, bbox_inches="tight")
plt.close(fig)

region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=5
)

binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100)

# impute count table
imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)

normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True
)

# compute differentially accessible peaks using imputed peaks
contrast_high = f"{CELL_TYPE}-High"
contrast_low = f"{CELL_TYPE}-Low"

markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='Cell_cond',
    var_features=variable_regions,
    contrasts = [[[contrast_high], [contrast_low]]],  
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=N_CPU_DAR,
    _temp_dir='/tmp',
    split_pattern = '-',
    ignore_reinit_error=True
)

# export peaks
os.makedirs(os.path.join(out_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "DARs_High_Low_each_celltype"), exist_ok = True)

for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for cell_type in markers_dict:
    region_names_to_coordinates(
        markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "DARs_High_Low_each_celltype", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )

pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)