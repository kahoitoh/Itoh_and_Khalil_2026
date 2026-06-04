#!/usr/bin/env python

# ==============================================================================
# Collect AP-1-high-associated topic peaks and DARs
# ------------------------------------------------------------------------------
# This script collects visually selected AP-1-high topic peak sets and
# differentially accessible regions from previously generated cisTopic outputs.
#
# Main steps:
#   1. Specify topics enriched in AP-1-high cells based on the topic heatmap.
#   2. Copy corresponding Otsu and top-3k topic BED files.
#   3. Copy AP-1-high vs AP-1-low DAR BED files.
#   4. Save selected region sets into a curated output directory.
#
# Notes:
#   - High-score topics should be selected manually by inspecting the heatmap.
#   - This script does not rerun cisTopic, topic binarization, or DAR detection.
# ==============================================================================

import os
import shutil
from glob import glob

# config -----------------------------------------------------------------------
OUT_DIR = "outs"
CELL_TYPE = "BA" # choose cell type of interest.

# Fill this after checking:
# outs/BA_topic_heatmap_Cell_cond.png
HIGH_SCORE_TOPICS = [
    # "Topic1",
    # "Topic7",
]

TOPIC_HEATMAP = os.path.join(
    OUT_DIR,
    f"outs/topic_heatmap_Cell_cond.png"
)

REGION_SET_DIR = os.path.join(OUT_DIR, "region_sets")

TOPIC_OTSU_DIR = os.path.join(REGION_SET_DIR, "Topics_otsu")
TOPIC_TOP3K_DIR = os.path.join(REGION_SET_DIR, "Topics_top_3k")
DAR_DIR = os.path.join(REGION_SET_DIR, "DARs_High_Low_each_celltype")

SELECTED_DIR = os.path.join(
    REGION_SET_DIR,
    f"{CELL_TYPE}_AP1_high_selected_region_sets"
)

SELECTED_OTSU_DIR = os.path.join(SELECTED_DIR, "Topics_otsu_High_score")
SELECTED_TOP3K_DIR = os.path.join(SELECTED_DIR, "Topics_top_3k_High_score")
SELECTED_DAR_DIR = os.path.join(SELECTED_DIR, "DARs_High_vs_Low")

os.makedirs(SELECTED_OTSU_DIR, exist_ok=True)
os.makedirs(SELECTED_TOP3K_DIR, exist_ok=True)
os.makedirs(SELECTED_DAR_DIR, exist_ok=True)


# helper -----------------------------------------------------------------------
def copy_if_exists(src, dst_dir):
    if os.path.exists(src):
        shutil.copy2(src, os.path.join(dst_dir, os.path.basename(src)))
        print(f"Copied: {src}")
    else:
        print(f"[warning] File not found: {src}")


# copy selected topic peaks ----------------------------------------------------
for topic in HIGH_SCORE_TOPICS:
    
    copy_if_exists(
        os.path.join(TOPIC_OTSU_DIR, f"{topic}.bed"),
        SELECTED_OTSU_DIR
    )
    
    copy_if_exists(
        os.path.join(TOPIC_TOP3K_DIR, f"{topic}.bed"),
        SELECTED_TOP3K_DIR
    )


# copy DAR peaks ---------------------------------------------------------------
dar_files = glob(os.path.join(DAR_DIR, f"*{CELL_TYPE}*.bed"))

if len(dar_files) == 0:
    print(f"[warning] No DAR BED files found for {CELL_TYPE} in {DAR_DIR}")

for dar_file in dar_files:
    copy_if_exists(dar_file, SELECTED_DAR_DIR)
