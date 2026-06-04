#!/usr/bin/env python

# modify genome annotation chromosome format -----------------------------------
# Workaround for chromosome naming compatibility across genome annotation,
# chromsizes, and peak files. SCENIC+ downstream steps expect chromosome names
# to use the "chr" prefix.

import pandas as pd

ga = pd.read_csv("genome_annotation.tsv", sep="\t")

ga["Chromosome"] = ga["Chromosome"].astype(str)

def fix_chr(x: str) -> str:
    if x.startswith("chr"):
        return x
    if x in ["MT", "Mt", "M"]:
        return "chrM"
    return "chr" + x

ga["Chromosome"] = ga["Chromosome"].map(fix_chr)

ga.to_csv("genome_annotation.tsv", sep="\t", index=False)