#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to select the high quality MAGs
Created on 16/3/21
@author: V.R.Marcelino
"""


import pandas as pd

gtdb_res="gtdb_bac_and_arch.tsv"
checkm_res="checkm.results.high_quality.tsv"


# make it a pandas df, then use groupby to get the more complete ones....
df = pd.read_csv(gtdb_res,sep='\t', index_col="binID")

# loop trough checkm results and store results in the df
with open(checkm_res) as c:
    for line in c:
        split_line = line.split()
        bin = split_line[0]
        completeness = split_line[12]
        df.at[bin, 'Completeness'] = float(completeness)


df.to_csv("all_bins_with_compl.csv")


# group by Classification (note - classification with captal "C" in older versions of gtdbtk)
grouped_df = df.sort_values('Completeness', ascending=False).drop_duplicates(['classification'])

# save to file:
grouped_df.to_csv("HQ_bins_with_compl.csv")
