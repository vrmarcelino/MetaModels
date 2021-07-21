#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Given a specific metabolite, calculate the contribution of different MAGs/species
to its production and consumption across samples

Created on 29/6/21
@author: V.R.Marcelino
"""

import pandas as pd
import os
import re

import numpy as np
import matplotlib.pyplot as plt

# note - healthy people have more fumarate.

metab = "EX_acald_e"
exchanges_folder = "0_MICOM_tradeoffs/healthy_test"
bin2sp_map_fp = "HQ_bins_with_compl.csv"
output_fp = "all_samples_" + metab + ".csv"


### make a dictionary of binID 2 lineage:
bin2sp_map = {}
with open (bin2sp_map_fp, 'r') as bins:
    next(bins)
    for line in bins:
        binID = line.split(",")[0]
        tax = line.split(",")[1]
        bin2sp_map[binID] = tax


### read exchange files and merge results

# create new dataframe:
all_samples = pd.DataFrame()

for file in os.listdir(exchanges_folder):
    if file.endswith(".csv"):
        sample_name = re.split(r'_|.csv', file)[3]
        result_fp = os.path.join(exchanges_folder, file)

        df = pd.read_csv(result_fp, sep=',', index_col=0)
        wanted_metab = df[metab]

        wanted_metab.rename(sample_name, inplace=True)
        all_samples = pd.concat([all_samples, wanted_metab], sort=True, axis=1)


# remove species that were not present in all individuals?
# keep only the top contributors/consumers?!
#all_samples.to_csv("check.csv")

## remove rows where all values are zero or NaN:
all_samples = all_samples.fillna(0) # convert NaNs to zeros:
all_samples = all_samples.loc[(all_samples!=0).any(axis=1)] # removes zeros

## sort values
all_samples["net_prod_cons"] = all_samples.sum(axis=1, skipna=True)
all_samples = all_samples.sort_values('net_prod_cons')

## add taxon:
all_samples['lineage'] = all_samples.index.map(bin2sp_map) # full lineage
#all_samples[['Domain','Phylum','Class','Order','Family','Genus','Species']] = all_samples['lineage'].str.split(';', expand=True)


## save to file:
all_samples.to_csv(output_fp)


