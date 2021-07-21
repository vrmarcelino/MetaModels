#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Identify the impact of potential probiotics on modulating the production of
specific metabolites

Created on 1/7/21
@author: V.R.Marcelino
"""

import pandas as pd
import os
import re

import numpy as np
import matplotlib.pyplot as plt

# note - healthy people have more fumarate.

metab = "EX_acald_e"
elasticities_folder = "diseased"
bin2sp_map_fp = "HQ_bins_with_compl.csv"
output_fp = "elasticity_all_samples_" + metab + ".csv"

### make a dictionary of binID 2 lineage:
bin2sp_map = {}
with open (bin2sp_map_fp, 'r') as bins:
    next(bins)
    for line in bins:
        binID = line.split(",")[0]
        tax = line.split(",")[1]
        bin2sp_map[binID] = tax

### read elasticity files and merge results

# create new dataframe:
all_samples = pd.DataFrame()

for file in os.listdir(elasticities_folder):
    if file.endswith(".csv"):
        sample_name = re.split(r'_|.csv', file)[2]
        result_fp = os.path.join(elasticities_folder, file)

        df = pd.read_csv(result_fp, sep=',', index_col=0, names=["reaction","taxon","effector","direction","elasticity","type"])

        # keep only bacteria
        df = df[df['type'] == "abundance"]

        # keep only forward direction (check - done in Diener script,
        # probably because we want to assess the effects of increasing the nutrient/bacteria availability (not decrease it))
        df = df[df['direction'] == "forward"]

        # sum the elasticity of a particular MAG across members of the microbiome:
        grouped_df = df.groupby(["effector"]).elasticity.sum()

        grouped_df.rename(sample_name, inplace=True)
        all_samples = pd.concat([all_samples, grouped_df], sort=True, axis=1)



## replace NaNs with zeros:
all_samples = all_samples.fillna(0)

## add taxon:
all_samples['lineage'] = all_samples.index.map(bin2sp_map) # full lineage

## save to file:
all_samples.to_csv(output_fp)