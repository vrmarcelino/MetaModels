#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Identify the impact of B. longum on modulating the production of
all metabolites

Created on 21/7/21
@author: V.R.Marcelino
"""

import pandas as pd
import os
import re


elasticities_folder = "0_MICOM_elast_filt"
bigg_models_fp="bigg_models_w_classes_curated.tsv"
output_fp = "elasticity_all_samples.csv"


### make a dictionary of metabolite_ID 2 metabolite name:
met2name_map = {}
with open (bigg_models_fp, 'r') as mets:
    next(mets)
    for line in mets:
        biggID_exchanges = "EX_"+ line.split("\t")[0]
        biggID_medium = biggID_exchanges + "_m"
        name = line.split("\t")[1]
        met2name_map[biggID_exchanges] = name
        met2name_map[biggID_medium] = name


### read elasticity files and merge results

# create new dataframe:
all_samples = pd.DataFrame()

for file in os.listdir(elasticities_folder):
    if file.endswith(".csv"):
        sample_name = re.split(r'_', file)[1]
        result_fp = os.path.join(elasticities_folder, file)

        df = pd.read_csv(result_fp, sep=',', index_col=0, names=["reaction","taxon","effector","direction","elasticity","type"])

        # sum the elasticity of a particular reaction across taxa affected by B. longum:
        grouped_df = df.groupby(["reaction"]).elasticity.sum()

        grouped_df.rename(sample_name, inplace=True)
        all_samples = pd.concat([all_samples, grouped_df], sort=True, axis=1)



## replace NaNs with zeros:
all_samples = all_samples.fillna(0)

## add taxon:
all_samples['met_name'] = all_samples.index.map(met2name_map) # full lineage

## save to file:
all_samples.to_csv(output_fp)


