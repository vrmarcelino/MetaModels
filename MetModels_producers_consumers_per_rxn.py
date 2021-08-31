#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to calculate the number of producers and consumers of each metabolite,
in order to estimate the metabolite's importance in the community and infer
community stability.

Created on 1/9/21
@author: V.R.Marcelino
"""

from argparse import ArgumentParser
import os, glob
import pandas as pd

parser = ArgumentParser()
parser.add_argument('-f', '--folder_w_exchange_files', help="""Path to the folder containing exchange files produced by MICOM""", required=True)
parser.add_argument('-o', '--output', help="""name of output file. Default = producers_consumers.csv""", required=False, default="producers_consumers.csv")

args = parser.parse_args()
in_path = args.folder_w_exchange_files
out_file = args.output

#in_path = "2_exchanges"
#out_file = "producers_consumers.csv"

## merge exchange files
all_files = glob.glob(os.path.join(in_path, "exchanges_grow_*.csv"))
df_from_each_file = (pd.read_csv(f, sep=',') for f in all_files)
df_merged = pd.concat(df_from_each_file, ignore_index=True)

## remove media:
exch_df = df_merged.loc[df_merged['taxon'] != "medium"]

## separate consumption and production fluxes into two columns:
pd.options.mode.chained_assignment = None  # disables the warning, default='warn'
mask = exch_df['flux'] < 0
exch_df['flux_production'] = exch_df['flux'].mask(mask)
exch_df['flux_consumption'] = exch_df['flux'].mask(~mask)
exch_df = exch_df.fillna(0) # replace NaNs with zeros


# calculate number of producers/consumers per metabolite, and the sum of their fluxes:
count_producers = exch_df.groupby(["sample_id", "metabolite"]).apply(lambda df: df.flux_production.astype(bool).sum(axis=0))
sum_product_flux = exch_df.groupby(["sample_id", "metabolite"]).apply(lambda df: sum(df.flux_production))

count_consumers = exch_df.groupby(["sample_id", "metabolite"]).apply(lambda df: df.flux_consumption.astype(bool).sum(axis=0))
sum_consump_flux = exch_df.groupby(["sample_id", "metabolite"]).apply(lambda df: sum(df.flux_consumption))

# merge the dataframes:
prod_con_summary = pd.concat([count_producers, sum_product_flux,count_consumers, sum_consump_flux], axis=1)
prod_con_summary.columns = ["n_producers", "flux_producers","n_consumers","flux_consumers"]

# remove metabolites that are not consumed by anyone:
prod_con_summary = prod_con_summary[prod_con_summary['n_consumers'] > 0]

## save it:
prod_con_summary.to_csv(out_file)

print ("Done!")