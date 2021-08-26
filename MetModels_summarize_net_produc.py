#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to merge the output of MICOM_grow_wf (exchanges_grow_xxxx.csv) and produce
a table with net production / consumption of metabolites by the microbiome
->> these are the exchanges with the media (“_m”), as they indicate the “excess" of
metabolites that are released or consumed from the environment.

The values are already corrected for species' relative abundance.

Created on 26/8/21
@author: V.R.Marcelino
"""

from argparse import ArgumentParser
import os, glob
import pandas as pd

parser = ArgumentParser()
parser.add_argument('-f', '--folder_w_exchange_files', help="""Path to the folder containing exchange files produced by MICOM""", required=True)
parser.add_argument('-o', '--output', help="""name of output file. Default = net_produc_consump_merged.csv""", required=False, default="net_produc_consump_merged.csv")

args = parser.parse_args()
in_path = args.folder_w_exchange_files
out_file = args.output

#in_path = "2_exchanges"
#out_file = "net_produc_consump_merged.csv"

## merge files
all_files = glob.glob(os.path.join(in_path, "exchanges_grow_*.csv"))
df_from_each_file = (pd.read_csv(f, sep=',') for f in all_files)
df_merged = pd.concat(df_from_each_file, ignore_index=True)

## keep only media:
media_rows = df_merged.loc[df_merged['taxon'] == "medium"]

# reshape table to have sample as rows and metab as columns:
net_df = media_rows.pivot(index = "sample_id", columns="metabolite", values="flux")
net_df = net_df.fillna(0) # replace NaNs with zeros

## save it
net_df.to_csv(out_file, index=True)

print ("Done!")

