#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to merge the output of MICOM_coop_tradeoff (minimal_fluxes_exchange_xxx)
Needed when running MICOM_coop_tradeoff for individual samples separately

Created on 19/5/21
@author: V.R.Marcelino
"""

from argparse import ArgumentParser
import os, glob
import pandas as pd

parser = ArgumentParser()
parser.add_argument('-f', '--folder_w_exchange_files', help="""Path to the folder containing exchange files produced by MICOM""", required=True)
parser.add_argument('-o', '--output', help="""name of the merged minimal_fluxes_exchange table""", required=False, default="minimal_fluxes_exchange_merged.csv")

args = parser.parse_args()
in_path = args.folder_w_exchange_files
out_file = args.output

## merge files
all_files = glob.glob(os.path.join(in_path, "minimal_fluxes_exchange_*.csv"))
df_from_each_file = (pd.read_csv(f, sep=',') for f in all_files)
df_merged = pd.concat(df_from_each_file, ignore_index=True)

# re-order column names
df_merged = df_merged.reindex(sorted(df_merged.columns), axis=1)
df_merged = df_merged[['compartment'] + [ col for col in df_merged.columns if col != 'compartment' ]] # move column 'compartment' to the first position

# fill NaNs with zeors
df_merged = df_merged.fillna(0)

# remove ".pickle" from sample names:
df_merged["sample"] = df_merged["sample"].str.replace(".pickle", "")

df_merged.to_csv(out_file, index=False)


