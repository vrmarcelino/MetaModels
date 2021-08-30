#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to merge the output of MICOM_grow_wf (exchanges_grow_xxxx.csv) and produce
tables with Total production Total consumption of metabolites by the microbiome

--> corrected for species Absolute abundances

Created on 26/8/21
@author: V.R.Marcelino
"""

from argparse import ArgumentParser
import os, glob
import pandas as pd

parser = ArgumentParser()
parser.add_argument('-f', '--folder_w_exchange_files', help="""Path to the folder containing exchange files produced by MICOM""", required=True)
parser.add_argument('-kma', '--kma', help="""Path to the merged kma results at species level, produced with  0_PCAs_MAGs_species_level.R (2_ccm_otus_clean_aggregated.csv)""", required=True)
parser.add_argument('-b', '--sp2bin', help="file indicating the highest quality bin for each species, produced with select_MAGs.py (HQ_bins_with_compl.csv)", required=True)
parser.add_argument('-op', '--output_production', help="""name of file to save total production. Default = total_production.csv""", required=False, default="total_production.csv")
parser.add_argument('-oc', '--output_consumption', help="""name of file to save consumption. Default = total_consumption.csv""", required=False, default="total_consumption.csv")

args = parser.parse_args()
in_path = args.folder_w_exchange_files
kma_res_fp = args.kma
sp2bin_fp=args.sp2bin
out_file_prod = args.output_production
out_file_cons = args.output_consumption

#in_path = "2_exchanges"
#out_file_prod = "total_production.csv"
#kma_res_fp = "2_ccm_otus_clean_aggregated.csv"
#sp2bin_fp="HQ_bins_with_compl.csv"


###### Store the preferred bin as dict:
sp2bin ={}
with open(sp2bin_fp) as sp:
    next(sp) # skip first line
    for line in sp:
        split_line = line.split(',')
        binID = split_line[0]
        tax = split_line[1].replace(" ", "_")
        sp2bin[tax] = binID

###### read and parse KMA results (aggregated @ species level) - adding the representative bins to the table.
kma_df = pd.read_csv(kma_res_fp,index_col = 0, encoding='latin1')
kma_df['representative_bin'] = kma_df.index.map(sp2bin)
kma_df2 = kma_df.set_index('representative_bin')
kma_df2.rename_axis(None, inplace=True)

# For each sample/taxon - get their absolute abundance
kma_columns = kma_df2.stack().reset_index().rename(columns={'level_0':'rep_bin','level_1':'Sample', 0:'Abs_abundance'})


## merge exchange files
all_files = glob.glob(os.path.join(in_path, "exchanges_grow_*.csv"))
df_from_each_file = (pd.read_csv(f, sep=',') for f in all_files)
df_merged = pd.concat(df_from_each_file, ignore_index=True)

## remove media:
non_media_rows = df_merged.loc[df_merged['taxon'] != "medium"]

## Add absolute abundance:
exch_w_abs_abund = pd.merge(non_media_rows, kma_columns,  how='left', left_on=['taxon','sample_id'], right_on = ['rep_bin','Sample'])

## calculated product/consump corrected for absolute abundances:
exch_w_abs_abund['flux_weighted'] = exch_w_abs_abund['flux'] * exch_w_abs_abund['Abs_abundance']

#### Summarize production:
exports = exch_w_abs_abund.loc[exch_w_abs_abund['direction'] == "export"]

# reshape table to have sample as rows and metab as columns:
production_df = exports.pivot_table(index = "sample_id", columns="metabolite", values="flux_weighted", aggfunc='sum')
production_df = production_df.fillna(0) # replace NaNs with zeros

# save it
production_df.to_csv(out_file_prod, index=True)

#### Summarize consumption:
imports = exch_w_abs_abund.loc[exch_w_abs_abund['direction'] == "import"]

# reshape table to have sample as rows and metab as columns:
consumption_df = imports.pivot_table(index = "sample_id", columns="metabolite", values="flux_weighted", aggfunc='sum')
consumption_df = consumption_df.fillna(0) # replace NaNs with zeros

# save it
consumption_df.to_csv(out_file_cons, index=True)

print ("Done!")

