#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to identify the core and accessory edges (i.e. metabolic exchanges) across samples

Note -- direction of interaction (import / export) is not taken into consideration to define core edges
but it is recorded (neg for import/pos for exp) for further analyses

Created on 26/8/21
@author: V.R.Marcelino
"""

import pandas as pd
from argparse import ArgumentParser
import os, glob


parser = ArgumentParser()

parser.add_argument('-f', '--folder_w_exchange_files', help="""Path to the folder containing exchange files produced by MICOM""", required=True)
parser.add_argument('-c', '--core', help="""Definition of core - min percentage of samples where the edge should be present to be considered core.
                    Default = 90""", required=False, default=90)
parser.add_argument('-oa', '--output_all', help="""name of the file to store all edges. Default = edges_all.csv""", required=False, default="edges_all.csv")
parser.add_argument('-oc', '--output_core', help="""name of the file to store core edges. Default = edges_core.csv""", required=False, default="edges_core.csv")

args = parser.parse_args()
in_path = args.folder_w_exchange_files
core_def = int(args.core)
out_file_all = args.output_all
out_file_core = args.output_core

#in_path = "2_exchanges"
#out_file = "something.csv"


## merge files
all_files = glob.glob(os.path.join(in_path, "exchanges_grow_*.csv"))
df_from_each_file = (pd.read_csv(f, sep=',') for f in all_files)
df_merged = pd.concat(df_from_each_file, ignore_index=True)

## remove media:
exch_df = df_merged.loc[df_merged['taxon'] != "medium"]

## make source_target column
exch_df['bin_met'] = exch_df['taxon'] + "_" + exch_df['metabolite']
print ("\nWarning is issued, but shouldn't be a problem here. Ignore.\n")

## reshape table to have sample as rows and source_target (bin_met) as columns:
exch_df_samples_as_rows = exch_df.pivot(index = "sample_id", columns="bin_met", values="flux")
exch_df_samples_as_rows = exch_df_samples_as_rows.fillna(0) # replace NaNs with zeros

## save entire table to file
exch_df_samples_as_rows.to_csv(out_file_all, index=True)

##################################
### stats on core interactions ###

###  loose core definition (x% of individuals)
n_samples = len(exch_df_samples_as_rows)
#number of individuals that must have the bin_met interaction for it to be considered core:
threshold = core_def * n_samples / 100

# count non-zeros for each column
count_non_zeros = exch_df_samples_as_rows.astype(bool).sum(axis=0)
wanted_edges = count_non_zeros > threshold

core_X = exch_df_samples_as_rows[wanted_edges.index[wanted_edges]]
print ("number of interactions found in %i%% of individuals: %i" %(core_def, len(core_X.columns)))
core_X.to_csv(out_file_core)

### strict core  (100% of individuals)
core_100 = exch_df_samples_as_rows.loc[:, (exch_df_samples_as_rows != 0).all()]
#core_100.to_csv("core100.csv")
print ("number of interactions found in 100%% of individuals: %i" %(len(core_100.columns)))

print ("Done!")