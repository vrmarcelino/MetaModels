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
parser.add_argument('-oa_ex', '--output_all_exports', help="""name of the file to store export exchanges file. Default = met_exchanges_all_exports_spp.csv""", required=False, default="met_exchanges_all_exports_spp.csv")
parser.add_argument('-oa_im', '--output_all_imports', help="""name of the file to store all import exchanges file. Default = met_exchanges_all_imports_spp.csv""", required=False, default="met_exchanges_all_imports_spp.csv")

parser.add_argument('-oc', '--output_core', help="""basename of the file to store core met exchanges. Default = met_exchanges_core_spp""", required=False, default="met_exchanges_core_spp")

args = parser.parse_args()
in_path = args.folder_w_exchange_files
core_def = int(args.core)
out_file_all_exports = args.output_all_exports
out_file_all_imports = args.output_all_imports
out_file_core = args.output_core

#in_path = "2_exchanges"
#core_def = 90
#out_file_core = "met_exchanges_core_spp"
#out_file_all_exports = "met_exchanges_all_exports_spp.csv"
#out_file_all_imports = "met_exchanges_all_imports_spp.csv"


## merge files
all_files = glob.glob(os.path.join(in_path, "exchanges_grow_*.csv"))
df_from_each_file = (pd.read_csv(f, sep=',') for f in all_files)
df_merged = pd.concat(df_from_each_file, ignore_index=True)

## remove media:
exch_df = df_merged.loc[df_merged['taxon'] != "medium"]

## make source_target column
exch_df['bin_met'] = exch_df['taxon'] + "_" + exch_df['metabolite']
print ("\nWarning is issued, but shouldn't be a problem here. Ignore.\n")


## separate imports and exports in different tables
exports_df = exch_df.loc[exch_df['direction'] == "export"]
imports_df = exch_df.loc[exch_df['direction'] == "import"]


## reshape table to have sample as rows and source_target (bin_met) as columns:
exports_df_samples_as_rows = exports_df.pivot(index = "sample_id", columns="bin_met", values="flux")
exports_df_samples_as_rows = exports_df_samples_as_rows.fillna(0) # replace NaNs with zeros

imports_df_samples_as_rows = imports_df.pivot(index = "sample_id", columns="bin_met", values="flux")
imports_df_samples_as_rows = imports_df_samples_as_rows.fillna(0) # replace NaNs with zeros


## save entire table to file
exports_df_samples_as_rows.to_csv(out_file_all_exports, index=True)
imports_df_samples_as_rows.to_csv(out_file_all_imports, index=True)


##################################
### stats on core interactions ###

def core_stats(in_df, core_def, out_file_core):
    n_samples = len(in_df)
    # number of individuals that must have the bin_met interaction for it to be considered core:
    threshold = core_def * n_samples / 100

    total_n_rxn = len(in_df.columns)
    print("Total number of species x metabolite exchanges (edges of the interactome): %i" %(total_n_rxn))

    # count non-zeros for each column
    count_non_zeros = in_df.astype(bool).sum(axis=0)
    wanted_edges = count_non_zeros > threshold

    core_X = in_df[wanted_edges.index[wanted_edges]]
    print ("Number of edges found in %i%% of individuals: %i" %(core_def, len(core_X.columns)))
    core_X.to_csv(out_file_core)

    ### strict core  (100% of individuals)
    core_100 = in_df.loc[:, (in_df != 0).all()]
    #core_100.to_csv("core100.csv")
    print ("Number of edges found in 100%% of individuals: %i" %(len(core_100.columns)))

# run it for exports and imports
out_file_core_exports = out_file_core + "_exports.csv"
out_file_core_imports = out_file_core + "_imports.csv"

print ("\nStats for export reactions:")
core_stats(exports_df_samples_as_rows,core_def,out_file_core_exports)

print ("\nStats for import reactions:")
core_stats(imports_df_samples_as_rows,core_def,out_file_core_imports)

print ("\nDone!\n")



