#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to identify the core and accessory metabolic exchanges across samples
(at the metabolite level)

Note -- direction of interaction (import / export) is not taken into consideration to define core edges
but it is recorded (neg for import/pos for exp) for further analyses

Created on 13 / Sep / 2021
@author: V.R.Marcelino
"""

import pandas as pd
from argparse import ArgumentParser
import os, glob


parser = ArgumentParser()

parser.add_argument('-f', '--folder_w_exchange_files', help="""Path to the folder containing exchange files produced by MICOM""", required=True)
parser.add_argument('-c', '--core', help="""Definition of core - min percentage of samples where the edge should be present to be considered core.
                    Default = 90""", required=False, default=90)
parser.add_argument('-oa_ex', '--output_all_exports', help="""name of the file to store export exchanges file. Default = met_exchanges_all_exports.csv""", required=False, default="met_exchanges_all_exports.csv")
parser.add_argument('-oa_im', '--output_all_imports', help="""name of the file to store all import exchanges file. Default = met_exchanges_all_imports.csv""", required=False, default="met_exchanges_all_imports.csv")

parser.add_argument('-oc', '--output_core', help="""basename of the file to store core met exchanges. Default = met_exchanges_core""", required=False, default="met_exchanges_core")

args = parser.parse_args()
in_path = args.folder_w_exchange_files
core_def = int(args.core)
out_file_all_exports = args.output_all_exports
out_file_all_imports = args.output_all_imports
out_file_core = args.output_core

#in_path = "2_exchanges"
#core_def = 90
#out_file_core = "met_exchanges_core"
#out_file_all_exports = "met_exchanges_all_exports.csv"
#out_file_all_imports = "met_exchanges_all_imports.csv"


## merge files
all_files = glob.glob(os.path.join(in_path, "exchanges_grow_*.csv"))
df_from_each_file = (pd.read_csv(f, sep=',') for f in all_files)
df_merged = pd.concat(df_from_each_file, ignore_index=True)

## remove media:
exch_df = df_merged.loc[df_merged['taxon'] != "medium"]

## separate imports and exports in different tables
exports_df = exch_df.loc[exch_df['direction'] == "export"]
imports_df = exch_df.loc[exch_df['direction'] == "import"]

## aggregate all fluxes for each sample
exports_df_agg = exports_df.groupby(["sample_id","metabolite"]).agg({"flux": "sum"}).reset_index()
imports_df_agg = imports_df.groupby(["sample_id","metabolite"]).agg({"flux": "sum"}).reset_index()


## reshape table to have sample as rows and metabolites as columns:
exports_df_samples_as_rows = exports_df_agg.pivot(index = "sample_id", columns="metabolite", values="flux")
exports_df_samples_as_rows = exports_df_samples_as_rows.fillna(0) # replace NaNs with zeros

imports_df_samples_as_rows = imports_df_agg.pivot(index = "sample_id", columns="metabolite", values="flux")
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
    print("Total number of exchanged metabolites: %i" %(total_n_rxn))

    # count non-zeros for each column
    count_non_zeros = in_df.astype(bool).sum(axis=0)
    wanted_edges = count_non_zeros > threshold

    core_X = in_df[wanted_edges.index[wanted_edges]]
    print ("Number of metabolites being exchanged in %i%% of individuals: %i" %(core_def, len(core_X.columns)))
    core_X.to_csv(out_file_core)

    ### strict core  (100% of individuals)
    core_100 = in_df.loc[:, (in_df != 0).all()]
    #core_100.to_csv("core100.csv")
    print ("Number of metabolites being exchanged in 100%% of individuals: %i" %(len(core_100.columns)))

# run it for exports and imports
out_file_core_exports = out_file_core + "_exports.csv"
out_file_core_imports = out_file_core + "_imports.csv"

print ("\nStats for export reactions:")
core_stats(exports_df_samples_as_rows,core_def,out_file_core_exports)

print ("\nStats for import reactions:")
core_stats(imports_df_samples_as_rows,core_def,out_file_core_imports)

print ("\nDone!\n")

