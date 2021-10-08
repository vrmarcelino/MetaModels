#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to identify the core and accessory metabolic exchanges across samples
(at the metabolite level)

outputs a summary of core vs. accessory metabolic exchanges across phenotypes

Created on 13 / Sep / 2021
Modified 22 / Sep / 2021

@author: V.R.Marcelino
"""

import pandas as pd
from argparse import ArgumentParser
import os, glob, re

parser = ArgumentParser()

parser.add_argument('-f', '--folder_w_exchange_files', help="""Path to the folder containing exchange files produced by MICOM""", required=True)
parser.add_argument('-c', '--core', help="""Definition of core - min percentage of samples where the edge should be present to be considered core.
                    Default = 95""", required=False, default=95)
parser.add_argument('-m', '--metadata', help="""Path to the metadata file""",required=True )
parser.add_argument('-g', '--grouping_header', help="""category to analyse core exchanges, must be a header in the metadata file (Diagnosis or HD)""",required=True)
parser.add_argument('-o', '--output_core', help="""file to store the summary of core met exchanges. Default = summary_met_exchanges_core.csv""", required=False, default="summary_met_exchanges_core.csv")

args = parser.parse_args()
in_path = args.folder_w_exchange_files
core_def = int(args.core)
metad_fp = args.metadata
grouping_header = args.grouping_header
out_file_core = args.output_core

#in_path = "2_exchanges"
#core_def = 95
#out_file_core = "summary_met_exchanges_core.csv"
#metad_fp = "metadata_rewiring_microbiome.csv"
#grouping_header = "HD"



def core_stats(in_df, core_def):
    """ Returns a list with:
    Total_[production|consumption], Core_[production|consumption], Proportion_core_[production|consumption]"""
    n_samples = len(in_df)
    # number of individuals that must have the bin_met interaction for it to be considered core:
    threshold = core_def * n_samples / 100

    total_n_rxn = len(in_df.columns)
    print("Total number of exchanged metabolites: %i" % (total_n_rxn))

    # count non-zeros for each column
    count_non_zeros = in_df.astype(bool).sum(axis=0)
    wanted_edges = count_non_zeros > threshold

    core_X = in_df[wanted_edges.index[wanted_edges]]
    number_core_rxn = len(core_X.columns)
    print("Number of metabolites being exchanged in %i%% of individuals: %i" % (core_def, number_core_rxn))

    # proportion:
    prop_core = round(number_core_rxn * 100 / total_n_rxn, 2)
    results = [total_n_rxn, number_core_rxn, prop_core]

    return (results)


## read metadata and identify groups
all_files = glob.glob(os.path.join(in_path, "exchanges_grow_*.csv"))

replace_string = in_path + "/|exchanges_grow_|_cat|.csv"
names_2_path_dict = {}
for i in (all_files):
    sample_id = re.sub(replace_string, '', i)
    names_2_path_dict[sample_id] = i

metad_all_samples = pd.read_csv(metad_fp)
metad_all_samples['file_path'] = metad_all_samples['Sample'].map(names_2_path_dict)
metad = metad_all_samples[metad_all_samples['file_path'].notna()] # keep only metadata samples for which we have samples

wanted_categories = list(set((metad[grouping_header])))

#process data in batches corresponding to phenotype (or healthy/disease)
summary_results = pd.DataFrame(["n_samples","Total_production","Core_production","Proportion_core_prod",
"Total_consumption","Core_consumption","Proportion_core_cons"], columns=['description'])


for cat in wanted_categories:
    print ("\n\nProcessing %s" %(cat))
    wanted_metad_df = metad[metad[grouping_header] == cat]
    wanted_files = list(wanted_metad_df['file_path'])

    ## merge files for each phenotype
    df_from_each_file = (pd.read_csv(f, sep=',') for f in wanted_files)
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


    ##################################
    ### stats on core interactions ###

    print ("\nproduction:")
    results_production = core_stats(exports_df_samples_as_rows,core_def)
    print ("\nconsumption:")
    results_consumption = core_stats(imports_df_samples_as_rows, core_def)

    wanted_info = [len(imports_df_samples_as_rows)] + results_production + results_consumption
    summary_results[cat] = wanted_info

# save it to file:
summary_results.to_csv(out_file_core, index=False)
print ("\n\nSummary saved to %s" %(out_file_core))
print ("\nDone!\n")

