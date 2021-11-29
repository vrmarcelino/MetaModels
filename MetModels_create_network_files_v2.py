#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create the node and edge files for network analysis in R
Based on large dataset of exchange files (use the entire exchange folder as input)

using the healthy population only!

Created on 29/11/21
@author: V.R.Marcelino
"""
from argparse import ArgumentParser
import os, glob, csv
import pandas as pd


parser = ArgumentParser()
parser.add_argument('-f', '--folder_w_exchange_files', help="""Path to the folder containing exchange files produced by MICOM""", required=True)
parser.add_argument('-m', '--metadata', help="""Path to the metadata file - required to select healthy individuals""",required=True)
parser.add_argument('-s', '--spp_classification', help="""Path to the spp_classification.tsv file""",required=True)
parser.add_argument('-b', '--bigg_models', help="""Path to the bigg_models_w_classes_curated.tsv file""",required=True)
parser.add_argument('-o', '--output', help="""basename of output file. Default = healthy""", required=False, default="healthy")

args = parser.parse_args()
in_path = args.folder_w_exchange_files
metad_fp = args.metadata
binID2spp = args.spp_classification
in_bigg = args.bigg_models
out_nodes = args.output + "_nodes.csv"
out_edges = args.output + "_edges.csv"


#in_path = "2_exchanges"
#metad_fp = "metadata_rewiring_microbiome.csv"
#binID2spp = "wanted_spp_classification.tsv"
#in_bigg = "bigg_models_w_classes_curated.tsv"
#out_nodes = "nodes_biome3.csv"
#out_edges = "edges_biome3.csv"


## read metadata
metad_all_samples = pd.read_csv(metad_fp)
metad_all_samples.rename(columns={'Sample':'sample_id'}, inplace=True)

## read bins 2 spp classification map:
binID2taxa = {}
with open(binID2spp, mode='r') as inp:
    reader = csv.reader(inp, delimiter="\t")
    for row in reader:
        binID = row[0].replace(".","_")
        phylum = row[1].split(";")[1].replace("p__","")
        binID2taxa[binID] = phylum

## read bigg models (table to get the metabolite names)
metID2name = {}
with open (in_bigg, mode='r') as bg:
    reader = csv.reader(bg, delimiter="\t")
    for row in reader:
        metID = row[0]
        met_name = row[1]
        metID2name[metID] = met_name


## merge exchange files
all_files = glob.glob(os.path.join(in_path, "exchanges_grow_*.csv"))
df_from_each_file = (pd.read_csv(f, sep=',') for f in all_files)
df_merged = pd.concat(df_from_each_file, ignore_index=True)


### Clean exchange table
exch_df = df_merged.loc[df_merged['taxon'] != "medium"] # remove media

## remove the "cat" from sample names:
pd.options.mode.chained_assignment = None  # disables the warning, default='warn'
exch_df["sample_id"] = exch_df["sample_id"].str.replace("_cat", "")

## add 'HD' info and keep only healthy individuals:
exch_df_metad = pd.merge(exch_df, metad_all_samples[["sample_id","HD"]], on="sample_id", how="left")
exch_df_healthy = exch_df_metad[exch_df_metad.HD == "healthy"]

### additional table cleaning:
exchanges = exch_df_healthy.drop(columns=['sample_id', 'tolerance','flux','abundance','Unnamed: 0'])
exchanges = exchanges.drop_duplicates() #remove duplicate links (links that occur in >1 sample appearin multiple lines)
exchanges['reaction'] = exchanges['reaction'].str.replace('EX_','') # remove "Ex_"


### Order table to make a directional graph, add to new table
producers = exchanges[exchanges.direction == 'export']
producers = producers.rename(columns={"taxon": "source", "reaction": "target"})

consumers = exchanges[exchanges.direction == 'import']
consumers = consumers.rename(columns={"taxon": "target", "reaction": "source"})

edges = pd.concat([producers, consumers])


##################
## format nodes ##
##################

# get unique MAGs and Metabolites:
column_values = edges[["source", "target"]].values.ravel()
unique_values = pd.unique(column_values)

nodes = pd.DataFrame(unique_values, columns=["node"])
for index, row in nodes.iterrows():
    if "_e" in row['node']:
        nodes.at[index,'type'] = "metab"
        nodes.at[index,'type_numb'] = 1
        nodes.at[index,'name'] = metID2name[row['node']]
    else:
        nodes.at[index,'type'] = "mag"
        nodes.at[index,'type_numb'] = 0
        nodes.at[index, 'name'] = binID2taxa[row['node']]

nodes = nodes.sort_values(["type_numb","name"])

# save to file:
edges.to_csv(out_edges, index=False)
nodes.to_csv(out_nodes, index=False)


print ("""\n\nDONE!!
    Edges saved as %s
    Nodes saved as %s \n""" %(out_edges, out_nodes))
