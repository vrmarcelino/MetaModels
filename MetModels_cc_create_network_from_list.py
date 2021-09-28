#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create the node and edge files for network analysis in R
Takes as input a folder with  micom exchange results and a list of prefixes

Created on 27/09/21
@author: V.R.Marcelino
"""

import pandas as pd
from argparse import ArgumentParser
import os, glob

parser = ArgumentParser()

parser.add_argument('-if', '--in_folder', help="""Path to the micom exchange results (e.g. 2_exchanges)""", required=True)
parser.add_argument('-s', '--samples', help="""text file indicating the prefixes of the samples to be analysed """, required=True)
parser.add_argument('-m', '--metabolite', help="""generate nodes/edges for this specific metabolite (e.g. h2s_e). Default = all""", required=False, default="all")
parser.add_argument('-sp', '--spp_classifications', help="""Path to tab-sep file containing binID in one column and spp classification in teh other(e.g. wanted_spp_classification.tsv)""", required=True)

parser.add_argument('-on', '--out_nodes', help="output file name for nodes", required=True)
parser.add_argument('-oe', '--out_edges', help="output file name for edges", required=True)


args = parser.parse_args()
in_folder = args.in_folder
in_samples = args.samples
out_nodes = args.out_nodes
out_edges = args.out_edges
metabolite = args.metabolite
sp_class = args.spp_classifications

#in_folder = "2_exchanges"
#in_samples = "wanted_prefix.txt"
#out_nodes = "nodes.csv"
#out_edges = "edges.csv"
#metabolite = "h2s_e"
#sp_class = "wanted_spp_classification.tsv"


####################################
### store bins2spp file into dict ###
####################################

bins2spp_dict = {}
with open(sp_class) as f:
  for line in f:
    (binID, lineage)=line.split('\t')
    binID = binID.replace('.',"_")
    binID_produc = binID + "_p"
    binID_consum = binID + "_c"
    bins2spp_dict[binID_produc]=lineage
    bins2spp_dict[binID_consum] = lineage


####################################
### merge wanted exchange tables ###
####################################

# wanted_files list
with open(in_samples) as f:
    samples_prefix = f.read().splitlines()

wanted_samples_fp = []
for s in samples_prefix:
    sample_wildcard = "exchanges_grow_" + s + "*csv"
    fp_card = os.path.join(in_folder, sample_wildcard)
    fp = glob.glob(fp_card)
    #check if sample was found:
    if len(fp) > 0:
        wanted_samples_fp.append(fp[0])
    else:
        print ("\n WARNING: sample %s not found. Skipping... \n" %(s))

## merge exchange files
df_from_each_file = (pd.read_csv(f, sep=',') for f in wanted_samples_fp)
df_merged = pd.concat(df_from_each_file, ignore_index=True)

## keep only wanted metabolite - if given:
if metabolite != "all":
    df_merged = df_merged[df_merged.metabolite == metabolite]

### Clean exchange table
exchanges_all = df_merged[df_merged.taxon != 'medium'] # remove medium
exchanges_all = exchanges_all.drop(columns=['tolerance'])


### Calculate weighted flux (flux * rel_abundance) - abs converts negative to positive values
exchanges_all['flux_weighted'] = abs(exchanges_all['flux']) *  exchanges_all['abundance']



##################
## format edges ##
##################

### group dataframe by taxon, metabolite and direction, calculating average/sum...
# note that abundance may vary for producers and consumers (depends on their mean abundance within producers/consumers)
exchanges = exchanges_all.groupby(["taxon","metabolite","direction"]).agg({'abundance':['mean'],'flux':['mean', 'count'], 'flux_weighted':['sum']})
exchanges.columns = ['_'.join(col) for col in exchanges.columns.values]
exchanges.reset_index(level=('taxon','metabolite','direction'), inplace=True) # convert direction/taxon and metabolite into columns
exchanges = exchanges.rename(columns={"flux_count": "occurrences"})
exchanges = exchanges.rename(columns={"abundance_mean": "rel_abundance_mean"})

#exchanges_grouped.to_csv("check.csv")

### Order table to make a directional graph, add to new table
# silence SettingWithCopyWarning warning:
pd.options.mode.chained_assignment = None  # default='warn'

producers = exchanges[exchanges.direction == 'export']
producers['taxon'] = producers['taxon'] + "_p"
producers = producers.rename(columns={"taxon": "source", "metabolite": "target"})

consumers = exchanges[exchanges.direction == 'import']
consumers['taxon'] = consumers['taxon'] + "_c"
consumers = consumers.rename(columns={"taxon": "target", "metabolite": "source"})

edges = pd.concat([producers, consumers])


##################
## format nodes ##
##################

# get unique MAGs and Metabolites
column_values = edges[["source", "target"]].values.ravel()
unique_values = pd.unique(column_values)

###### make dictionaries with mean abundances and etc (to determine node sizes / labels)
producers_abund = pd.Series(producers.rel_abundance_mean.values,index=producers.source).to_dict()
consumers_abund = pd.Series(consumers.rel_abundance_mean.values,index=consumers.target).to_dict()

producers_flux_mean = pd.Series(producers.flux_mean.values,index=producers.source).to_dict()
consumers_flux_mean = pd.Series(consumers.flux_mean.values,index=consumers.target).to_dict()

producers_occurrences = pd.Series(producers.occurrences.values,index=producers.source).to_dict()
consumers_occurrences = pd.Series(consumers.occurrences.values,index=consumers.target).to_dict()

producers_flux_weighted_sum = pd.Series(producers.flux_weighted_sum.values,index=producers.source).to_dict()
consumers_flux_weighted_sum = pd.Series(consumers.flux_weighted_sum.values,index=consumers.target).to_dict()




# identify mags and metabolites, producers and consumers
nodes = pd.DataFrame(unique_values, columns=["node"])
for index, row in nodes.iterrows():
    if "_e" in row['node']:
        nodes.at[index,'type'] = "metab"
        nodes.at[index,'type_numb'] = 1
        nodes.at[index, 'lineage'] = 'metab'
        nodes.at[index, 'prod_cons'] = 'metab'
        nodes.at[index, 'rel_abundance_mean'] = 10 # random number!! - make it the largest ball!
        nodes.at[index, 'flux_mean'] = 10 # random number!! - make it the largest ball!
        nodes.at[index, 'occurrences'] = 10 # random number!! - make it the largest ball!
        nodes.at[index, 'flux_weighted_sum'] = 10 # random number!! - make it the largest ball!


    else:
        nodes.at[index,'type'] = "mag"
        nodes.at[index,'type_numb'] = 0

        bin_id = nodes.at[index, 'node']
        nodes.at[index, 'lineage'] = bins2spp_dict[bin_id]

        if "_p" in row['node']:
            nodes.at[index, 'prod_cons'] = 'p'
            nodes.at[index,'rel_abundance_mean'] = producers_abund[bin_id]
            nodes.at[index, 'flux_mean'] = producers_flux_mean[bin_id]
            nodes.at[index, 'occurrences'] = producers_occurrences[bin_id]
            nodes.at[index, 'flux_weighted_sum'] = producers_flux_weighted_sum[bin_id]
        elif "_c" in row['node']:
            nodes.at[index, 'prod_cons'] = 'c'
            nodes.at[index, 'rel_abundance_mean'] = consumers_abund[bin_id]
            nodes.at[index, 'flux_mean'] = consumers_flux_mean[bin_id]
            nodes.at[index, 'occurrences'] = consumers_occurrences[bin_id]
            nodes.at[index, 'flux_weighted_sum'] = consumers_flux_weighted_sum[bin_id]
        else:
            print ("\n\nCan't tell if this is a consumer or producer, check!\n\n")


## consider adding an extra column to flag species that are both producers and consumers

# save to file:
edges.to_csv(out_edges, index=False)
nodes.to_csv(out_nodes, index=False)

print ("\nDone!\n")

