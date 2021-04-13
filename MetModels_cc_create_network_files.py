#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create the node and edge files for network analysis in R
Created on 30/3/21
@author: V.R.Marcelino
"""

import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('-bigg', '--bigg_models', help="""Path to the file bigg_models_w_classes_curated.tsv""", required=True)
parser.add_argument('-i', '--in_micom', help="""Path to the micom exchnage results (e.g. exchanges_grow_sample_1.csv)""", required=True)
parser.add_argument('-r', '--in_reactions', help="""Path to tab-sep file containing the total number of reactions per GEM (e.g. all_reactions.txt)""", required=True)

parser.add_argument('-on', '--out_nodes', help="output file name for nodes", required=True)
parser.add_argument('-oe', '--out_edges', help="output file name for edges", required=True)


args = parser.parse_args()
in_bigg = args.bigg_models
in_micom = args.in_micom
in_n_react = args.reactions
out_nodes = args.out_nodes
out_edges = args.out_edges


#in_bigg = "bigg_models_w_classes_curated.tsv"
#in_micom = "exchanges_grow_sample_3.csv"
#in_n_react = "all_reactions.txt"
#out_nodes = "nodes_biome3.csv"
#out_edges = "edges_biome3.csv"

### Import tables
exchanges = pd.read_csv(in_micom, index_col=0)
bigg_models = pd.read_csv(in_bigg, sep='\t')
reactions = pd.read_csv(in_n_react, sep='\t', names=["taxon","n_reactions"])

##################
## format edges ##
##################


### Clean exchange table
exchanges = exchanges[exchanges.taxon != 'medium'] # remove medium
exchanges = exchanges.drop(columns=['sample_id', 'tolerance'])
exchanges['reaction'] = exchanges['reaction'].str.replace('EX_','') # remove "Ex_"

### Add number of reactions:
exchanges = exchanges.merge(reactions, how = 'inner', on = ['taxon'])


### Order table to make a directional graph, add to new table
producers = exchanges[exchanges.direction == 'export']
producers = producers.rename(columns={"taxon": "source", "reaction": "target"})

consumers = exchanges[exchanges.direction == 'import']
consumers = consumers.rename(columns={"taxon": "target", "reaction": "source"})

edges = pd.concat([producers, consumers])


### add classes info to the table:
bigg_models = bigg_models.rename(columns={"bigg_id":"metabolite", "SuperClass":"superclass", 
                                          "Class":"class","SubClass":"subclass"})

edges = pd.merge(edges, bigg_models, on="metabolite")

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
    else:
        nodes.at[index,'type'] = "mag"
        nodes.at[index,'type_numb'] = 0



# save to file:
edges.to_csv(out_edges, index=False)
nodes.to_csv(out_nodes, index=False)

