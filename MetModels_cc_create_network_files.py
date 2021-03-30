#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create the node and edge files for network analysis in R
Created on 30/3/21
@author: V.R.Marcelino
"""

import pandas as pd


in_bigg = "bigg_models_w_classes.tsv"
in_micom = "result_exchanges_example.csv"

out_nodes = "nodes.csv"
out_edges = "edges.csv"

### Import tables
exchanges = pd.read_csv(in_micom, index_col=0)
bigg_models = pd.read_csv(in_bigg, sep='\t')

##################
## format edges ##
##################


### Clean exchange table
exchanges = exchanges[exchanges.taxon != 'medium'] # remove medium
exchanges = exchanges.drop(columns=['sample_id', 'tolerance'])
exchanges['reaction'] = exchanges['reaction'].str.replace('EX_','') # remove "Ex_"


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


