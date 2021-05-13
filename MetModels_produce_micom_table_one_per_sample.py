#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create the input tables for micom - one table per sample
Created on 13/05/21
@author: V.R.Marcelino
"""

import pandas as pd
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-kma', '--kma', help="""Path to the merged kma results at species level, produced with  0_PCAs_MAGs_species_level.R (2_ccm_otus_clean_aggregated.csv)""", required=True)
parser.add_argument('-GEMs', '--GEMs', help="""Path to the folder containing the Genome scale metabolic models""", required=True)
parser.add_argument('-b', '--sp2bin', help="file indicating the highest quality bin for each species, produced with select_MAGs.py (HQ_bins_with_compl.csv)", required=True)

parser.add_argument('-of', '--output_folder', default = '5_MICOM/0_MAGs_tables',
                    help='Path to the output folder where tables will be stored. Default = 5_MICOM/0_MAGs_tables', required=False)

args = parser.parse_args()
#input
kma_res_fp = args.kma
GEMs_fp = args.GEMs
sp2bin_fp=args.sp2bin
#output
community_types_fp = args.output_folder

#kma_res_fp = "2_ccm_otus_clean_aggregated_test.csv"
#GEMs_fp = "1_GEMs"
#sp2bin_fp="HQ_bins_with_compl.csv"
#community_types_fp = "5_MICOM/0_MAGs_tables"

if not os.path.exists(community_types_fp):
    os.makedirs(community_types_fp)


###### Store the preferred bin as another dict:
sp2bin ={}
with open(sp2bin_fp) as sp:
    next(sp) # skip first line
    for line in sp:
        split_line = line.split(',')
        binID = split_line[0]
        tax = split_line[1].replace(" ", "_")
        sp2bin[tax] = binID


###### read and parse KMA results - adding the HQ bins to the table.
df = pd.read_csv(kma_res_fp, index_col = 0, encoding='latin1')
df['binID'] = df.index.map(sp2bin)
#df.to_csv("CCM_table_with_bins.csv")
columns = df.columns.tolist()[:-1] # sample names


### create one table per sample:
# to stop warning:
pd.set_option('mode.chained_assignment', None)

for i in columns:
    new_df = df[['binID']].rename(columns={"binID": "id"})
    new_df['species'] = new_df.index
    new_df['sample_id'] = i
    new_df['file'] = GEMs_fp + "/" + df['binID'] + ".xml"
    new_df['abundance'] = df[i]

    # remove rows with zero counts:
    new_df = new_df.drop(new_df[new_df.abundance == 0].index)

    #save
    file_name =  community_types_fp + "/" + i + ".csv"
    new_df.to_csv(file_name, index=False)




print ("Done. Happy simulations.")