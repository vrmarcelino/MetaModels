#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create the input tables for micom - one table per sample
Created on 13/05/21
Updated on 18/09/21
@author: V.R.Marcelino
"""

import pandas as pd
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-kma', '--kma', help="""Path to the merged kma results at species level, produced with  MetModels_parse_KMA.py (1_merged_kma_res.csv)""", required=True)
parser.add_argument('-GEMs', '--GEMs', help="""Path to the folder containing the Genome scale metabolic models""", required=True)

parser.add_argument('-of', '--output_folder', default = '5_MICOM/0_MAGs_tables',
                    help='Path to the output folder where tables will be stored. Default = 5_MICOM/0_MAGs_tables', required=False)

args = parser.parse_args()
#input
kma_res_fp = args.kma
GEMs_fp = args.GEMs
#output
community_types_fp = args.output_folder

#kma_res_fp = "../2_Abundance_Calc/1_merged_kma_res.csv"
#GEMs_fp = "1_GEMs"
#community_types_fp = "5_MICOM/0_MAGs_tables"

if not os.path.exists(community_types_fp):
    os.makedirs(community_types_fp)



###### read and parse KMA results - adding the HQ bins to the table.
df = pd.read_csv(kma_res_fp, index_col = 0, encoding='latin1')
columns = df.columns.tolist()[:-1] # sample names
df.index.rename("binID", inplace=True)

### create one table per sample:
# to stop warning:
pd.set_option('mode.chained_assignment', None)


for i in columns:
    new_df = pd.DataFrame(index=df.index)
    new_df.index.rename("id", inplace=True)
    new_df['species'] = df[['Taxonomy']]
    new_df['sample_id'] = i
    new_df['file'] = GEMs_fp + "/" + new_df.index + ".xml"
    new_df['abundance'] = df[[i]]

    # remove rows with zero counts:
    new_df = new_df.drop(new_df[new_df.abundance == 0].index)

    #save
    file_name =  community_types_fp + "/" + i + ".csv"
    new_df.to_csv(file_name, index=True)


print ("Done. Happy simulations.")