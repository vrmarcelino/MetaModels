#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to copy GEMs to one folder per community type,
also creates a table indicating to which community type each bin belongs,
Created on 8/3/21
@author: V.R.Marcelino
"""

import pandas as pd
import shutil
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-kma', '--kma', help="""Path to the merged kam results, produced with MetModels_parse_KMA.py.""", required=True)
parser.add_argument('-dmm', '--dmm', help="""Path to the DMM sample assignments (1_sample_assignments_DMM_species.csv)""", required=True)
parser.add_argument('-GEMs', '--GEMs', help="""Path to teh folder containing the Genome scale metabolic models""", required=True)

parser.add_argument('-ot', '--output_table', default = 'bins2biomes.csv',
                    help='Path to the output table indicating to which community type each bin belongs. Default = bins2biomes.csv', required=False)
parser.add_argument('-of', '--output_folder', default = '2_DMM_GEMs',
                    help='Path to the output folder where GEMs will be stored. Default = 2_DMM_GEMs', required=False)

args = parser.parse_args()
#input
kma_res_fp = args.kma
commtypes_fp = args.dmm
GEMs_fp = args.GEMs
#output
out_matrix = args.output_table
community_types_fp = args.output_folder

#kma_res_fp = "1_merged_kma_res_test.csv"
#commtypes_fp = "1_sample_assignments_DMM_species.csv"
#GEMs_fp = "1_GEMs_test"
#out_matrix = "bins2biomes.csv"
#community_types_fp = "2_DMM_GEMs"

if not os.path.exists(community_types_fp):
    os.mkdir(community_types_fp)


###### Store community types as a dictionary:
comm_types ={}
with open(commtypes_fp) as ct:
    next(ct) # skip first line
    for line in ct:
        sample = line.split(",")[0].replace("_straindb.ccm.csv", "")
        sample = sample.replace('"', "") # remove the quotes
        DMM = line.split(",")[1].strip()
        comm_types[sample] = DMM


###### read and parse KMA results
df = pd.read_csv(kma_res_fp, index_col = 0, encoding='latin1')
columns = df.columns.tolist()[:-1] # sample names

for index, row in df.iterrows():
    for c in columns:
        if df.at[index,c] > 0:
            df.at[index,c] = comm_types[c]

pd.DataFrame.to_csv(df, out_matrix, index=True)

###### copy models to new folder

# create folders:
unique_comm_types = list(set(comm_types.values()))
for biome in unique_comm_types:
    new_dir = community_types_fp + "/" + biome
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)

# copy models:
for index, row in df.iterrows():
    for c in columns:
        assigned_bin = int(df.at[index, c])
        if assigned_bin == 0:
            pass
        else:
            bin_fp = community_types_fp + "/" + str(assigned_bin) + "/" + index + ".xml"
            GEM_fp = GEMs_fp + "/" + index + ".xml"
            shutil.copyfile(GEM_fp, bin_fp)

print ("Done!")

