#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to copy GEMs - one per species! - to one folder per sample,
retaining only the TOP 30 most abundant bins in each sample.
also creates a table indicating to which community type each bin belongs,
Created on 23/3/21
@author: V.R.Marcelino
"""

import pandas as pd
import shutil
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-kma', '--kma', help="""Path to the merged kam results at species level, produced with  0_PCAs_MAGs_species_level.R.""", required=True)
parser.add_argument('-dmm', '--dmm', help="""Path to the DMM sample assignments (1_sample_assignments_DMM_species.csv)""", required=True)
parser.add_argument('-GEMs', '--GEMs', help="""Path to teh folder containing the Genome scale metabolic models""", required=True)
parser.add_argument('-b', '--sp2bin', help="file indicating the highest quality bin for each species, produced with select_MAGs.py", required=True)

parser.add_argument('-ot', '--output_table', default = 'bins2biomes.csv',
                    help='Path to the output table indicating to which community type each bin belongs. Default = bins2biomes.csv', required=False)
parser.add_argument('-of', '--output_folder', default = '2_DMM_GEMs',
                    help='Path to the output folder where GEMs will be stored. Default = 3_DMM_GEMs', required=False)

args = parser.parse_args()
#input
kma_res_fp = args.kma
commtypes_fp = args.dmm
GEMs_fp = args.GEMs
sp2bin_fp=args.sp2bin
#output
out_matrix = args.output_table
community_types_fp = args.output_folder

#kma_res_fp = "2_ccm_otus_clean_aggregated.csv"
#commtypes_fp = "1_sample_assignments_DMM_species.csv"
#GEMs_fp = "1_GEMs_test"
#sp2bin_fp="HQ_bins_with_compl.csv"
#out_matrix = "bins2biomes.csv"
#community_types_fp = "3_DMM_GEMs"

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

###### Store the preferred bin as another dict:
sp2bin ={}
with open(sp2bin_fp) as sp:
    next(sp) # skip first line
    for line in sp:
        split_line = line.split(',')
        binID = split_line[0]
        tax = split_line[1].replace(" ", "_")
        sp2bin[tax] = binID


###### read and parse KMA results, keeping only the 30 most abundant bins!
df = pd.read_csv(kma_res_fp, index_col = 0, encoding='latin1')
df['binID'] = df.index.map(sp2bin)

#identify the 30 most abundat bins in each sample, store as a dict:
columns = df.columns.tolist()[:-1] # sample names
top30 = {}
for c in columns:
    top30_df = df.nlargest(30, c)
    selected_bins = top30_df['binID'].tolist()
    top30[c] = selected_bins


# flag the bins that made the top 30
for index, row in df.iterrows():
    for c in columns:
        if df.at[index,c] > 0:
            if df.at[index,'binID'] in top30[c]:
                df.at[index,c] = comm_types[c]
            else:
                df.at[index, c] = 0 # remove bin from sample by giving it a zero abundance

pd.DataFrame.to_csv(df, out_matrix, index=True)

###### copy models to new folder

# create folders:
for sample in comm_types.keys():
    new_dir = community_types_fp + "/" + comm_types[sample] + "/" + sample
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

# copy models:
for index, row in df.iterrows():
    for c in columns:
        assigned_bin = int(df.at[index, c])
        if assigned_bin == 0:
            pass
        else:
            bin_ID = df.at[index, 'binID']
            sample_fp = community_types_fp + "/" + str(assigned_bin) + "/" + c + "/" + bin_ID + ".xml"

            GEM_fp = GEMs_fp + "/" + bin_ID + ".xml"
            print (sample_fp)
            print (GEM_fp)
            shutil.copyfile(GEM_fp, sample_fp)

print ("Done!")

