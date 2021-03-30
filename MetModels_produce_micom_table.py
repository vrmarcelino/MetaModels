#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create the input tables for micom - one per community type.
Created on 30/3/21
@author: V.R.Marcelino
"""

import pandas as pd
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-kma', '--kma', help="""Path to the merged kam results at species level, produced with  0_PCAs_MAGs_species_level.R.""", required=True)
parser.add_argument('-dmm', '--dmm', help="""Path to the DMM sample assignments (1_sample_assignments_DMM_species.csv)""", required=True)
parser.add_argument('-GEMs', '--GEMs', help="""Path to the folder containing the Genome scale metabolic models""", required=True)
parser.add_argument('-b', '--sp2bin', help="file indicating the highest quality bin for each species, produced with select_MAGs.py", required=True)

parser.add_argument('-of', '--output_folder', default = '5_MICOM/0_MAGs_tables',
                    help='Path to the output folder where tables will be stored. Default = 5_MICOM/0_MAGs_tables', required=False)

args = parser.parse_args()
#input
kma_res_fp = args.kma
commtypes_fp = args.dmm
GEMs_fp = args.GEMs
sp2bin_fp=args.sp2bin
#output
community_types_fp = args.output_folder

#kma_res_fp = "2_ccm_otus_clean_aggregated.csv"
#commtypes_fp = "1_sample_assignments_DMM_species.csv"
#GEMs_fp = "1_GEMs"
#sp2bin_fp="HQ_bins_with_compl.csv"
#community_types_fp = "5_MICOM/0_MAGs_tables"

if not os.path.exists(community_types_fp):
    os.makedirs(community_types_fp)

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


###### read and parse KMA results
df = pd.read_csv(kma_res_fp, index_col = 0, encoding='latin1')
df['binID'] = df.index.map(sp2bin)

columns = df.columns.tolist()[:-1] # sample names


######  make a MAG table with all community types and all samples:
d = {} # the dictionary to pass to pandas dataframe (needed because appending a series in the loop is too slow)
i=0 # a counter to use to add entries to "dict"
for index, row in df.iterrows():
    for c in columns:
        if df.at[index, c] > 0:
            bin_ID = df.at[index, 'binID']
            GEM_fp = GEMs_fp + "/" + bin_ID + ".xml"
            species = index
            abund = df.at[index, c]
            comm = 'sample_' + str(comm_types[c])
            d[i] = {"id":bin_ID, "species":species, "abundance":abund, "sample_id":comm, "file":GEM_fp}
            i = i + 1

MAG_table_all = pd.DataFrame.from_dict(d, "index")
#MAG_table_all.to_csv("all_samples_and_bins.csv") check file if needed


###### get average abundances
MAG_table_all_DMMs = MAG_table_all.groupby(['id','species','sample_id','file']).mean().reset_index()
MAG_table_all_DMMs = MAG_table_all_DMMs.round(2) # keep only 2 decimal places
#MAG_table_all_DMMs.to_csv("MAG_table_all_DMMs.csv") check file if needed


###### Separate communities and save it to different files
unique_comm_types = MAG_table_all_DMMs.sample_id.unique()

for biome in unique_comm_types:
    biome_table = MAG_table_all_DMMs[MAG_table_all_DMMs.sample_id == biome]
    table_fp = community_types_fp + "/" + biome + ".csv"
    biome_table.to_csv(table_fp, index=False)


print ("Done. Happy simulations.")