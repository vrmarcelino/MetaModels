#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to parse the mapping results of KMA, when the ref database is a set of VAMB bins

Created on 5/3/21
@author: V.R.Marcelino
"""

import pandas as pd
from argparse import ArgumentParser
import os


parser = ArgumentParser()
parser.add_argument('-i', '--input_fp', help="""Path to the folder containing KMA res files.""", required=True)
parser.add_argument('-o', '--output_fp', default = 'merged_samples.csv',
                    help='Path to the output file. Default = merged_samples', required=False)

args = parser.parse_args()
in_folder = args.input_fp
output = args.output_fp


#in_folder = "01_KMA"
#output = "merged_samples.csv"


### define how we will aggregate the info
f4agg = {}
f4agg['Depth']= 'sum'
f4agg['Taxonomy']= 'first'

# create new dataframe to store all samples:
all_samples = pd.DataFrame()

# read input files and merge results
for file in os.listdir(in_folder):
    if file.endswith(".res"):
        sample_name = file.split("_mags.res")[0]
        result_fp = os.path.join(in_folder, file)
        out_sample_fp = result_fp + ".summary.csv"

        df = pd.read_csv(result_fp, sep='\t', encoding='latin1')

        # Make 3 new columns for binID, contigIDs and taxonomy
        df[["BinID","contigID","Taxonomy"]] = df["#Template"].str.split('|', expand=True)

        depth_by_tax = df.groupby(by="BinID").agg(f4agg)
        depth_by_tax.rename(columns={'Depth': sample_name}, inplace=True)

        #save this file for reference
        pd.DataFrame.to_csv(depth_by_tax, out_sample_fp, index=True)

        all_samples = pd.concat([all_samples, depth_by_tax], sort=True, axis=1)

# Group taxon ranks
all_samples = all_samples.groupby(by=all_samples.columns, axis = 1).first()

# Fill NaN with zeros:
all_samples = all_samples.fillna(0)

# Save to file
pd.DataFrame.to_csv(all_samples, output, index=True)

