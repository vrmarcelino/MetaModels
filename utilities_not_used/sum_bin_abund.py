#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to summarize the abundance of bins (generated with VAMB) across metagenome samples
Run after MetModels_filter_clusters.py to use only the HQ bins

Retuns a Bin X Sample table (similar to a OTU table)

This script is not being used in our pipeline

Created on 2/3/21
@author: V.R.Marcelino
"""

import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-jgi', '--jgi_abund_fp', help="""Path to the jgi.abundance.dat file, 
                        produced with vamb""", required=True)
parser.add_argument('-c', '--clusters_fp', help="""Path to the clusters_filtered.tsv, 
                        produced with MetModels_filter_clusters.py""", required=True)
parser.add_argument('-o', '--output_fp', default = '3_Models/0_bin_table.tsv',
                    help='Path to the output file. Default = 3_Models/0_bin_table.tsv', required=False)

args = parser.parse_args()
jgi_abund_fp = args.jgi_abund_fp
wanted_clusters_fp = args.clusters_fp
processed_abund_fp = args.output_fp

#jgi_abund_fp="2_MAGs/1_jgi_matrix/jgi.abundance_test.dat.txt"
#wanted_clusters_fp = "2_MAGs/2_vamb/clusters_filtered_test.tsv"
#processed_abund_fp = "3_Models/0_bin_table.tsv"

## make a dictionary with contigs 2 clusters data:
ctg2bin = {}
with open (wanted_clusters_fp) as clusters:
    for line in clusters:
        bin = line.split("\t")[0]
        ctg = line.split("\t")[1].rstrip()
        ctg2bin[ctg] = bin


###### produce a filtered jgi_abundance file:
processed_abund_all_ctgs = processed_abund_fp + "_all_ctgs.tsv"
bin_table_temp = open(processed_abund_all_ctgs, 'w')

with open(jgi_abund_fp, 'r') as jgi:

    # write header
    first_line = jgi.readline()
    header = "binID\t" + first_line
    bin_table_temp.write(header)

    for line in jgi:
        contig_ID = line.split("\t")[0]
        if contig_ID in ctg2bin.keys():
            new_line = ctg2bin[contig_ID] + "\t" + line
            bin_table_temp.write(new_line)

bin_table_temp.close()

###### Merge abundances
### add option to filter low abundance as well?!
df = pd.read_csv(processed_abund_all_ctgs, sep='\t')

# merge abundances:
new_df = df.groupby(by="binID").agg(sum)

### Save
pd.DataFrame.to_csv(new_df,processed_abund_fp, sep='\t')

print ("Done")

