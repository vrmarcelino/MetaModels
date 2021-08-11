#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that calculates the number AND quality of "donor" links for each bacteria,
based on the mininal_fluxes_exchange file from MICOM (tradeoffs workflow)

The idea is to test if healthy microbiomes are enriched with bacteria with a higher number of 'donor' links.

Quality of 'donor' links can be measured by: export flux from bacteria * centrality of metabolite
Centrality of metabolite takes into consideration fwd fluxes only (reflecting their uptake by other bacteria).

### Example:
python 3_MetModels_calc_donors.py -i 0_input/mininal_fluxes_exchange_merged.csv -m 0_input/metadata_diseased_healthy_tests.csv -o 2_Donor_links/1_donor_links.csv
###

Created on 5/8/21
@author: V.R.Marcelino
"""

from argparse import ArgumentParser
import pandas as pd

parser = ArgumentParser()
parser.add_argument('-i', '--minimal_fluxes', help="""Path to the mininal_fluxes_exchange file produced with MetModels_merge_exchange_tables.py""", required=True)
parser.add_argument('-m', '--metadata', help="""path to the metadata, indicating if samples are healthy or diseased""", required=True)
parser.add_argument('-o', '--output', help="""name of the output file""", required=False, default="Donor_links.csv")

args = parser.parse_args()
min_fluxes_fp = args.minimal_fluxes
metadata_fp = args.metadata
out_file = args.output

# min_fluxes_fp = "0_input/mininal_fluxes_exchange_merged.csv"
# metadata_fp = "0_input/metadata_diseased_healthy_tests.csv"
# out_file = "Donor_links.csv"

################################################
################  Functions     ################

def calc_met_centrality(sampleX_df):
# function takes as input a pd.dataframe with the exchanges within a single sample, and
# calculates the centrality OF EACH METABOLITE
# by counting the occurrence of negative fluxes (i.e. how many organisms uptake the metabolite).
    centrality_dict ={}
    for metab in sampleX_df.columns:
        if metab.startswith("EX_"): # check if this is a metabolite:
            number_uptakes = (sampleX_df[metab] < 0).sum()
            if number_uptakes > 0: # if it is uptaken by any microbe, store in dict
                centrality_dict[metab] = number_uptakes
    return (centrality_dict)


def calc_donor_links (sampleX_df, sample_name, health_status):
# function takes as input a pd.dataframe with the exchanges within a single sample,
# calculates n. of donor links and the sum of centrality scores of all metabolites produced by each bin (== BIN Centrality)
# returns a dataframe with the info required for statistical analyses
    df = pd.DataFrame()
    metab_centrality = calc_met_centrality(sampleX_df)
    count_bins = 0

    for index, row in sampleX_df.iterrows():
        binID = row['compartment']
        count_donor_links = 0
        calc_bin_centrality = [] # sum of links
        calc_weighted_bin_centrality = [] # sum of links multiplied by flux
        #count donor links:
        for met in row.index:
            if met in metab_centrality.keys(): # check if metabolite is consumed by other species
                flux = row[met]
                if flux > 0: #check if metabolite is exported (positive flux) by this bin
                    count_donor_links += 1
                    calc_bin_centrality.append(metab_centrality[met])

                    weighted_centrality = metab_centrality[met] * flux
                    calc_weighted_bin_centrality.append(weighted_centrality)

        # calculate bin centrality and number of bins:
        bin_centrality = sum(calc_bin_centrality)
        weighted_bin_centrality = sum(calc_weighted_bin_centrality)
        count_bins += 1
        new_row = [[sample_name, binID,count_donor_links,bin_centrality, weighted_bin_centrality, health_status]]

        df = df.append(new_row)
    df['n_bins'] = count_bins
    return (df)




##################################################
################  Run analyses    ################

# read files
df_fluxes_all = pd.read_csv(min_fluxes_fp , sep=',')
metadata = pd.read_csv(metadata_fp, sep=',')

# add metadata to dict:
meta_dic = pd.Series(metadata.HD.values,index=metadata.Sample).to_dict()

result_df = pd.DataFrame()

## loop trough fluxes dataframe (via dict), subseting one sample for processing at a time:
# slow - took 30sec for 42 samples
for sample,value in meta_dic.items():

    wanted_rows = df_fluxes_all.loc[df_fluxes_all['sample'] == sample]
    health_status = value

    sampleX_links = calc_donor_links(wanted_rows,sample,health_status)

    result_df = result_df.append(sampleX_links)

# rename columns
result_df.columns = ["sample", "binID", "n_donor_links","bin_centrality","weighted_bin_centrality","HD","n_bins"]

# calculate 'donor score' (n donor links X bin centrality)
result_df['donor_score'] = result_df['n_donor_links'] * result_df['bin_centrality']
result_df['weighted_donor_score'] = result_df['n_donor_links'] * result_df['weighted_bin_centrality'] # taking into acocunt export fluxes

result_df.to_csv(out_file, index=False)

print ("Done!")

