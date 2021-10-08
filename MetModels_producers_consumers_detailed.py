#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to identify the producers and consumers of each metabolite
-- focus on the healthy population only!

Created on 6/10/21
@author: V.R.Marcelino
"""

from argparse import ArgumentParser
import os, glob, csv
import pandas as pd

parser = ArgumentParser()
parser.add_argument('-f', '--folder_w_exchange_files', help="""Path to the folder containing exchange files produced by MICOM""", required=True)
parser.add_argument('-m', '--metadata', help="""Path to the metadata file""",required=True)
parser.add_argument('-s', '--spp_classification', help="""Path to the wanted_spp_classification.tsv file""",required=True)
parser.add_argument('-o', '--output', help="""name of output file. Default = producers_consumers_detailed.csv""", required=False, default="producers_consumers.csv")

args = parser.parse_args()
in_path = args.folder_w_exchange_files
metad_fp = args.metadata
binID2spp = args.spp_classification
out_file = args.output


#in_path = "2_exchanges"
#out_file = "producers_consumers_detailed.csv"
#metad_fp = "metadata_rewiring_microbiome.csv"
#binID2spp = "wanted_spp_classification.tsv"

## read metadata
metad_all_samples = pd.read_csv(metad_fp)
metad_all_samples.rename(columns={'Sample':'sample_id'}, inplace=True)

## read bins 2 spp classification map:
binID2taxa = {}
with open(binID2spp, mode='r') as inp:
    reader = csv.reader(inp, delimiter="\t")
    binID2taxa = {rows[0].replace(".","_"):rows[1] for rows in reader}



## merge exchange files
all_files = glob.glob(os.path.join(in_path, "exchanges_grow_*.csv"))
df_from_each_file = (pd.read_csv(f, sep=',') for f in all_files)
df_merged = pd.concat(df_from_each_file, ignore_index=True)


## remove media:
exch_df = df_merged.loc[df_merged['taxon'] != "medium"]

## remove the "cat" from sample names:
pd.options.mode.chained_assignment = None  # disables the warning, default='warn'
exch_df["sample_id"] = exch_df["sample_id"].str.replace("_cat", "")

## add 'HD' info and keep only healthy individuals:
exch_df_metad = pd.merge(exch_df, metad_all_samples[["sample_id","HD"]], on="sample_id", how="left")
exch_df_healthy = exch_df_metad[exch_df_metad.HD == "healthy"]


## separate consumption and production fluxes into two columns:
mask = exch_df_healthy['flux'] < 0
exch_df_healthy['flux_production'] = exch_df_healthy['flux'].mask(mask)
exch_df_healthy['flux_consumption'] = exch_df_healthy['flux'].mask(~mask)
exch_df_healthy = exch_df_healthy.fillna(0) # replace NaNs with zeros

## also separate consumers and producer MAGs into two columns:
exch_df_healthy['taxon_production'] = exch_df_healthy['taxon'].mask(mask)
exch_df_healthy['taxon_consumption'] = exch_df_healthy['taxon'].mask(~mask)


# calculate number of producers/consumers per metabolite, and the sum of their fluxes:
count_producers = exch_df_healthy.groupby(["sample_id", "metabolite"]).apply(lambda df: df.flux_production.astype(bool).sum(axis=0))
sum_product_flux = exch_df_healthy.groupby(["sample_id", "metabolite"]).apply(lambda df: sum(df.flux_production))
producer_MAGs = exch_df_healthy.groupby(["sample_id", "metabolite"])['taxon_production'].agg(['unique']).squeeze() #the .squeeze converts dataframe to series


count_consumers = exch_df_healthy.groupby(["sample_id", "metabolite"]).apply(lambda df: df.flux_consumption.astype(bool).sum(axis=0))
sum_consump_flux = exch_df_healthy.groupby(["sample_id", "metabolite"]).apply(lambda df: sum(df.flux_consumption))
consumer_MAGs = exch_df_healthy.groupby(["sample_id", "metabolite"])['taxon_consumption'].agg(['unique']).squeeze() #the .squeeze converts dataframe to series


## merge the dataframes:
prod_con_summary = pd.concat([count_producers, sum_product_flux,count_consumers, sum_consump_flux,producer_MAGs, consumer_MAGs], axis=1)
prod_con_summary.columns = ["n_producers", "flux_producers","n_consumers","flux_consumers","producer_MAGs", "consumer_MAGs"]


## calculate importance score:
prod_con_summary['importance_score'] = 2 * ((prod_con_summary['n_producers'] * prod_con_summary['n_consumers']) / (prod_con_summary['n_producers'] + prod_con_summary['n_consumers']))


######## group by metabolite, keeping track of the producer and consumer bins:

# prepare dataframe:
prod_con_summary.reset_index(inplace=True) # convert multi-index to columns

prod_con_summary['producer_MAGs_str'] = ["".join(item) for item in prod_con_summary['producer_MAGs'].astype(str)] # convert to string
prod_con_summary['consumer_MAGs_str'] = ["".join(item) for item in prod_con_summary['consumer_MAGs'].astype(str)]


# aggregate dataset by metabolite
fun4agg = {'importance_score':['mean', 'std'], 'flux_producers':['mean', 'std'], 'flux_consumers':['mean', 'std'],
           'n_producers': ['sum', 'mean', 'std'], 'n_consumers':['sum', 'mean', 'std'],
           'producer_MAGs_str':['sum'], 'consumer_MAGs_str': ['sum']}

prod_con_samples_agg = prod_con_summary.groupby('metabolite').agg(fun4agg)

# merge multilevel column names
prod_con_samples_agg.columns = ['_'.join(col) for col in prod_con_samples_agg.columns]

####### clean and convert consumer_MAG column to a meaningful list:
# first replace all patterns with a space (otherwise we may concatenate bin names), then replace 2+ spaces with a single space
pattern2replace = '|'.join(["'", "\[", "\]", "nan"])

### producers:
prod_con_samples_agg['producer_MAGs_str_sum'] = prod_con_samples_agg['producer_MAGs_str_sum'].str.replace(pattern2replace, ' ')
prod_con_samples_agg['producer_MAGs_str_sum'] = prod_con_samples_agg.producer_MAGs_str_sum.str.replace('\s{2,}', ' ') # replaces any number of spaces by just one

# convert to list:
prod_con_samples_agg['producer_MAGs_str_sum'] = prod_con_samples_agg.producer_MAGs_str_sum.apply(lambda x: x[1:-1].split(' '))

# replace empty lists with None (otherwise they are counted late ron)
prod_con_samples_agg['producer_MAGs_str_sum'] = prod_con_samples_agg.producer_MAGs_str_sum.apply(lambda x: list(filter(None, x)))

# remove duplicates
prod_con_samples_agg['producer_MAGs_unique'] = prod_con_samples_agg['producer_MAGs_str_sum'].apply(lambda x: list(set(x)))


### consumers:
prod_con_samples_agg['consumer_MAGs_str_sum'] = prod_con_samples_agg['consumer_MAGs_str_sum'].str.replace(pattern2replace, ' ')
prod_con_samples_agg['consumer_MAGs_str_sum'] = prod_con_samples_agg.consumer_MAGs_str_sum.str.replace('\s{2,}', ' ') # replaces any number of spaces by just one

# convert to list:
prod_con_samples_agg['consumer_MAGs_str_sum'] = prod_con_samples_agg.consumer_MAGs_str_sum.apply(lambda x: x[1:-1].split(' '))

# replace empty lists with None (otherwise they are counted late ron)
prod_con_samples_agg['consumer_MAGs_str_sum'] = prod_con_samples_agg.consumer_MAGs_str_sum.apply(lambda x: list(filter(None, x)))

# remove duplicates
prod_con_samples_agg['consumer_MAGs_unique'] = prod_con_samples_agg['consumer_MAGs_str_sum'].apply(lambda x: list(set(x)))


### identify flexible MAGs (those that can be producers & consumers)

# helper functions:
def find_duplicates(l1, l2):
    """returns MAGs that can be producers and consumers"""
    dups = list(set(l1).intersection(l2))
    return(dups)

def find_exclusives(lx, l_dups):
    """returns MAGs that are only present in lx, not in the duplicate MAG list"""
    exclusive_in_lx = list(set(lx).symmetric_difference(set(l_dups)))
    return(exclusive_in_lx)

# calc:
prod_con_samples_agg['flexible_MAGs'] = prod_con_samples_agg.apply(lambda x: find_duplicates(x.producer_MAGs_unique, x.consumer_MAGs_unique), axis=1)
prod_con_samples_agg['exclusive_producers'] = prod_con_samples_agg.apply(lambda x: find_exclusives(x.producer_MAGs_unique, x.flexible_MAGs), axis=1)
prod_con_samples_agg['exclusive_consumers'] = prod_con_samples_agg.apply(lambda x: find_exclusives(x.consumer_MAGs_unique, x.flexible_MAGs), axis=1)


### calculate % of each category (makes it easier to plot in R):
def calc_total(l1,l2,l3):
    total = len(l1) + len(l2) + len(l3)
    return (total)

def calc_proportion(lx, n_total):
    prop = len(lx) * 100 / n_total
    return (prop)



prod_con_samples_agg['total_n_unique_MAGs'] = prod_con_samples_agg.apply(lambda x: calc_total(x.flexible_MAGs, x.exclusive_producers, x.exclusive_consumers), axis=1)
prod_con_samples_agg['perc_flexible'] = prod_con_samples_agg.apply(lambda x: calc_proportion(x.flexible_MAGs, x.total_n_unique_MAGs), axis=1)
prod_con_samples_agg['perc_producers'] = prod_con_samples_agg.apply(lambda x: calc_proportion(x.exclusive_producers, x.total_n_unique_MAGs), axis=1)
prod_con_samples_agg['perc_consumers'] = prod_con_samples_agg.apply(lambda x: calc_proportion(x.exclusive_consumers, x.total_n_unique_MAGs), axis=1)


### add species classification:
def translate_bin2spp(l1, dict):
    spp_list = [dict[x] for x in l1]
    return (spp_list)

prod_con_samples_agg['flexible_MAGs_taxa'] = prod_con_samples_agg.apply(lambda x: translate_bin2spp(x.flexible_MAGs,binID2taxa), axis=1)
prod_con_samples_agg['exclusive_producers_taxa'] = prod_con_samples_agg.apply(lambda x: translate_bin2spp(x.exclusive_producers,binID2taxa), axis=1)
prod_con_samples_agg['exclusive_consumers_taxa'] = prod_con_samples_agg.apply(lambda x: translate_bin2spp(x.exclusive_consumers,binID2taxa), axis=1)

## save large file:
prod_con_samples_agg.to_csv(out_file)

######## save a simplified/shorter version of the table:
out_simplified = out_file.replace(".csv", "_simplified_meanMIS_grater_than_5.csv")

#remove rows (metabolites) with mean importance score smaller than
prod_con_samples_agg_simpler = prod_con_samples_agg[prod_con_samples_agg['importance_score_mean'] > 5]
prod_con_samples_agg_simpler = prod_con_samples_agg_simpler[['importance_score_mean','importance_score_std','flexible_MAGs_taxa', 'exclusive_producers_taxa','exclusive_consumers_taxa', 'perc_flexible','perc_producers','perc_consumers']]


## save large file:
prod_con_samples_agg_simpler.to_csv(out_simplified)


print ("""\n\nDONE!!
    Detailed output saved as %s
    Simplified output saved as %s \n""" %(out_file, out_simplified))


