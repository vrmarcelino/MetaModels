#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to extract minimal growth media and optimize the cooperative tradeoff
Run this after MICOM_build_comm_models.py
Created on 30/3/21
edited: 19/05/21
@author: V.R.Marcelino

Script adapted from C. Diener's scripts:
https://github.com/micom-dev/paper
"""

#imports
import sys
import os
from argparse import ArgumentParser
import pandas as pd
from micom import load_pickle
from micom.media import minimal_medium
from micom.workflows import workflow
from micom.logger import logger
import time


parser = ArgumentParser()
parser.add_argument('-sl', '--sample_list', help="""path to a txt file containing sample.pickle names""", required=False)
parser.add_argument('-sn', '--sample_name', help="""sample name (e.g. SRR413605.pickle). Required when doing analyses with one sample at time (e.g. with snakemake)""", required=False)
parser.add_argument('-p', '--pickles_path', help="""path to the folder containing the pickle files. Default = 1_communities""", required=False, default="1_communities")
parser.add_argument('-trad', '--trade_off', help="""trade_off (fraction) to use in the cooperative tradeoff""", required=False, default=0.5)
parser.add_argument('-t', '--threads', help="""number of threads. Default = 2""", required=False, default=2)
parser.add_argument('-o', '--out_folder', help="""output_folder. Default = 2_TradeOffs""", required=False, default="2_TradeOffs")


args = parser.parse_args()

samples_list_fp = args.sample_list
samples_name_fp = args.sample_name
pickles_path = args.pickles_path
trade_off = args.trade_off # fraction
max_procs = int(args.threads)
out_dir=args.out_folder


#samples_fp='1_communities/individual_ERR_samples_test.txt'
#pickles_path = '1_communities'
#trade_off = 0.5
#max_procs = 2
#out_dir = "2_TradeOffs"

# check if sample list or name were given:
if samples_list_fp != None:
    ## read sample names to a list
    with open(samples_list_fp) as f:
        samples = f.read().splitlines()
    if 'samples.txt' in samples:
        samples.remove('samples.txt') #there if the list was created with ls -1 > samples.txt

elif samples_name_fp != None:
    samples = [samples_name_fp]
    print (samples)

else:
    print ("Error: please provide --sample_list or --sample_name.")
    sys.exit()


## make out dir
if not os.path.exists(out_dir):
    os.mkdir(out_dir)



## function that optimizes the cooperative tradeoff, first using the western media for upper boundaries,
## then using the minimal media to get the metabolic exchanges
def media_and_gcs(sam):

    com = load_pickle(pickles_path +"/"+ sam)

    # Get growth rates
    try:
        sol = com.cooperative_tradeoff(fraction=trade_off)
        rates = sol.members["growth_rate"].copy()
        rates["community"] = sol.growth_rate
        rates.name = sam
    except Exception:
        logger.warning("Could not solve cooperative tradeoff for %s." % sam)
        return None

    # Get the minimal medium
    med = minimal_medium(com, 0.95 * sol.growth_rate, exports=True)
    med.name = sam

    # Apply medium and reoptimize
    com.medium = med[med > 0]
    sol = com.cooperative_tradeoff(fraction=0.5, fluxes=True, pfba=False) # uses the 'classic' FBA instead of the parsimonious FBA
    fluxes = sol.fluxes
    fluxes["sample"] = sam
    return {"medium": med, "gcs": rates, "fluxes": fluxes}


gcs = pd.DataFrame()
media = pd.DataFrame()
fluxes = pd.DataFrame()

#run multiple samples in parallel with micom.workflow
print ("\n simulating cooperative trade-off\n")
results = workflow(media_and_gcs, samples, max_procs)

for r in results:
    gcs = gcs.append(r["gcs"])
    media = media.append(r["medium"])
    fluxes = fluxes.append(r["fluxes"])

# get timestamp:
ts = str(round(time.time()))
print ("timestamp: %s"%(ts))


gcs_fp = out_dir + "/growth_rates_" + ts + ".csv"
media_fp = out_dir + "/minimal_imports_" + ts + ".csv"
fluxes_fp = out_dir + "/minimal_fluxes_all_" + ts + ".csv"

gcs.to_csv(gcs_fp)
media.to_csv(media_fp)
fluxes.to_csv(fluxes_fp)

### Get only the flux of the exchange_reactions:
# exchange reactions are reactions that move metabolites across in silico compartments.
ex_fluxes_fp = out_dir + "/minimal_fluxes_exchange_" + ts + ".csv"
ex_flux = fluxes.filter(regex='^EX_') # get only exchanges (starts with 'EX_')
ex_flux = ex_flux.filter(regex='e$') # remove media (media ends with 'e_m', so I want the ones that end with 'e' only)
ex_flux = ex_flux.fillna(0) #fill NANs with zeros

ex_flux = ex_flux.loc[:, (ex_flux != 0).any(axis=0)] #remove columns with all zeros
ex_flux['sample'] = fluxes['sample']
ex_flux = ex_flux.drop(index='medium') # remove medium

ex_flux.to_csv(ex_fluxes_fp)

print ("\nDONE! Exchange fluxes saved to %s\n"%(ex_fluxes_fp))





