#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to run MICOM's grow workflow for one community type (a.k.a. sample or biome)

Works with MICOM v 0.25.1

Created on 30/3/21
Last updated on 25/8/21

@author: V.R.Marcelino
"""

import pandas as pd
from micom import load_pickle
from micom.workflows import grow
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument('-c', '--comm_fp', help="""Path to folder containing community pickles""", required=False, default="1_communities/")
parser.add_argument('-s', '--sample', help="""sample, or community type, to be analysed""", required=True)
parser.add_argument('-t', '--trade_off', help="""trade_off to use in the grow workflow. Default here is 0.5, but the original default is 1""", required=False, default=0.5)
parser.add_argument('-th', '--threads', help="""threads to use""", required=False, default=1)
parser.add_argument('-o', '--out_folder', help="""output_folder. Default = 2_Exchanges""", required=False, default="2_TradeOffs")


args = parser.parse_args()

# community type and file_paths:
pickles_path = args.comm_fp
sample = args.sample
trade_off = args.trade_off
th = int(args.threads)
out_dir=args.out_folder

#pickles_path = '1_communities'
#sample = 'ERR589448'
#trade_off = 0.5
#th = 2

# Load community:
comm_fp = pickles_path +"/"+ sample + '.pickle'
com = load_pickle(comm_fp)


#############################
#### SIMULATE GROWTH
print ("\n simulating growth...\n")

#generate a manifest:
comm_filename = sample + '.pickle'
manifest = pd.DataFrame({"sample_id": sample, "file": comm_filename}, index=[0])

# Medium (western diet, after diluting nutrients absorbed in the small intestine)
med = pd.DataFrame(com.medium.items(), columns=['reaction', 'flux'])

# grow!! using parsimonious FBA
res = grow(manifest, pickles_path, medium=med, tradeoff=trade_off, threads=th, strategy="pFBA")
res.exchanges

## add sample name:
res.exchanges['sample_id'] = sample

## save to file:
out_fp = out_dir +"/"+ "exchanges_grow_" + sample + ".csv"
outfile=open(out_fp,"w")
res.exchanges.to_csv(outfile)
outfile.close()

print ("\nDONE!!\n")
