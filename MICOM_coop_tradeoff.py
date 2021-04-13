#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to run MICOM for one community type (a.k.a. sample or biome)
Specifically - run the 'grow' workflow
Later - test the cooperative trade-off

Created on 30/3/21
edited: 08/03/21
@author: V.R.Marcelino
"""

from micom import Community
import pandas as pd
from micom.qiime_formats import load_qiime_medium
from micom.workflows import grow
from micom.workflows import build
from micom import load_pickle
import os
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument('-f', '--folder_fp', help="""Path to folder containing MAG tables""", required=False, default="0_MAGs_tables/")
parser.add_argument('-s', '--sample', help="""sample, or community type, to be analysed""", required=True)
parser.add_argument('-t', '--trade_off', help="""trade_off (fration) to use in the cooperative tradeoff""", required=False, default=0.5)

args = parser.parse_args()


# community type and file_paths:
in_folder = args.folder_fp
sample = args.sample
trade_off = args.trade_off # fraction


#in_folder = '0_MAGs_tables/'
#sample = 'sample_1'

# folder to store intermediate files
model_path = "1_interm_files"
#os.mkdir(model_path)


### import MAGs table containing genome-scale model paths
## in the tutorial, they call this table 'taxonomy'
fp = in_folder + sample + ".csv"
mag_tb = pd.read_csv(fp)

#convert to relative abundances (check whether it is necessary?!)
mag_tb['abundance'] = mag_tb['abundance'].div(mag_tb['abundance'].sum())


##############################
# BUILD COMMUNITY MODELS

# In order to convert the specification in a community model we will use the Community class from micom
# which derives from the cobrapy Model class.
# this took ~45 sec for 10 MAGs

#com = Community(mag_tb)
#print("Built a community with a total of {} reactions.".format(len(com.reactions)))



# add media - using the western diet as a media:
#medium = load_qiime_medium("0_diet/western_diet_gut_carveme.qza")
# replace _m with _e_m (see what needs to change with the --fbc2 models)
#medium = medium.replace(regex=r'_m$', value='_e_m')

# check if the names match & add to the community object:
#ex_ids = [r.id for r in com.exchanges]
#medium['reaction'].isin(ex_ids).sum() # must be a large number of reactions (>100)
#medium.index = medium.index.str.replace(r'_m$','_e_m', regex=True)
#med = medium[medium.index.isin(ex_ids)] # exclude medium items not used by the microbiome

# check if it is the flux that should be used as media in the community object.
# This only affects the com.cooperative_tradeoff (the grow workflow is not affected)
#com.medium = med['flux']


# save this community to file (Can I save one pickle file per GEM?)
comm_fp = model_path + "/" + sample + "_community.pickle"
#com.to_pickle(comm_fp)
#load it:
com = load_pickle(comm_fp)
# see exchanges:
#com.exchanges

print("Loaded a community with a total of {} reactions.".format(len(com.reactions)))


#############################
# Cooperative trade-off

print ("\n simulating cooperative trade-off\n")

sol = com.cooperative_tradeoff(fraction=trade_off, fluxes=True, pfba=True)
sol

####### from Basile tutorial - save file:
matrix=sol.fluxes
matrix1=matrix.filter(regex='^EX_')
out_fp = "coop_tradeoff_fluxes_" + sample + ".csv"
outfile=open(out_fp,"w")
matrix1.to_csv(outfile)
outfile.close()



print ("\nDONE!!\n")





