#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to run MICOM and build community models.
Also adds a media to the community and saves it as a pickle file.
Using here the Western Diet csv file, formated to work with carveme the MICOM website,
which also includes the dilution factors that lower import flux bounds for metabolites that are commonly
absorbed in the small intestine

Created on 06/05/21
@author: V.R.Marcelino
"""
import os
from micom import Community
import pandas as pd
from micom.qiime_formats import load_qiime_medium
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument('-s', '--sample', help="""sample, or community type, to be analysed. Must be present in the 
folder containing MAGs (or community) tables""", required=True)

parser.add_argument('-f', '--tables_fp', help="""Path to folder containing MAGs (or community) tables""", required=False, default="0_MAGs_tables")
parser.add_argument('-p', '--pickles', help="""Path to folders to store the community models""", required=False, default="1_communities")

parser.add_argument('-m', '--media', help="path to the media file", required=False, default="0_diet/western_diet_gut_carveme.qza")


args = parser.parse_args()


# community type and file_paths:
in_folder = args.tables_fp
comm_folder = args.pickles
sample = args.sample
w_media = args.media

#in_folder = '0_MAGs_tables'
#comm_folder = '1_communities'
#sample = 'test_sample'
#w_media = '0_diet/western_diet_gut_carveme.qza'

if not os.path.exists(comm_folder):
    os.mkdir(comm_folder)

### import MAGs table containing genome-scale model paths
## in the tutorial, they call this table 'taxonomy'
fp = in_folder + "/" + sample + ".csv"
mag_tb = pd.read_csv(fp)

### import and parse media - western diet
medium = load_qiime_medium(w_media)


##############################
# BUILD COMMUNITY MODELS

# In order to convert the specification in a community model we will use the Community class from micom
# which derives from the cobrapy Model class.
# this took ~45 sec for 10 MAGs
print ("\nBuilding community, be patient...\n")
com = Community(mag_tb)
print("Done. Built a community with a total of {} reactions.\n".format(len(com.reactions)))


##### Add media to the file
print ("Adding media and saving to pickle\n")
# check if the names match & add to the community object:
ex_ids = [r.id for r in com.exchanges]

# replace _m with _e_m (see what needs to change with the --fbc2 models)
medium = medium.replace(regex=r'_m$', value='_e_m')

medium['reaction'].isin(ex_ids).sum() # must be a large number of reactions (>100)
medium.index = medium.index.str.replace(r'_m$','_e_m', regex=True)
med = medium[medium.index.isin(ex_ids)] # exclude medium items not used by the microbiome

diet = med.flux * med.dilution # dilute nutrients absorbed in the small intestine


# This only affects the com.cooperative_tradeoff (the grow workflow is not affected)
com.medium = diet


# save this community to file (Can I save one pickle file per GEM?)
comm_fp = comm_folder + "/" + sample + ".pickle"
com.to_pickle(comm_fp)

print ("\nDone building community model for %s!\n"%(sample))


