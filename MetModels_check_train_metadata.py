#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scrip to get some stats and filter samples acording to max number of samples per study

Produces a filtered metadata table

Created on11/07/21
@author: V.R.Marcelino
"""

import pandas as pd

in_metadata_fp = "train_datatset_to_filter.csv"
out_fp = "filtered_metadata.csv"

in_metad = pd.read_csv(in_metadata_fp)


### print number of samples per study:
studies = in_metad.Author_Year.unique()
n_studies = len(studies)
print ("number of studies: %i" %(n_studies))


# get a table with samples per health status:
in_metad['health_status'].value_counts()

# get a table with samples per study:
in_metad['Author_Year'].value_counts()






### make a dictionary with run_acc to BioSample
in_acc.index = in_acc['SAMD00114718']
in_acc_d = in_acc.T.to_dict('list') # needed to transpose the dataframe

### Filter the dictionary based on the wanted accesion lists, save BioSample:
biosample =[]

for key,value in in_acc_d.items():
    for acc in value:
        if acc in wanted_samples:
            biosample.append(key)

# remove duplicates:
list_set = set(biosample)
biosample_uniq= list(list_set)

### Filter the metadata table:
filt_metad = in_metad[in_metad['BioSample'].isin(biosample_uniq)]
n_samples_remaining = len(filt_metad.index)
print ("number of samples after filtering by minimum number of reads: %i" %(n_samples_remaining))

### Filter overweight and underweight (intermediary categories between healthy and disease)
filt_metad = filt_metad[filt_metad['health_status'] != "overweight"]
filt_metad = filt_metad[filt_metad['health_status'] != "underweight"]
n_samples_remaining = len(filt_metad.index)
print ("number of samples after filtering overweight and underweight individuals: %i" %(n_samples_remaining))


### print some stats:
studies = filt_metad.Author_Year.unique()
n_studies = len(studies)
print ("number of studies: %i" %(n_studies))


# get a table with samples per health status:
filt_metad['health_status'].value_counts()

# get a table with samples per study:
filt_metad['Author_Year'].value_counts()

####################################
### save a list of wanted files:
wanted_biosamples = filt_metad['BioSample'].tolist()
filtered_files_prefix = {}

for key,value in in_acc_d.items():
    for acc in value:
        if acc in wanted_samples:
            if key in wanted_biosamples:
                filtered_files_prefix[acc] = key


# add this info to the metadata table:
filtered_files_prefix_inv = {v: k for k, v in filtered_files_prefix.items()} # invert keys/items
filt_metad['file_prefix'] = filt_metad['BioSample'].map(filtered_files_prefix_inv)

# save to file:
filt_metad.to_csv(out_fp, index=False)

