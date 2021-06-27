#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scrip to filter samples according to number of reads
and unwanted health status (overweight and underweight)

Produces a filtered metadata table

Created on 24/6/21
@author: V.R.Marcelino
"""

import pandas as pd

in_metadata_fp = "0_wanted_samples_all.csv"
in_accessions_fp = "accession2filenames_batches123.csv"
in_lines = "sample_lines.csv" # file with number of lines (reads*4) per sample
out_fp = "filtered_metadata.csv"

min_reads = 15000000 # 15M

in_metad = pd.read_csv(in_metadata_fp)
in_acc = pd.read_csv(in_accessions_fp)
in_reads = pd.read_csv(in_lines)

#### Filter sample table according to depth threshold:
in_reads_filt = in_reads.loc[in_reads['Lines'] > min_reads * 4]
wanted_samples = in_reads_filt['Sample'].tolist()
# remove the '_cat' from run_acc:
wanted_samples = [x.replace("_cat", "") for x in wanted_samples]


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

