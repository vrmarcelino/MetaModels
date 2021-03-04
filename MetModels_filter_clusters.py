#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to get a list of High Quality MAGs (bins) and filter the cluster.tsv file from VAMB
Created on 2/3/21
@author: V.R.Marcelino
"""
import os


bins_HQ_fp = "2_MAGs/3_bins_HQ"
clusters_fp = "2_MAGs/2_vamb/clusters.tsv"
wanted_clusters_fp = "2_MAGs/2_vamb/clusters_filtered.tsv"

# list the high quality bins (use f.path if want the full directory path)
hq_bins = [ f.name for f in os.scandir(bins_HQ_fp) if f.is_dir() ]


# save wanted lines in a new file:
wanted_clusters = open(wanted_clusters_fp, 'w')

with open(clusters_fp,'r') as clusters:
    for line in clusters:
        if line.split("\t")[0] in hq_bins:
            wanted_clusters.write(line)

wanted_clusters.close()

print ("done!")

