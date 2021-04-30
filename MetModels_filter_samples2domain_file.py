#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to filter teh samples2domain file, so we run carveme on the wanted bins only
(i.e. one bin per species)

Created on 26/4/21
@author: V.R.Marcelino
"""

in_HQ_bins = "HQ_bins_with_compl.csv"
in_samples2domain = "samples2domain.tsv"
out_samples2domain = "samples2domain_one_per_sp.tsv"

### get a list of wanted bins:
wanted_bins = []
with open(in_HQ_bins) as f:
    next(f) # skip first line
    for line in f:
        bin = line.split(",")[0]
        wanted_bins.append(bin)

### save wanted bins:
out = open(out_samples2domain, 'w')
with open(in_samples2domain) as s:
    for line in s:
        bin = line.split("\t")[0]
        if bin in wanted_bins:
            print (bin)
            out.write(line)

out.close()

