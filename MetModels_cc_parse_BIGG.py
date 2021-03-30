#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to classify chemical classes (CC) - part 1
Parse the bigg_models_metabolites_raw.txt file, extracting the HMDB ID for each metabolite
Created on 26/3/21
@author: V.R.Marcelino
"""

import re

in_file = "bigg_models_metabolites_raw.txt"
out_file = "bigg_models_simplified.tsv"

out = open(out_file, 'w')
out.write("bigg_id\tname\thmbdID\n") # header


with open (in_file, 'r') as bigg:
    next(bigg)
    for line in bigg:
        splitted_line = line.split('\t')
        metabolite_abrev = splitted_line[0]
        name = splitted_line[2]
        find_hmdb = re.search('hmdb/(.*?); ', line)
        if find_hmdb != None:
            hmdbID = find_hmdb.group(1)
        else:
            hmdbID = "no_hmbdID_found"

        new_line = metabolite_abrev + '\t' + name + '\t' + hmdbID + '\n'
        out.write(new_line)
out.close()
