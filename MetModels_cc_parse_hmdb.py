#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to classify chemical classes (CC) - part 2
Extract the chemical class for all metabolites with an HMDB
Created on 26/3/21
@author: V.R.Marcelino
"""

# first produce a smaller file, removing all lines with 5+ spaces:
#grep -v '     ' hmdb_metabolites.xml > hmdb_metabolites_clean.xml

# exclude every line with more than 5 spaces - DONE
# parse every line of the filtered xml file
# if accession is in list of wanted accession -> get the closest class info
# stop search if "accesion" is found before class.

import re
import pandas as pd

in_bigg = "bigg_models_simplified.tsv"
in_hmdb = "hmdb_metabolites_clean.xml"

out_annotated_bigg = "bigg_models_w_classes.tsv"

# open and store the hmbd identifiers:
bigg_db = pd.read_csv(in_bigg, sep='\t')
wanted_hmdb_ids = list(bigg_db.hmbdID.unique())
wanted_hmdb_ids.remove('no_hmbdID_found')
wanted_hmdb_ids_xml = [">" + id +'</' for id in wanted_hmdb_ids]

#create a dictionary with HMDB_IDs and class, superclass and subclass:
hmdbID2class = {}
hmdbID2class['no_hmbdID_found'] = ['no super class','no class','no sub class']

line_counter = 0
end_of_last_metabolite = 0

with open (in_hmdb) as db:
    for line in db:
        line_counter += 1

        # check if it has an wanted accession
        if any(substring in line for substring in wanted_hmdb_ids_xml):
            accession=re.split('>|<',line)[2]
            accession_position = line_counter

        #check when it reaches the end of the metabolite
        if line.startswith('<metabolite>'):
            end_of_last_metabolite = line_counter - 1
            accession_position = 5000000 # ridiculously high number

        # find next "class" info:
        if '<super_class>' in line:
            super_class_id = re.split('>|<',line)[2]

            # check if it has a class in the next line:
            line = next(db)
            if '<class>' in line:
                classID = re.split('>|<',line)[2]
                line = next(db) # move to next line
            else:
                classID = "no class"

            #get sub class:
            if '<sub_class>' in line:
                sub_classID = re.split('>|<',line)[2]
            else:
                sub_classID = "no sub class"

            class_position = line_counter

            # check if the class is assigned to one of the wanted accessions, if yes - save it to dic:
            if accession_position < class_position > end_of_last_metabolite:
                print(super_class_id, '\n', classID, '\n', sub_classID, '\n')
                hmdbID2class[accession] = [super_class_id, classID,sub_classID]



#### now add it to the table:
bigg_db['Classes'] = bigg_db['hmbdID'].map(hmdbID2class)

#separate classes info into 3 columns:
new_bigg_db = bigg_db[['bigg_id', 'name', 'hmbdID']]
new_bigg_db[['SuperClass','Class','SubClass']] = bigg_db.Classes.apply(pd.Series)

### Save to file:
new_bigg_db.to_csv(out_annotated_bigg, sep='\t', index=False)


print("done!!")