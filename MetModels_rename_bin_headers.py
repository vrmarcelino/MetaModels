#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to add binID and taxonomy to sequence headers to metagenome bins created with VAMB

Also converts to sequential fasta

Created on Mar 4 2021

@author: V.R.Marcelino
"""

from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument('-i', '--input_fasta', help="""Path to the input fasta file that will be renamed""", required=True)
parser.add_argument('-c', '--clusters_fp', help="""Path to the clusters_filtered.tsv, 
                        produced with MetModels_filter_clusters.py""", required=True)
parser.add_argument('-g', '--gtdb_taxonomy_fp', help="""Path to file containing binID to taxonomy map.""", required=True)
parser.add_argument('-o', '--output_fp', default = 'merged_renamed_db',
                    help="""Path to the output file. Default='renamed_db.fas'""", required=False)

args = parser.parse_args()
in_fasta = args.input_fasta
wanted_clusters_fp = args.clusters_fp
tax_fp = args.gtdb_taxonomy_fp
output = args.output_fp



#### input files:
#in_fasta = "merged_HQ_Bins.fna"
#wanted_clusters_fp = "clusters_filtered_test.tsv"
#tax_fp = "gtdb_bac_and_arch.tsv"
#output = "merged_renamed_db.fas"

### Define a class object for each ctg header, which stores binID and taxonomy
class CtgInfo():
    def __init__(self, contigID=None, full_ctg_name=None, binID=None, taxonomy=None):

        self.contigID = contigID
        self.full_ctg_name = full_ctg_name
        self.binID = binID
        self.taxonomy = taxonomy

### store info in a list of CtgInfo objects
ctgs_info = []
with open (wanted_clusters_fp) as clusters:
    for line in clusters:
        bin = line.split("\t")[0]
        full_ctg_name = line.split("\t")[1].rstrip()
        contigID = full_ctg_name.split("-")[0]

        ci = CtgInfo()
        ci.contigID = contigID
        ci.full_ctg_name = full_ctg_name
        ci.binID = bin
        ctgs_info.append(ci)

### Store taxonomy in ctgs_info
with open (tax_fp) as tax:
    next(tax) # skip header
    for line in tax:
        bin = line.split("\t")[0]
        taxonomy = line.split("\t")[1]
        taxonomy = taxonomy.replace(" ", "_") # removes spaces from species names.

        for c in ctgs_info:
            if c.binID == bin:
                c.taxonomy = taxonomy


### function that takes as input a ctg header, and outputs the wanted ctg header format
def beheader(old_header, list_of_ctgsInf, counter):
    for x in list_of_ctgsInf:
        if x.full_ctg_name == old_header:
            new_header = x.binID + "|" + x.contigID + "|" + x.taxonomy
            counter += 1
            break
        else:
            new_header = "not_renamed_" + old_header + "\n"
    return(new_header, counter)


### rename  contig headers, store in a new file:
renamed_db = open(output, 'w')
ctgs = 0
ctgs_renamed = 0
with open(in_fasta) as genome:
    for line in genome:
        if line.startswith(">"):
            ctgs += 1
            contig_name = line.split(">")[1].strip()
            rename_seq = beheader(contig_name, ctgs_info, ctgs_renamed)
            if ctgs == 1: # does not need a new line in the beginning of the file
                line = ">" + rename_seq[0]
            else:
                line = "\n>" + rename_seq[0]
            ctgs_renamed = rename_seq[1]

        else:
            line = line.replace("\n", "") # convert interleaved to sequential

        renamed_db.write(line)

renamed_db.close()

ctgs_not_found = ctgs - ctgs_renamed

print ("Done! %i contigs renamed, and %i contigs were not found" %(ctgs_renamed,ctgs_not_found))

