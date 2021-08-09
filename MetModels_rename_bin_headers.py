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

### make a dictionary with contigs 2 clusters data:
ctg2bin = {}
with open (wanted_clusters_fp) as clusters:
    for line in clusters:
        bin = line.split("\t")[0]
        ctg = line.split("\t")[1].rstrip()
        ctg2bin[ctg] = bin


### make another dict bin2tax:
bin2tax = {}
with open (tax_fp) as tax:
    next(tax) # skip header
    for line in tax:
        bin = line.split("\t")[0]
        taxonomy = line.split("\t")[1]
        taxonomy = taxonomy.replace(" ", "_") # removes spaces from species names.
        bin2tax[bin] = taxonomy


### rename  contig headers, store in a new file:
renamed_db = open(output, 'w')
ctgs_not_found = 0
ctgs_renamed = 0
ctgs = 0
with open(in_fasta) as genome:
    for line in genome:
        if line.startswith(">"):
            ctgs += 1
            contig_full_name = line.split(">")[1].strip()

            if contig_full_name in ctg2bin.keys():
                binID = ctg2bin[contig_full_name]
                tax = bin2tax[binID]
                contigID = contig_full_name.split("-")[0]
                new_header = binID + "|" + contigID + "|" + tax + "\n"
                ctgs_renamed += 1

                if ctgs == 1:  # does not need a new line in the beginning of the file
                    line = ">" + new_header
                else:
                    line = "\n>" + new_header
            else:
                ctgs_not_found += 1
                line = "\n>" + contig_full_name + "\n"

        else: # convert interleaved to sequential
            line = line.replace("\n", "")

        renamed_db.write(line)

renamed_db.close()

print ("Done! %i contigs renamed, and %i contigs were not found" %(ctgs_renamed,ctgs_not_found))
