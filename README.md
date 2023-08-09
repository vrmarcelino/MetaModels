# README

This repository contain scripts used (an many tested but not used) in the workflow associated with the manuscript [**Disease-specific loss of microbial cross-feeding interactions in the human gut**](https://www.biorxiv.org/content/10.1101/2023.02.17.528570v1), from quality control to visualization of the final results.


Intermediate files (e.g. metabolic exchanges obtained wth MICOM) and input files for the R scripts can be found in our Zenodo repository (doi:10.5281/zenodo.7582762).

The Zenodo repository also has a more organised structure - with scripts in folder named by the subset of analyses they were written for.


<br />

## Quality control

QC - Read quality with TrimGalore! (Krueger)
This will remove Illumina adapters and low-quality sequences.
Illumina’s Genome Analyser also uses Sanger quality encoding (Phred), so no need to worry. Minimum len = 80bp. and minimum phred = 25.

```
trim_galore --paired --length 80 --quality 25 --cores 4 --output_dir 01_CleanReads --fastqc 00_ENA_downloads/ERRxx.fastq.gz 00_ENA_downloads/ERRxx.fastq.gz
```

Remove human reads (use python script for batch analyses) - Bowtie2 (Langmead 2012)
```
bowtie2 -p $th -x $db -1 01_CleanReads/ERR1190655_1_val_1.fq.gz -2 01_CleanReads/ERR1190655_2_val_2.fq.gz --un-conc  02_Host_Removed/ERRxx > 02_Host_Removed/ERRxx.sam.temp
```

Rarefy sequences
```
seqtk sample -s 8 02_Host_Removed/ERRxx.1 15000000 > 03_Rarefied/ERRxx_R1.fq
```
<br />

## Metagenome assembly and binning

We performed sequence assembly on individual samples using Megahit (Li et al 2015)

```
megahit -t 12 --presets meta-sensitive --verbose --min-contig-len 1000 -1 {input.R1} -2 {input.R2} -o {params.out_folder}_tmp
```

See 'Assembly_Binning/1_megahit/' folder for the Snakemake workflow.

We then performed co-binning of samples in two batches (Sup. table S1) using the [workflow suggested for VAMB](https://github.com/RasmussenLab/vamb).

<br />

VAMB (Nissen et al 2021) was run with default parameters:
```
vamb --outdir 2_vamb --fasta {input.contigs} --jgi {input.jgi} -p 12 -o C -m 2000 --minfasta 500000
```

See 'Assembly_Binning/2_vamb/' folder for the Snakemake workflow.


<br />

## MAGs QC, classification and abundance:

CheckM (Parks 2015) was used to evaluate the completeness and contamination of the bins:
```
checkm lineage_wf -f  + out_file +  -t 8 -x fna  + in_folder +   + out_folder
```
<br />

The High quality bins (>90% completeness and < 5% contamination) were dereplicated with dRep (Olm et al 2017):
```
dRep dereplicate 1_drep_genomes -g 0_all_HQ_bins/*.fna --genomeInfo checkm_all.csv -p 24 -pa 0.95 --SkipSecondary
```
<br />
Taxonomic classification was performed with GTDBtk (Chaumeil et al 2019).
This step was performed in our Nectar instance (instead of our HPC) because HPC servers often have issues with pplacer (they think they need more memory than they actually need).
We placed the metagenome bins into folders containing 1000 bins each (temp_bins_xx), and run GTDBtk:

```
gtdbtk classify_wf --genome_dir 4_Classification/temp_bins_0 --out_dir 4_Classification/gtdb.outdir.0 --cpus 64
```
<br />
To estimate the abundance of the species-level MAGs in ecah sample, we used KMA (Clausen et al 2018).
```
#index database:
kma index -i 0_merged_bins_renamed.fas -o MAGs_db_sparse -NI -Sparse TG

#run KMA:
kma -ipe $in_file_f $in_file_r -o $out_file -t_db $db -t 3 -1t1 -mem_mode -and -apm p -ef -tmp KMA_temp
 
```

<br />
- **Genome-scale metabolic models**

We then built genome-scale metabolic models for each species-level MAG with CarveMe (Machado et al 2018), using their domain classification (Bacteria or Archaea) as parameter for their universe.

```
carve {input.bin} --gapfill western_diet_gut --mediadb {input.media} -u {params.domain} -v -o {output}
```

<br />
- **Community-wide modelling with MICOM (Diener et al 2020)**

First make a txt file with sample names in the 0_MAGs_tables folder (ls -1 > all_samples.txt). Remove the “.csv” from all lines & “all_samples” from first line.

```bash
screen -L -S micom
mkdir -p z_snakemake_logs
conda activate snakemake_cplex

snakemake --snakefile MICOM_Snakefile_grow.py --latency-wait 60 --cluster 'sbatch --output=z_snakemake_logs/%j.out --error=z_snakemake_logs/%j.out -t {resources.time_min} --mem={resources.mem_mb} -c {resources.cpus} -p short,comp' -j 900 -np

```
<br />

## Metabolite Exchange Scoring System for Interdependence (MESSI):

Calculate number of producers and consumers per metabolite

```bash
python3 MetModels_producers_consumers_per_rxn.py -f 2_exchanges -o 3_parsed_exchanges/producers_consumers.csv

```

Then process the output files in R with the scripts in folder ‘MES/Differences_in_MES’

MESSI == (2 x  ((n_produc * n_cons)/(n_produc+n_cons)))


<br />

## Figure 2: MESSI scores in Health and Disease

Panel a: The tree file was generated with GTDBtk de novo workflow (Chaumeil et al 2020), and visualized with iTOL (Letunic & Bork 2021).

Scripts to reproduce panels b and c are in 'MES/Differences_in_MES' folder.

See 'Figure2b_metabolite_imp_healthy.R' script to reproduce panel b

See 'Figure2c_MES_barplots_all.R' script to reproduce the panel c


<br />

## Crohn’s disease (CD)

Calculate flux considering species abundances (done in HPC)

```bash
smux new-session --time=6:00:00 --mem=200G --partition=short,comp

python3 MetModels_summarize_net_produc.py -f 2_exchanges -o 3_parsed_exchanges/net_produc_consump_merged.csv

python3 MetModels_summarize_total_produc_consump.py -f 2_exchanges -kma 1.1_merged_kma_simplified4summarize_production_consumption.csv -op 3_parsed_exchanges/total_production.csv -oc 3_parsed_exchanges/total_consumption.csv

```

Note that net production here is a table with net production / consumption of metabolites by the microbiome
-> these are the exchanges with the media (“_m”), as they indicate the “excess" of metabolites that are released or consumed from the environment. The values are already corrected for species' relative abundance.

Statistics for H2S production and consumption are detailed within the R scripts (scripts in CD+_focus/R_graphs/xxx.R)

<br />

## Network analyses for CD:

From metadata, produce 2 files containing prefixes of samples of CD and Healthy individuals .
It is important to have the same number of samples in healthy and diseased cohort here, so I am using all 38 samples from the healthy cohort that worked and 38 samples from the CD cohort (randomly deleted other samples)

In HPC, produce nodes and edges file:

```bash
./MetModels_create_global_network_from_list.py -if 2_exchanges/ -s He_CD_prefixes_rarefied.txt -m h2s_e -sp wanted_spp_classification.tsv -on 4_nodes_edges_He2017/He_CD_nodes.csv -oe 4_nodes_edges_He2017/He_CD_edges.csv

./MetModels_create_global_network_from_list.py -if 2_exchanges/ -s He_healthy_prefixes_all_that_worked.txt -m h2s_e -sp wanted_spp_classification.tsv -on 4_nodes_edges_He2017/He_healthy_nodes.csv -oe 4_nodes_edges_He2017/He_healthy_edges.csv
```

The output of these scripts can be processed with R (scripts in CD+_focus/R_graphs/xxx.R) to identify a consortium of microbes with promising therapeutic potential.

<br />
If you have any questions please get in touch: vmarcelino-at-unimelb.edu.au
<br />
