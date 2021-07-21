#!/usr/bin/env python

# set configurations
INDEX_SIZE = config.get("index_size", "12G")
MM_MEM = config.get("minimap_mem", 35000)
MM_PPN = config.get("minimap_ppn", 10)
VAMB_MEM = config.get("vamb_mem", 20000)
VAMB_PPN = config.get("vamb_ppn", 12)
SAMPLE_DATA = config.get("sample_data", "samples2data.txt")
CONTIGS = config.get("contigs", "contigs.txt")
VAMB_PARAMS = config.get("vamb_params", "-o C -m 2000 --minfasta 500000")
VAMB_PRELOAD = config.get("vamb_preload", "")

# parse if GPUs is needed #
#VAMB_split = VAMB_PPN.split(":") 
#VAMB_threads = VAMB_split[0]
VAMB_threads = VAMB_PPN

## read in sample information ##

# read in sample2path
IDS = []
sample2path = {}
fh_in = open(SAMPLE_DATA, 'r')
for line in fh_in:
    line = line.rstrip()
    fields = line.split('\t')
    IDS.append(fields[0])
    sample2path[fields[0]] = [fields[1], fields[2]]

# read in list of per-sample assemblies
contigs_list = []
fh_in = open(CONTIGS, 'r')
for line in fh_in:
    line = line.rstrip()
    contigs_list.append(line)


## start of snakemake rules ##

# targets
rule all:
    input:
        "1_jgi_matrix/jgi.abundance.dat",
        "2_vamb/clusters.tsv"

rule cat_contigs:
    input:
        contigs_list
    output:
        "contigs.flt.fna.gz"
    resources:
        time_min=480, mem_mb=10000, cpus=1
    threads:
        int(1)
    log:
        "log/contigs/catcontigs.log"
    shell:
        "concatenate.py {output} {input} -m 2000"

rule index:
    input:
        contigs = "contigs.flt.fna.gz"
    output:
        mmi = "contigs.flt.mmi"
    resources:
        time_min=480, mem_mb=230000, cpus=8
    threads:
        int(1)
    log:
        "log/contigs/index.log"
    benchmark:
        "benchmarks/index/"+"index.benchmark.txt"
    shell:
        "minimap2 -I {INDEX_SIZE} -d {output} {input} 2> {log}"

rule dict:
    input:
        contigs = "contigs.flt.fna.gz"
    output:
        dict = "contigs.flt.dict"
    resources:
        time_min=240, mem_mb=12000, cpus=1
    threads:
        int(1)
    log:
        "log/contigs/dict.log"
    shell:
        "samtools dict {input} | cut -f1-3 > {output} 2> {log}"

rule minimap:
    input:
        fq = lambda wildcards: sample2path[wildcards.sample],
        mmi = "contigs.flt.mmi",
        dict = "contigs.flt.dict"
    output:
        bam = temp("mapped/{sample}.bam")
    resources:
        time_min=360, mem_mb=MM_MEM, cpus=VAMB_PPN
    threads:
        int(MM_PPN)
    log:
        "log/map/{sample}.minimap.log"
    benchmark:
        "benchmarks/minimap/"+"{sample}.minimap.benchmark.txt"
    conda:
        "envs/minimap2.yaml"
    shell:
        '''minimap2 -t {threads} -ax sr {input.mmi} {input.fq} | grep -v "^@" | cat {input.dict} - | samtools view -F 3584 -b - > {output.bam} 2>{log}'''

rule sort:
    input:
        "mapped/{sample}.bam"
    output:
        temp("mapped/{sample}.sort.bam")
    resources:
        time_min=240, mem_mb=100000, cpus=2
    params:
        prefix="mapped/tmp.{sample}"
    threads:
        int(2)
    log:
        "log/map/{sample}.sort.log"
    benchmark:
        "benchmarks/sort/"+"{sample}.sort.benchmark.txt"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools sort {input} -T {params.prefix} --threads 1 -m 3G -o {output} 2>{log}"

rule jgi:
    input:
        bam = "mapped/{sample}.sort.bam"
    output:
        jgi = temp("1_jgi/{sample}.raw.jgi")
    resources:
        time_min=60, mem_mb=10000, cpus=1
    threads:
        int(1)
    log:
        "log/jgi/{sample}.jgi"
    benchmark:
        "benchmarks/jgi/"+"{sample}.jgi.benchmark.txt"
    shell:
        """
        echo 'running jgi...'
        jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth {output} {input} 2>{log}"""

rule cut_column1to3: 
    input:
        "1_jgi/%s.raw.jgi" % IDS[0] 
    output:
        "1_jgi/jgi.column1to3"
    resources:
        time_min=240, mem_mb=10000, cpus=1
    log:
        "log/jgi/column1to3"
    benchmark:
        "benchmarks/cut_column1to3/"+"cut_column1to3.benchmark.txt"
    shell: 
        "cut -f1-3 {input} > {output} 2>{log}"

rule cut_column4to5:
    input:
        "1_jgi/{sample}.raw.jgi"
    output:
        "1_jgi/{sample}.cut.jgi"
    resources:
        time_min=60, mem_mb=10000, cpus=1
    log:
        "log/jgi/{sample}.cut.log"
    benchmark:
        "benchmarks/cut_column4to5/"+"{sample}.cut_column4to5.benchmark.txt"
    shell: 
        "cut -f1-3 --complement {input} > {output} 2>{log}"

rule paste_abundances:
    input:
        column1to3="1_jgi/jgi.column1to3",
        data=expand("1_jgi/{sample}.cut.jgi", sample=IDS)
    output:
        "1_jgi_matrix/jgi.abundance.dat" 
    resources:
        time_min=60, mem_mb=10000, cpus=1
    log:
        "log/jgi/paste_abundances.log"
    shell: 
        "paste {input.column1to3} {input.data} > {output} 2>{log}" 

rule vamb:
    input:
        jgi = "1_jgi_matrix/jgi.abundance.dat",
        contigs = "contigs.flt.fna.gz"
    output:
        "2_vamb/clusters.tsv",
        "2_vamb/latent.npz",
        "2_vamb/lengths.npz",
        "2_vamb/log.txt",
        "2_vamb/model.pt",
        "2_vamb/mask.npz",
        "2_vamb/tnf.npz"
    resources:
        time_min=7200, cpus=VAMB_PPN, mem_mb=VAMB_MEM
    log:
        "log/vamb/vamb.log"
    benchmark:
        "benchmarks/vamb/"+"vamb.benchmark.txt"
    threads:
        int(VAMB_threads)
    conda:
        "envs/vamb.yaml"
    shell:
        "{VAMB_PRELOAD}"
        "rm -rf 2_vamb;"
        "vamb --outdir 2_vamb --fasta {input.contigs} --jgi {input.jgi} -p {VAMB_threads} {VAMB_PARAMS} 2>{log}"

