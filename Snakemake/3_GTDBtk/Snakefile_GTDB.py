configfile: "config_MegaHit.yaml"


#create a list of samples
samples = glob_wildcards('0_Quality_Control/03_Rarefied/{sample}_R1.fq.gz').sample


rule all:
    input:
        expand(config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{sample}/contigs.fasta.gz", sample = samples)


rule megahit:
    input:
        R1 = "0_Quality_Control/03_Rarefied/{sample}_R1.fq.gz",
        R2 = "0_Quality_Control/03_Rarefied/{sample}_R2.fq.gz"
    output:
       config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{sample}/contigs.fasta.gz"
    resources:
        time_min=1440, mem_mb=100000, cpus=config["cores"]["megahit"]
    benchmark:
        config["path"]["root"]+"/benchmarks/megahit/"+'{sample}.megahit.benchmark.txt'
    params:
        out_folder = config["path"]["root"]+"/"+config["folder"]["assemblies"]+"/{sample}"
    log:
        std_out = config["path"]["root"]+"/"+config["folder"]["logs"]+"/megahit/{sample}.log"
    shell:
        """
        echo "Running megahit ... "
        megahit -t {config[cores][megahit]} \
            --presets {config[params][assemblyPreset]} \
            --verbose \
            --min-contig-len {config[params][assemblyMin]} \
            -1 {input.R1} \
            -2 {input.R2} \
            -o {params.out_folder}_tmp > {log.std_out} 2>&1;
        echo "done with megahit... "

        echo "Renaming assembly ... "
        mv {params.out_folder}_tmp/final.contigs.fa {params.out_folder}_tmp/contigs.fasta
        
        echo "Fixing contig header names: replacing spaces with hyphens ... "
        sed -i 's/ /-/g' {params.out_folder}_tmp/contigs.fasta

        echo "Zipping and moving assembly ... "
        gzip {params.out_folder}_tmp/contigs.fasta
        mkdir -p {params.out_folder}
        mv {params.out_folder}_tmp/contigs.fasta.gz {output}
        rm -r {params.out_folder}_tmp
        echo "done. "
        """
