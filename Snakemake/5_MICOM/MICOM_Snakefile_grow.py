# build community models and obtain metabolic exchanges

configfile: "MICOM_config.yaml"
samples_fp = "0_MAGs_tables/all_samples.txt"


# Process samples file
#samples_list = open(config["path"]["root"]+config["samples"]["sam"]) # not working in PyCharm
with open(samples_fp) as f:
    samples = f.read().splitlines()


rule all:
    input:
        expand(config["path"]["root"]+"/"+config["folder"]["tradeoffs"]+"/minimal_fluxes_exchange_{sample}.csv", sample = samples)


rule build_community:
    input:
        sample_file = config["path"]["root"]+"/"+config["folder"]["tables_fp"]+"/{sample}.csv"
    params:
        smpl = "{sample}"
    output:
        config["path"]["root"]+"/"+config["folder"]["pickles"]+"/{sample}.pickle"
    resources:
        time_min=45, mem_mb=10000, cpus=config["cores"]["build_comm"]
    log:
        std_out = config["path"]["root"]+"/"+config["folder"]["logs"]+"/micom_build_comm/{sample}.log"
    benchmark:
        config["path"]["root"]+"/benchmarks/build_comm/"+'{sample}.benchmark.txt'
    shell:
        """     
        echo "Begin building commuities with MICOM ... "
        python3 MICOM_build_comm_models.py -s {params.smpl}
        
        echo "Done"
        """

rule tradeoff:
    input:
        config["path"]["root"]+"/"+config["folder"]["pickles"]+"/{sample}.pickle"
    params:
        smpl = "{sample}.pickle",
        out_folder = config["path"]["root"]+"/"+config["folder"]["tradeoffs"]
    output:
        config["path"]["root"]+"/"+config["folder"]["tradeoffs"]+"/minimal_fluxes_exchange_{sample}.csv"
    resources:
        time_min=360, mem_mb=20000, cpus=config["cores"]["tradeoff"]
    log:
        std_out = config["path"]["root"]+"/"+config["folder"]["logs"]+"/micom_tradeoff_comm/{sample}.log"
    benchmark:
        config["path"]["root"]+"/benchmarks/tradeoffs/"+'{sample}.benchmark.txt'
    shell:
        """
       
        echo "Begin tradeoff analyses with MICOM ... "

        python3 MICOM_coop_tradeoff.py -sn {params.smpl} -t {resources.cpus} -o {params.out_folder}
        
        echo "Done!"
        """

