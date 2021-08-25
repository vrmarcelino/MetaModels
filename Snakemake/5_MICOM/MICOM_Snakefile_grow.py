# build community models and obtain metabolic exchanges

configfile: "MICOM_config_grow.yaml"
samples_fp = "0_MAGs_tables/all_samples.txt"


# Process samples file
#samples_list = open(config["path"]["root"]+config["samples"]["sam"]) # not working in PyCharm
with open(samples_fp) as f:
    samples = f.read().splitlines()


rule all:
    input:
        expand(config["path"]["root"]+"/"+config["folder"]["exchanges"]+"/exchanges_grow_{sample}.csv", sample = samples)


rule build_community:
    input:
        sample_file = config["path"]["root"]+"/"+config["folder"]["tables_fp"]+"/{sample}.csv"
    params:
        smpl = "{sample}"
    output:
        config["path"]["root"]+"/"+config["folder"]["pickles"]+"/{sample}.pickle"
    resources:
        time_min=60, mem_mb=8000, cpus=config["cores"]["build_comm"]
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

rule grow_wf:
    input:
        config["path"]["root"]+"/"+config["folder"]["pickles"]+"/{sample}.pickle"
    params:
        smpl = "{sample}.pickle",
        out_folder = config["path"]["root"]+"/"+config["folder"]["exchanges"],
        pickles_fp = config["path"]["root"]+"/"+config["folder"]["pickles"]
    output:
        config["path"]["root"]+"/"+config["folder"]["exchanges"]+"/exchanges_grow_{sample}.csv"
    resources:
        time_min=60, mem_mb=12000, cpus=config["cores"]["exchanges"]
    log:
        std_out = config["path"]["root"]+"/"+config["folder"]["logs"]+"/micom_grow/{sample}.log"
    benchmark:
        config["path"]["root"]+"/benchmarks/exchanges/"+'{sample}.benchmark.txt'
    shell:
        """
       
        echo "Begin grow workflow to calculate metabolic exchanges with MICOM... "
        echo "using parsimonious FBA"

        python3 MICOM_grow_wf.py -c {params.pickles_fp} -s {params.smpl} -th {resources.cpus} -o {params.out_folder}
        
        echo "Done!"
        """
