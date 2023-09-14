configfile: "config/config.yaml"
configfile: "config/samples.yaml"
#$ conda env export -n pod5 > /lila/data/greenbaum/users/ahunos/apps/dorado_ont_wf/config/pod5.yaml


rule all:
    input: 
        expand("results/pod5/{samples}/{samples}.pod5", samples=config["samples"])

rule pod5:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
        "results/{rule}/{samples}/{samples}.pod5"
    log:
      "logs/{rule}/{samples}/{samples}.log"
    conda: config["pod5_env"]
    shell:
        """ 
pod5 convert from_fast5 {input} --output {output} --force-overwrite 2> {log}
        """


# this runs ok
# snakemake -s Snakefile.smk --cores 12 --forcerun -np #dry run with cores
# nohup snakemake -s Snakefile.smk --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 5:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3 &
# sh run_snakefile.sh 

#conda env create -f environment.yml
# /opt/common/CentOS_7/anaconda/anaconda3/bin/conda install -n base -c conda-forge mamba