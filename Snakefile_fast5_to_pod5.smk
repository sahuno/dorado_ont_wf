configfile: "config/config.yaml"
configfile: "config/samples_fast5_2_pod5_spectrum.yaml"
#$ conda env export -n pod5 > /lila/data/greenbaum/users/ahunos/apps/dorado_ont_wf/config/pod5.yaml


rule all:
    input: 
        expand("results/pod5/{samples}/{samples}.pod5", samples=config["samples"]),
        expand("results/pod5/{samples}/{samples}_pod5_stats.tsv", samples=config["samples"])

rule pod5:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
        out_pod5="results/{rule}/{samples}/{samples}.pod5",
        out_stats="results/{rule}/{samples}/{samples}_pod5_stats.tsv"
    log:
      "logs/{rule}/{samples}/{samples}.log"
    # conda: config["pod5_env"]
    shell:
        """ 
pod5 convert fast5 {input} --threads 64 --output {output.out_pod5} --force-overwrite 2> {log}
pod5 view -t 64 {output.out_pod5} --output {output.out_stats} --force-overwrite
        """

# this runs ok
# snakemake -s Snakefile_fast5_pod5.smk --cores 12 --forcerun -np #dry run with cores
# nohup snakemake -s Snakefile_fast5_pod5.smk --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 5:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3 &
# sh run_snakefile_pod5.sh 

# snakemake -s Snakefile_fast5_pod5.smk -np --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 5:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3
#conda env create -f environment.yml
# /opt/common/CentOS_7/anaconda/anaconda3/bin/conda install -n base -c conda-forge mamba