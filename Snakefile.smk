configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input: 
        expand("results/rule1/{samples}/{samples}.proccessed.rds", samples=config["samples"])

rule rule1:
    input:
        input_file_from_config=lambda wildcards: config["samples"][wildcards.samples]
    output:
        output_file_user_sets="results/{rule}/{samples}/{samples}.proccessed.rds"
    params:
        rscript_from_config=config["Rscript_config"]
    log:
        "logs/{rule}/{samples}/{samples}.log"
    shell:
        "Rscript {params.rscript_from_config} --sample {wildcards.samples} --input_file {input.input_file_from_config} --output_file {output.output_file_user_sets} 2> {log}"

# this runs ok
# cd snakemake_template
# snakemake -np #test run with
# snakemake -s Snakefile.smk --cores 12 --forcerun -np #dry run with cores
# nohup snakemake -s Snakefile.smk --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 5:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3 &
