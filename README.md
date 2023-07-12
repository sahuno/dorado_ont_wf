# template for running snakemake workflows
## date - June 6th 2023
## author - samuel ahuno, mskcc

## how to run snakemake workflow
```
### first activate conda
$ conda activate snakemake #we assume there's conda smk environment already installed on your device

### enter the directory of the snakemake
$cd snakemake_template

### test run with 
$snakemake -s Snakefile.smk -np #target filename may change
$snakemake -s Snakefile.smk --cores 12 --forcerun -np #dry run with cores

#run actual pipeline on the cluster
$nohup snakemake -s Snakefile.smk --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 1:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3 &

#alternatievely make a script that run snalmake
$ sh run_snakefile.sh 
$ cat run_snakefile.sh 
#!/bin/bash

# Run snakemake
snakemake --jobname 's.{jobid}.{rulename}' \
	--snakefile Snakefile_agg_stats_ONT.smk \
	--keep-going \
	--reason \
	--printshellcmds \
	--latency-wait 10 \
	--rerun-incomplete \
	--stats snakemake_$(date +"%Y%m%d_%H%M%S").stats \
	-j 500 \
	--cluster-config config/cluster.json \
	--cluster "bsub -q {cluster.queue} -n {cluster.threads} -W {cluster.time} -M{cluster.mem} -R\"span[hosts=1] select[mem>{cluster.mem}] rusage[mem={cluster.mem}]\" {cluster.extra} -o out.txt -e err.txt" 

```


## how to dump results files to single directory
```
#set your rscript such that files are written into a directory; let your rscript accept dir
#in the specific rule; use keyword `directory()` to specificy dir to save files 
rule all:
    input: 
        expand("results/per_read_aggregate/{samples}/chr{chr}/data/{samples}.chr{chr}.data_aggregate_stats.txt", samples=config["samples"], chr=config["chrs"]),
        expand("results/gene_promoters/{samples}/chr{chr}", samples=config["samples"], chr=config["chrs"])


rule gene_promoters:
    input:
        infile="results/per_read_aggregate/{samples}/chr{chr}/data/{samples}.chr{chr}.data_aggregate_stats.txt"
    output:
        directory("results/gene_promoters/{samples}/chr{chr}")
    params:
        promoter_scripts=config["promoter_plotsR"],
        # chromosomes=lambda wildcards: config["chrs"][wildcards.chrs],
        results_dir="results/gene_promoters/{samples}/chr{chr}"
    log:
       "logs/gene_promoters/{samples}/{samples}.chr{chr}.log"
    shell:
       """
Rscript {params.promoter_scripts} --input_file {input.infile} --dir {params.results_dir} --chrom {wildcards.chr} &> {log}
        """



```