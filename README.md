# template for running snakemake workflows
## date - June 6th 2023
## author - samuel ahuno, mskcc

## how to run snakemake workflow

### first activate conda
$ conda activate snakemake

### enter the directory of the snakemake
$cd snakemake_template
### test run with 
$snakemake -s Snakefile.smk -np #test run with
$snakemake -s Snakefile.smk --cores 12 --forcerun -np #dry run with cores

#run actual pipeline on the cluster
$nohup snakemake -s Snakefile.smk --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 1:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3 &