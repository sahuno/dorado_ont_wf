# ONT dorado + modkit snakemake workflow
## date - Sept 3rd 2023
## author - samuel ahuno, mskcc

## how to run snakemake workflow
```
### first activate conda
$ conda activate snakemake #we assume there's conda smk environment already installed on your device

### enter the directory of the snakemake
$cd dorado_ont_wf

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
    --use-conda \
	--keep-going \
	--reason \
	--printshellcmds \
	--latency-wait 10 \
	--rerun-incomplete \
	--stats snakemake_$(date +"%Y%m%d_%H%M%S").stats \
	-j 500 \
	--cluster-config config/cluster.json \
	--cluster "bsub -q {cluster.queue} -n {cluster.threads} -W {cluster.time} -M{cluster.mem} -R\"span[hosts=1] select[mem>{cluster.mem}] rusage[mem={cluster.mem}]\" {cluster.extra} -o out.txt -e err.txt" 

## notes on run time
1. pod5 conversion
you dont need lots of memory for pod5 conversion. Max mem=4GB used. Rather increase nthreads=64, nGPUs=4

2. dorado 
