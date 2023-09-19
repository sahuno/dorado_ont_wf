#!/bin/bash

# Run snakemake
#argument 1; snakemake file
#argument 2; cluster config file
snakemake --jobname 's.{jobid}.{rulename}' \
	--snakefile $1 \
	--keep-going \
	--reason \
	--printshellcmds \
	--latency-wait 10 \
	--rerun-incomplete \
	--stats snakemake_$(date +"%Y%m%d_%H%M%S").stats \
	-j 5000 \
	--cluster-config $2 \
	--cluster "bsub -q {cluster.queue} -n {cluster.threads} -W {cluster.time} -M{cluster.mem} -R\"span[hosts=1] select[mem>{cluster.mem}] rusage[mem={cluster.mem}]\" {cluster.extra} -o out.txt -e err.txt"