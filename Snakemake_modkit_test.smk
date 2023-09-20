configfile: "config/config.yaml"
configfile: "config/samples_sorted_bams.yaml"


rule all:
    input: 
        expand("results/modkit/{samples}/{samples}_modBase_summary.log", samples=config["samples"]),
        expand("results/modkit/{samples}/{samples}_modBase_summary.txt", samples=config["samples"])    

rule modkit:
    input:
        lambda wildcards: wildcards.samples,
        sample_name=lambda wildcards: wildcards.samples
    output:
         summary_log="results/{rule}/{samples}/{samples}_modBase_summary.log",
         summary_txt="results/{rule}/{samples}/{samples}_modBase_summary.txt"
        #  sample_prob_log="results/{rule}/{samples}/{samples}_modBase_sample_prob.log",
        #  sample_prob_tsv="results/{rule}/{samples}/{samples}_probabilities.tsv",
        #  sample_prob_txt="results/{rule}/{samples}/{samples}_probabilities.txt",
        #  sample_prob_thresh_tsv="results/{rule}/{samples}/{samples}_thresholds.tsv",
        #  modpileup_combined="results/{rule}/{samples}/{samples}_modpileup_combined.bed",
        #  modpileup_5mC="results/{rule}/{samples}/{samples}_modpileup_5mC.bed",
        #  modpileup_5mC_log="results/{rule}/{samples}/{samples}_modpileup_5mC.log",
        #  modpileup_combined_log="results/{rule}/{samples}/{samples}_modpileup_combined.log"
    params:
          reference_genome=config["ref_genome"],
          modkit_threads=32,
          modkit_prob_percentiles=0.1,
          outdir= "results/{rule}/{samples}/"
    log:
        "logs/{rule}/{samples}/{samples}.log"
    shell:
        """ 
module load modkit
echo "{input.sample_name}"
echo ""modkit pileup...""
modkit summary --threads 12 --only-mapped {input} --log-filepath  {output.summary_log} > {output.summary_txt}

# echo ""modkit pileup...""
# # its either C or 5mC(its either unmethylated or 5mC) but what happens to 5hmC when their prob is redistributed? 
# modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_5mC} --cpg --ref {params.reference_genome} --ignore h --log-filepath {output.modpileup_5mC_log}

# # its either C or Methylated(i don't care whether it's 5mC or 5hmC)
# modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_combined} --cpg --ref {params.reference_genome} --combine-mods --log-filepath {output.modpileup_combined_log}


# echo ""find the prob filtering threshold...""
# modkit sample-probs --threads {params.modkit_threads} \
# --percentiles {params.modkit_prob_percentiles} \
# --out-dir {params.outdir} \
# --prefix {input.sample_name} --hist \
# --log-filepath {output.sample_prob_log} \
# {input}
        """
