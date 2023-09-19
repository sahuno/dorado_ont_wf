configfile: "config/config.yaml"
configfile: "config/samples.yaml"


rule all:
    input: 
        expand("results/mod_bases/{samples}/{samples}_modBaseCalls_sorted_dedup.bam", samples=config["samples"]),
        expand("results/mod_bases/{samples}/{samples}_modBaseCalls_sorted_dedup.bam.bai", samples=config["samples"]),
        expand("results/mod_bases/{samples}/{samples}_seq_summary.txt", samples=config["samples"])


rule mod_bases:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
         mod_calls_bam="results/{rule}/{samples}/{samples}_modBaseCalls.bam",
         mod_calls_sorted_bam="results/{rule}/{samples}/{samples}_modBaseCalls_sorted.bam",
         mod_calls_sorted_bam_bai="results/{rule}/{samples}/{samples}_modBaseCalls_sorted.bam.bai",
         dedup_mod_calls_sorted_bam="results/{rule}/{samples}/{samples}_modBaseCalls_sorted_dedup.bam",
         dedup_mod_calls_sorted_bam_bai="results/{rule}/{samples}/{samples}_modBaseCalls_sorted_dedup.bam.bai",
         seq_summary="results/{rule}/{samples}/{samples}_seq_summary.txt"
    params:
        methyl_context="5mCG_5hmCG",
        basecall_model_file="/lila/data/greenbaum/users/ahunos/refs/dna_r10.4.1_e8.2_400bps_sup@v4.1.0",
        reference_genome=config["ref_genome"],
        samtools_threads=32
    log:
        "logs/{rule}/{samples}/{samples}.log"
    shell:
        """ 
module load samtools dorado gatk

echo ""running dorado...""
dorado basecaller {params.basecall_model_file} \
{input} \
--reference {params.reference_genome} \
--modified-bases {params.methyl_context} --verbose > {output.mod_calls_bam}

echo ""sorting..""
samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} {output.mod_calls_bam}

echo ""indexing..""
samtools index -@ {params.samtools_threads} {output.mod_calls_sorted_bam}

echo ""seq summary..""
dorado summary {output.mod_calls_sorted_bam} > {output.seq_summary}

echo ""marking duplictes..""
gatk MarkDuplicates \
--INPUT {output.mod_calls_sorted_bam} \
--OUTPUT {output.dedup_mod_calls_sorted_bam} \
--METRICS_FILE marked_dup_metrics.txt

echo ""indexing..""
samtools index -@ {params.samtools_threads} {output.dedup_mod_calls_sorted_bam}

echo ""delete tmp files..""
rm {output.mod_calls_bam} {output.mod_calls_sorted_bam} {{output.mod_calls_sorted_bam.bai}}
        """

rule modkit:
    input:
        "results/{rule}/{samples}/{samples}_modBaseCalls_sorted_dedup.bam",
          sample_name=lambda wildcards: wildcards.samples
    output:
         summary_log="results/{rule}/{samples}/{samples}_modBase_summary.log",
         summary_txt="results/{rule}/{samples}/{samples}_modBase_summary.txt",
         sample_prob_log="results/{rule}/{samples}/{samples}_modBase_sample_prob.log",
         sample_prob_tsv="results/{rule}/{samples}/{samples}_probabilities.tsv",
         sample_prob_txt="results/{rule}/{samples}/{samples}_probabilities.txt",
         sample_prob_thresh_tsv="results/{rule}/{samples}/{samples}_thresholds.tsv",
         modpileup_combined="results/{rule}/{samples}/{samples}_modpileup_combined.bed",
         modpileup_5mC="results/{rule}/{samples}/{samples}_modpileup_5mC.bed",
         modpileup_5mC_log="results/{rule}/{samples}/{samples}_modpileup_5mC.log",
         modpileup_combined_log="results/{rule}/{samples}/{samples}_modpileup_combined.log"

    params:
          reference_genome=config["ref_genome"],
          modkit_threads=32,
          modkit_prob_percentiles=0.1,
          outdir= "results/{rule}/{samples}/"
    log:
        "logs/{rule}/{samples}/{samples}.log"
    shell:
        """ 
mododule load modkit

echo ""modkit pileup...""
modkit summary --threads 12 --only-mapped {input} --log-filepath  {output.summary_log} > {output.summary_txt}

echo ""modkit pileup...""
# its either C or 5mC(its either unmethylated or 5mC) but what happens to 5hmC when their prob is redistributed? 
modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_5mC} --cpg --ref {params.reference_genome} --ignore h --log-filepath {output.modpileup_5mC_log}

# its either C or Methylated(i don't care whether it's 5mC or 5hmC)
modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_combined} --cpg --ref {params.reference_genome} --combine-mods --log-filepath {output.modpileup_combined_log}


echo ""find the prob filtering threshold...""
modkit sample-probs --threads {params.modkit_threads} \
--percentiles {params.modkit_prob_percentiles} \
--out-dir {params.outdir} \
--prefix {input.sample_name} --hist \
--log-filepath {output.sample_prob_log} \
{input}
        """


