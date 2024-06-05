parent_dir = "/home/ahunos/apps/dorado_ont_wf/"
configfile: parent_dir + "config/config.yaml"
configfile: "/home/ahunos/apps/dorado_ont_wf/config/samples_mergeBam_modkit.yaml"
set_species = "human"

rule all:
    input: 
        expand("results/merge_bams/{samples}/{samples}_modBaseCalls_dedup_sorted.bam", samples=config["samples"]),
        expand("results/merge_bams/{samples}/{samples}_modBaseCalls_dedup_sorted.bam.bai", samples=config["samples"]),
        expand("results/mark_duplicates/{samples}/{samples}_modBaseCalls_sorted_dup.bam", samples=config["samples"]),
        expand("results/mark_duplicates/{samples}/{samples}_modBaseCalls_sorted_dup.bai", samples=config["samples"]),
        expand("results/mark_duplicates/{samples}/{samples}_marked_dup_metrics.txt", samples=config["samples"]),
        expand("results/modkit/{samples}/{samples}_modBase_summary.log", samples=config["samples"]),
        expand("results/modkit/{samples}/{samples}_modBase_summary.txt", samples=config["samples"]),
        expand("results/modkit/{samples}/{samples}_modpileup_combined.bed", samples=config["samples"]),    
        expand("results/modkit/{samples}/{samples}_modpileup_5mC.bed", samples=config["samples"]),
        expand("results/modkit/{samples}/{samples}_modpileup_5mC.log", samples=config["samples"]),   
        expand("results/modkit/{samples}/{samples}_modpileup_combined.log", samples=config["samples"]),
        expand("results/modkit/{samples}/{samples}_modBase_sample_prob.log", samples=config["samples"]), 
        expand("results/modkit/{samples}/{samples}_probabilities.tsv", samples=config["samples"]), 
        expand("results/modkit/{samples}/{samples}_probabilities.txt", samples=config["samples"]), 
        expand("results/modkit/{samples}/{samples}_thresholds.tsv", samples=config["samples"]) 




rule merge_bams:
    input:
        bamsAll=lambda wildcards: config["samples"][wildcards.samples]
        # bam1=lambda wildcards: config["samples"][wildcards.samples][0]
    output:
         merged_bam="results/merge_bams/{samples}/{samples}_modBaseCalls_dedup_sorted.bam",
         merged_bai="results/merge_bams/{samples}/{samples}_modBaseCalls_dedup_sorted.bam.bai"
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        samtools_threads=16
    log:
        "logs/merge_bams/{samples}/{samples}.log"
    benchmark:
        "benchmarks/merge_bams/{samples}/{samples}.benchmark.txt"
    shell:
        """ 
        echo "Merging bams for {wildcards.samples}"
        samtools merge --threads {params.samtools_threads} -o {output.merged_bam} {input.bamsAll} && samtools index -@ 32 {output.merged_bam} 2> {log}
        """
#SILENT
#echo "first bam file {input.bam1} of sample {wildcards.samples}"
rule mark_duplicates:
    input:
        bams="results/merge_bams/{samples}/{samples}_modBaseCalls_dedup_sorted.bam"
    output:
         markdup_bam="results/mark_duplicates/{samples}/{samples}_modBaseCalls_sorted_dup.bam",
         markdup_bam_bai="results/mark_duplicates/{samples}/{samples}_modBaseCalls_sorted_dup.bai",
         markdup_stats="results/mark_duplicates/{samples}/{samples}_marked_dup_metrics.txt"
    log:
        "logs/mark_duplicates/{samples}/{samples}.log"
    benchmark:
        "benchmarks/mark_duplicates/{samples}/{samples}.benchmark.txt"
    shell:
        """ 
        module load java
        gatk MarkDuplicates --INPUT {input.bams} --OUTPUT {output.markdup_bam} --METRICS_FILE {output.markdup_stats} --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT 2> {log}
        """

rule modkit:
    input:
        "results/mark_duplicates/{samples}/{samples}_modBaseCalls_sorted_dup.bam"
    output:
         summary_log="results/modkit/{samples}/{samples}_modBase_summary.log",
         summary_txt="results/modkit/{samples}/{samples}_modBase_summary.txt",
         sample_prob_log="results/modkit/{samples}/{samples}_modBase_sample_prob.log",
         sample_prob_tsv="results/modkit/{samples}/{samples}_probabilities.tsv",
         sample_prob_txt="results/modkit/{samples}/{samples}_probabilities.txt",
         sample_prob_thresh_tsv="results/modkit/{samples}/{samples}_thresholds.tsv",
         modpileup_combined="results/modkit/{samples}/{samples}_modpileup_combined.bed",
         modpileup_5mC="results/modkit/{samples}/{samples}_modpileup_5mC.bed",
         modpileup_5mC_log="results/modkit/{samples}/{samples}_modpileup_5mC.log",
         modpileup_combined_log="results/modkit/{samples}/{samples}_modpileup_combined.log"
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        modkit_threads=4,
        modkit_prob_percentiles=0.1,
        outdir= "results/modkit/{samples}/"
    log:
        "logs/modkit/{samples}/{samples}.log"
    benchmark:
        "benchmarks/modkit/{samples}/{samples}.benchmark.txt"
    conda: config["pod5_env"]
    shell:
        """ 
################################################        
echo ""modkit pileup...""
modkit summary --threads 12 --only-mapped {input} --log-filepath  {output.summary_log} > {output.summary_txt} 2> {log}

################################################
echo ""modkit pileup...""
# its either C or 5mC(its either unmethylated or 5mC) but what happens to 5hmC when their prob is redistributed? 
modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_5mC} --cpg --ref {params.reference_genome} --ignore h --log-filepath {output.modpileup_5mC_log} --only-tabs 2> {log}
echo ""done modkit pileup...""

################################################
echo ""modkit combine 5mc 5hmc pileup...""
# its either C or Methylated(i don't care whether it's 5mC or 5hmC)
modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_combined} --cpg --ref {params.reference_genome} --combine-mods --log-filepath {output.modpileup_combined_log} --only-tabs 2> {log}
echo ""done modkit combine 5mc 5hmc pileup...""

################################################
echo ""find the prob filtering threshold...""
modkit sample-probs --threads {params.modkit_threads} \
--percentiles {params.modkit_prob_percentiles} \
--out-dir {params.outdir} \
--prefix {wildcards.samples} --hist \
--log-filepath {output.sample_prob_log} \
{input} 2> {log}

###########5mc & 5hmc, bedgraph outputs ############
modkit pileup --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
--prefix {wildcards.samples} --cpg --combine-mods --ref {params.reference_genome}

#these files will be generated
# prefix_C_CG0_positive.bedgraph
# prefix_C_CG0_negative.bedgraph

###########5mc & 5hmc, bedgraph outputs############
modkit pileup --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
--prefix {wildcards.samples} --cpg --ref {params.reference_genome}

#these files will be generated
# prefix_h_CG0_negative.bedgraph
# prefix_m_CG0_negative.bedgraph
# prefix_h_CG0_positive.bedgraph
# prefix_m_CG0_positive.bedgraph

echo ""done modkit sample-probs...""
        """
# snakemake -s /home/ahunos/apps/dorado_ont_wf/mergeBams_mdup_modkit.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs 10 --cores all --use-conda --keep-going --forceall -np