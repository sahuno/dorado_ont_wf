configfile: "config/config.yaml"
configfile: "config/samples_basecalling_modkit.yaml"


rule all:
    input: 
        expand("results/mod_bases/{samples}/{samples}_modBaseCalls_sorted_dedup.bam", samples=config["samples"]),
        expand("results/mod_bases/{samples}/{samples}_modBaseCalls_sorted_dedup.bai", samples=config["samples"]),
        expand("results/mod_bases/{samples}/{samples}_seq_summary.txt", samples=config["samples"])


rule mod_bases:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
         mod_calls_bam="results/{rule}/{samples}/{samples}_modBaseCalls.bam",
         mod_calls_sorted_bam="results/{rule}/{samples}/{samples}_modBaseCalls_sorted.bam",
         mod_calls_sorted_bam_bai="results/{rule}/{samples}/{samples}_modBaseCalls_sorted.bam.bai",
         dedup_mod_calls_sorted_bam="results/{rule}/{samples}/{samples}_modBaseCalls_sorted_dedup.bam",
         dedup_mod_calls_sorted_bai="results/{rule}/{samples}/{samples}_modBaseCalls_sorted_dedup.bai",
         seq_summary="results/{rule}/{samples}/{samples}_seq_summary.txt"
    params:
        methyl_context="5mCG_5hmCG",
        basecall_model_file="/lila/data/greenbaum/users/ahunos/refs/dna_r10.4.1_e8.2_400bps_sup@v4.1.0",
        reference_genome=config["reference_genome"],
        samtools_threads=32
    log:
        "logs/{rule}/{samples}/{samples}.log"
    shell:
        """ 
module load samtools dorado gatk

echo ""running dorado...""
dorado basecaller {params.basecall_model_file} \
/data/greenbaum/users/ahunos/TRI_EPI_DIVYA/pod5/D-A-1/ \
--reference {params.reference_genome} \
--modified-bases {params.methyl_context} --verbose > {output.mod_calls_bam}

echo ""sorting..""
samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} {output.mod_calls_bam}

echo ""indexing..""
samtools index -@ {params.samtools_threads} {output.mod_calls_sorted_bam}

echo ""seq summary..""
dorado summary {output.mod_calls_sorted_bam} > {output.seq_summary}

echo ""marking duplictes..""
gatk MarkDuplicates --CREATE_INDEX true \
--INPUT {output.mod_calls_sorted_bam} \
--OUTPUT {output.dedup_mod_calls_sorted_bam} \
--METRICS_FILE marked_dup_metrics.txt

echo ""delete tmp files..""
rm {output.mod_calls_bam} {output.mod_calls_sorted_bam} {output.mod_calls_sorted_bam.bai}
        """

# rule modkit:
#     input:
#         lambda wildcards: config["samples"][wildcards.samples]
#     output:
#     mod_calls_bam="results/{rule}/{samples}/{samples}_modBaseCalls.bam"
#     mod_calls_sorted_bam="results/{rule}/{samples}/{samples}_modBaseCalls_sorted.bam"
#     mod_calls_sorted_bam_bai="results/{rule}/{samples}/{samples}_modBaseCalls_sorted.bam.bai"
#     CpG_5mC_5hmC_merged_bed="results/{rule}/{samples}/{samples}_5mCG_5hmCG_merged.bed"
#     CpG_5mC_bed="results/{rule}/{samples}/{samples}_5mCG.bed"
#     params:
#         dorado_script=config["dorado"],
#         methyl_context="5mCG_5hmCG",
#         basecall_model_file="/lila/data/greenbaum/users/ahunos/refs/dna_r10.4.1_e8.2_400bps_sup@v4.1.0",
#         reference_genome=
#     log:
#         "logs/{rule}/{samples}/{samples}.log"
#     shell:
#         """ 
#         source {params.dorado_script} {input} \
#         {params.methyl_context} \
#         {params.basecall_model_file} {params.reference_genome} \
#         {output.modified_calls_bam} {output.modified_calls_sorted_bam} {output.modified_calls_sorted_bam_bai} \
#         {output.CpG_5mC_5hmC_merged_bed} {output.CpG_5mC_bed} 2> {log}
#         """