#author: Samuel Ahuno
#purpose: run pycoQC on sequencing_summary.txt files from `guppy`. sequencing metrics was generated by megaloodon pipeline  

import sys
import os
import glob
import subprocess

#how to run 
## $ conda activate pod5


#set paths to meagalodon data
#dir_parent = sys.argv[1] #get path from command line
dir_parent = '/data/greenbaum/projects/ONT/projects/RERUN/TRI_EPIGENETIC/organized/FAST5/'
dest_dir =  "/data/greenbaum/users/ahunos/TRI_EPI_DIVYA" #output project dir name
#dest_dir = sys.argv[2] #get path from command line


#get stem names of fast5 files
fast5_dirs_list = os.listdir(dir_parent)
sample_list = []


#test if dir has fast5 files
for s in fast5_dirs_list:
    for fname in os.listdir(f'{dir_parent}{s}/'):
        if fname.endswith('.fast5'):
            # do stuff on the file
            os.system(f'mkdir -p {dest_dir}/pod5/{s}')
            os.system(f'mkdir -p {dest_dir}/mod_bases/{s}')
            os.system(f'mkdir -p {dest_dir}/scripts')
            os.system(f'mkdir -p {dest_dir}/logs')
            os.system(f'touch {dest_dir}/scripts/{s}.job') #create job submission script
            sample_list.append(s)
            break
    else:
        # do stuff if a file .true doesn't exist.
        print(f'directory {s} doesnt contain .fast5 files')



#dorado options
methyl_context="5mCG_5hmCG"
basecall_model_file="/lila/data/greenbaum/users/ahunos/refs/dna_r10.4.1_e8.2_400bps_sup@v4.1.0"
reference_genome="/data/greenbaum/database/mm10/mm10.fa"


####################
#create lsf job submission scripts to run fast5 to pod5
#################
sample_small=list(filter(lambda k: '-0-' in k, sample_list))
for sample in sample_small:
    job_file = os.path.join(f'{dest_dir}/scripts/{sample}.job')
    output = os.path.join(f'{dest_dir}/pod5/', sample)
    sample_data = os.path.join(f'{dir_parent}', sample)
    # Create sample directories
    with open(job_file, 'w') as fh:  # Added 'w' to open the file for writing
        fh.writelines("#!/bin/bash\n")
        fh.writelines(f'#BSUB -J pod5_{sample}\n')
        fh.writelines("#BSUB -q gpuqueue\n")
        fh.writelines("#BSUB -gpu num=4\n")  # Corrected the quotation marks
        fh.writelines("#BSUB -n 64\n")
        fh.writelines("#BSUB -R rusage[mem=4]\n")  # Corrected the quotation marks
        fh.writelines("#BSUB -W 5:00\n")
        fh.writelines(f'#BSUB -e {dest_dir}/logs/pod5_{sample}_%J_%I.err\n')
        fh.writelines(f'#BSUB -o {dest_dir}/logs/pod5_{sample}_%J_%I.out\n')
        fh.writelines(f'\n')
        fh.writelines(f'source /data/greenbaum/software/miniconda/etc/profile.d/conda.sh\n')
        fh.writelines(f'conda activate pod5\n')
        # fh.writelines(f'module load dorado\n')
        # fh.writelines(f'module load modkit\n')
        # fh.writelines(f'module load samtools\n')
        fh.writelines(f'pod5 convert from_fast5 {sample_data} --output {output}/{sample}.pod5 --force-overwrite\n')
        fh.writelines(f'\n')
        fh.writelines(f'echo ""running dorado...""\n')
        # fh.writelines(f'dorado basecaller {basecall_model_file} {output}/{sample}.pod5 --reference {reference_genome} --modified-bases {methyl_context} --verbose | samtools view --threads 32 -O BAM -o {dest_dir}/mod_bases/{sample}/{sample}_5mCG_5hmCG_calls.bam\n')
        # fh.writelines(f'echo ""sorting..""\n')
        # fh.writelines(f'samtools sort --threads 32 -o {dest_dir}/mod_bases/{sample}/{sample}_5mCG_5hmCG_calls_sorted.bam {dest_dir}/mod_bases/{sample}/{sample}_5mCG_5hmCG_calls.bam\n')
        # fh.writelines(f'echo ""indexing..""\n')
        # fh.writelines(f'samtools index {dest_dir}/mod_bases/{sample}/{sample}_5mCG_5hmCG_calls_sorted.bam\n')
        # fh.writelines(f'\n')
        # fh.writelines(f'\n')
        # fh.writelines(f'echo ""modkit pileup...""\n')
        # fh.writelines(f'modkit pileup {dest_dir}/mod_bases/{sample}/{sample}_5mCG_5hmCG_calls_sorted.bam {dest_dir}/mod_bases/{sample}/{sample}_CpG_5mC_bed --cpg --ref {reference_genome} --ignore h\n')
        # fh.writelines(f'modkit pileup {dest_dir}/mod_bases/{sample}/{sample}_5mCG_5hmCG_calls_sorted.bam {dest_dir}/mod_bases/{sample}/{sample}_CpG_5mC_5hmC_merged_bed --cpg --ref {reference_genome} --combine-mods\n')
    os.system(f'bsub < {job_file}')
#os.system(f'conda activate pod5')
#l /data/greenbaum/users/ahunos/TRI_EPI_DIVYA/logs/pod5_D-C-319626774_0.err conda init bash
#fh.writelines(f'echo ""pod5 convert from_fast5 {sample_data} --output {output}/{sample}.pod5 --force-overwrite\n""')
        #fh.writelines(f'pod5 -h\n')



    

#sort bams and index



# ## methylation bed generation



sample_small=list(filter(lambda k: '-A-' in k, sample_list))
for sample in sample_small:
    os.system(f'touch {dest_dir}/scripts/{sample}_dorado.job') #create job submission script
    dorado_job_file = os.path.join(f'{dest_dir}/scripts/{sample}_dorado.job')
    output = os.path.join(f'{dest_dir}/mod_bases/{sample}')
    sample_data = os.path.join(f'{dest_dir}/pod5/{sample}')
    # Create sample directories
    with open(dorado_job_file, 'w') as fh:  # Added 'w' to open the file for writing
        fh.writelines("#!/bin/bash\n")
        fh.writelines(f'#BSUB -J pod5_{sample}\n')
        fh.writelines("#BSUB -q gpuqueue\n")
        fh.writelines("#BSUB -gpu num=2\n")  # Corrected the quotation marks
        fh.writelines("#BSUB -n 32\n")
        fh.writelines("#BSUB -R rusage[mem=32]\n")  # Corrected the quotation marks
        fh.writelines("#BSUB -W 10:00\n")
        fh.writelines(f'#BSUB -e {dest_dir}/logs/dorado_{sample}_%J_%I.err\n')
        fh.writelines(f'#BSUB -o {dest_dir}/logs/dorado_{sample}_%J_%I.out\n')
        fh.writelines(f'\n')
        fh.writelines(f'source /data/greenbaum/software/miniconda/etc/profile.d/conda.sh\n')
        fh.writelines(f'module load dorado\n')
        fh.writelines(f'module load modkit\n')
        fh.writelines(f'module load samtools\n')
        fh.writelines(f'\n')
        fh.writelines(f'echo ""running dorado...""\n')
        fh.writelines(f'dorado basecaller {basecall_model_file} {sample_data}/ --reference {reference_genome} --modified-bases {methyl_context} --verbose | samtools view --threads 32 -O BAM -o {output}/{sample}_5mCG_5hmCG_calls.bam\n')
        fh.writelines(f'echo ""sorting..""\n')
        fh.writelines(f'samtools sort --threads 32 -o {output}/{sample}_5mCG_5hmCG_calls_sorted.bam {output}/{sample}_5mCG_5hmCG_calls.bam\n')
        fh.writelines(f'echo ""indexing..""\n')
        fh.writelines(f'samtools index -@ 32 {output}/{sample}_5mCG_5hmCG_calls_sorted.bam\n')
        fh.writelines(f'\n')
        fh.writelines(f'\n')
        fh.writelines(f'echo ""modkit pileup...""\n')
        fh.writelines(f'modkit pileup {output}/{sample}_5mCG_5hmCG_calls_sorted.bam {output}/{sample}_CpG_5mC_bed --cpg --ref {reference_genome} --ignore h\n')
        fh.writelines(f'modkit pileup {output}/{sample}_5mCG_5hmCG_calls_sorted.bam {output}/{sample}_CpG_5mC_5hmC_merged_bed --cpg --ref {reference_genome} --combine-mods\n')
    os.system(f'bsub < {dorado_job_file}')