#!/bin/bash
#SBATCH -a 1-12
#SBATCH --mem 32G

shopt -s expand_aliases
alias hisat2="apptainer exec /usr/local/biotools/h/hisat2:2.2.1--h1b792b2_3 hisat2"
alias samtools="apptainer exec /usr/local/biotools/s/samtools:1.18--h50ea8bc_1 samtools"

proj=~/vocal-learning
ref=~/ref/lonStrDom2/lonStrDom2_hs2
clean_data=${proj}/rnaseq_clean_data
samples=($(ls ${clean_data} | grep -e "D$" | sort -V))
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

mkdir -p ${proj}/rnaseq_bam/${sample}
hisat2 -x ${ref} \
  -1 ${clean_data}/${sample}/${sample}_clean_1.fq.gz \
  -2 ${clean_data}/${sample}/${sample}_clean_2.fq.gz \
  -p 4 |\
  samtools sort -O BAM - -o ${proj}/rnaseq_bam/${sample}/${sample}_hs2.bam
samtools index ${proj}/rnaseq_bam/${sample}/${sample}_hs2.bam
