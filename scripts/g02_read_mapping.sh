#!/bin/bash
#SBATCH -a 1-59
#SBATCH --mem 16G
#SBATCH -o /dev/null


samples=($(ls clean_data))
samples=($(ls ~/raw_data/songbird2 | sort -V))
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

shopt -s expand_aliases
alias bwa="apptainer exec /usr/local/biotools/b/bwa:0.7.17--h5bf99c6_8 bwa"
alias samtools="apptainer exec /usr/local/biotools/s/samtools:1.18--h50ea8bc_1 samtools"

reference=~/ref/lonStrDom2/GCF_005870125.1.fa
proj=~/vocal-learning
cd ${proj}

mkdir -p bam/${sample}

bwa mem -t 4 -M ${reference} \
  clean_data/${sample}/${sample}_clean_1.fq.gz clean_data/${sample}/${sample}_clean_2.fq.gz |
  samtools view -b -F 4 > bam/${sample}/${sample}.bam
samtools collate -Ou bam/${sample}/${sample}.bam | \
  samtools fixmate - - -mu | \
  samtools sort - -u | \
  samtools markdup - -r bam/${sample}/${sample}.cfsm.bam
samtools index bam/${sample}/${sample}.cfsm.bam

samtools flagstat -O tsv bam/${sample}/${sample}.bam > bam/${sample}/${sample}.stat
samtools flagstat -O tsv bam/${sample}/${sample}.cfsm.bam > bam/${sample}/${sample}.cfsm.stat
samtools coverage bam/${sample}/${sample}.cfsm.bam > bam/${sample}/${sample}.cov

rm bam/${sample}/${sample}.bam
