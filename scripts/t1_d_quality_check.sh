#!/bin/bash
#SBATCH -a 1-12
#SBATCH --mem 16G


proj=~/vocal-learning
cd ${proj}
raw_data=${proj}/raw_data/rnaseq

samples=($(ls ${raw_data} | grep -e "D$" | sort -V))
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

shopt -s expand_aliases
alias fastp="apptainer exec /usr/local/biotools/f/fastp:0.23.4--h5f740d0_0 fastp"

[ ! -e uncleaned/${sample} ] && mkdir -p uncleaned/${sample}
[ ! -e rnaseq_clean_data/${sample} ] && mkdir -p rnaseq_clean_data/${sample}
[ ! -e out/qc_rna ] && mkdir -p out/qc_rna

paired_1=($(ls ${raw_data}/${sample}/*1.fastq.gz))
paired_2=($(ls ${raw_data}/${sample}/*2.fastq.gz))

fastp \
  -i ${paired_1} \
  -I ${paired_2} \
  -o rnaseq_clean_data/${sample}/${sample}_clean_1.fq.gz \
  -O rnaseq_clean_data/${sample}/${sample}_clean_2.fq.gz \
  -h /dev/null \
  -j out/qc_rna/${sample}_qc.json \
  --disable_quality_filtering
