#!/bin/bash
#SBATCH -a 1-12
#SBATCH --mem 16G

shopt -s expand_aliases
alias stringtie="apptainer exec /usr/local/biotools/s/stringtie:2.2.3--h43eeafb_0 stringtie"

proj=~/vocal-learning
ref_gff=~/ref/lonStrDom2/GCF_005870125.1.gff
bamdir=${proj}/rnaseq_bam
samples=($(ls ${bamdir} | grep -e "D$" | sort -V))
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

stringtie -e -p 4 \
  -G ${ref_gff} \
  -o ${bamdir}/${sample}/${sample}_stringtie.gtf \
  ${bamdir}/${sample}/${sample}_hs2.bam
