#!/bin/bash
#SBATCH -a 1-97
#SBATCH --mem 63G
#SBATCH -o /dev/null
#SBATCH -e /dev/null


shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"

reference=~/ref/lonStrDom2/GCF_005870125.1.fa
proj=~/vocal-learning
cd ${proj}

samples=($(ls bam | sort -V))
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

bcftools mpileup -f ${reference} bam/${sample}/${sample}.cfsm.bam --max-depth 500 --no-BAQ | \
  bcftools call -vm -Oz -o bam/${sample}/${sample}.vcf.gz
bcftools index bam/${sample}/${sample}.vcf.gz
bcftools stats bam/${sample}/${sample}.vcf.gz > bam/${sample}/${sample}.vcf.stat

mode_depth=$(cat bam/${sample}/${sample}.vcf.stat | awk -F '\t' '$1 == "DP"' | sort -t $'\t' -k7,7nr | awk -F '\t' 'NR==1 {print $3}')
max_depth=$((${mode_depth}*3))
min_depth=$((${mode_depth}/3))

bcftools filter \
  -i "QUAL>30 && INFO/DP>=${min_depth} && MQ>20" \
  --set-GTs . \
  bam/${sample}/${sample}.vcf.gz \
  -Oz -o bam/${sample}/${sample}.q.vcf.gz
bcftools index bam/${sample}/${sample}.q.vcf.gz
bcftools stats bam/${sample}/${sample}.q.vcf.gz > bam/${sample}/${sample}.q.vcf.stat
