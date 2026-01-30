#!/bin/bash
#SBATCH -o /dev/null
#SBATCH -e /dev/null


shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"

proj=~/vocal-learning
cd ${proj}

poppy summary -i bam -o out/flag_summary.tsv --mode flag
poppy summary -i bam -o out/mapping_summary.tsv --mode bam
poppy summary -i bam -o out/vcf_summary.tsv --mode vcf --suffix q.vcf.stat
python3 src/identify_sex.py -i bam -c	NC_042595.1 -o out/sex_chrom_depth.tsv

samples=($(ls bam | sort -V | grep -v -e "SBM"))
bam_file_list=()
for s in ${samples[@]}; do bam_file_list+=("bam/$s/$s.q.vcf.gz"); done

workdir=${proj}/vcf
vcf=${workdir}/lonchura.vcf.gz
snpvcf=${workdir}/lonchura.snp.vcf.gz

[ ! -e ${workdir} ] && mkdir ${workdir}

bcftools merge "${bam_file_list[@]}" -0 -Oz > ${vcf}
bcftools index ${vcf}

bcftools view -v snps --min-ac 1 -m2 -M2 -Oz ${vcf} > ${snpvcf}
bcftools index ${snpvcf}
