#!/bin/bash
#SBATCH -o /dev/null
#SBATCH -e /dev/null


workdir=~/vocal-learning/snpeff
vcf=~/vocal-learning/vcf/lonchura.vcf.gz
snpvcf=~/vocal-learning/vcf/lonchura.snp.vcf.gz

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias bedtools="apptainer exec /usr/local/biotools/b/bedtools:2.31.0--h468198e_0 bedtools"

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}


# Around HTR2C; NC_042570.1: 19500001..19860000

gene_id=HTR2C
region=NC_042570.1:19500001-19860000
bcftools view --min-ac 1 --regions ${region} -m2 -M2 -Ov ${vcf} | snpEff lonStrDom2 -stats ${gene_id} > ${gene_id}.annot.vcf
rm ${gene_id}


# Around SLC44A5; NC_042574.1: 30840000..30981424

gene_id=SLC44A5
region=NC_042574.1:30840000-30981424
bcftools view --min-ac 1 --regions ${region} -m2 -M2 -Ov ${vcf} | snpEff lonStrDom2 -stats ${gene_id} > ${gene_id}.annot.vcf
rm ${gene_id}


# Z_FST top 0.1%

bed=~/vocal-learning/snpeff/zfst001.bed
bcftools view --min-af 0.05 -m2 -M2 --regions-file ${bed} -Ov ${vcf} | snpEff lonStrDom2 - -stats zfst001 > zfst001.maf005.annot.vcf
rm tmp.vcf zfst001
