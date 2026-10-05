#!/bin/bash
#SBATCH -o /dev/null
#SBATCH -e /dev/null


shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"

proj=~/lonchura
vcf=${proj}/vcf/lonchura.snp.vcf.gz

workdir=${proj}/selection
[ ! -e ${workdir} ] && mkdir -p ${workdir}
cd ${workdir}

bcftools query -l ${vcf} | grep -E 'WRM' > wrm.txt
bcftools query -l ${vcf} | grep -v -E 'SBM|WRM' > bf.txt

[ ! -e ${workdir}/daf ] && mkdir ${workdir}/daf

bcftools view --regions NC_042570.1:19400000-20000000 -Oz ${vcf} > ${workdir}/daf/chr4a.vcf.gz
bcftools index ${workdir}/daf/chr4a.vcf.gz

bcftools view --regions NC_042574.1:30800000- -Oz ${vcf} > ${workdir}/daf/chr8.vcf.gz
bcftools index ${workdir}/daf/chr8.vcf.gz

vcfutil daf --vcf ${workdir}/daf/chr4a.vcf.gz --site 0 --pop1 bf.txt --pop2 wrm.txt > ${workdir}/daf/bf_vs_wrm.chr4a.daf
vcfutil daf --vcf ${workdir}/daf/chr8.vcf.gz --site 0 --pop1 bf.txt --pop2 wrm.txt > ${workdir}/daf/bf_vs_wrm.chr8.daf
