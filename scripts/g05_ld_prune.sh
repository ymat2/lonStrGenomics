#!/bin/bash
#SBATCH -o /dev/null
#SBATCH -e /dev/null

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias plink2="apptainer exec /usr/local/biotools/p/plink2:2.00a5--h4ac6f70_0 plink2"

proj=~/vocal-learning
workdir=${proj}/structure
vcf=${proj}/vcf/lonchura.snp.vcf.gz
prefix=lonchura

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

plink2 --vcf ${vcf} \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --indep-pairwise 50 10 0.1 \
  --out ${prefix}

plink2 --vcf ${vcf} \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --extract ${prefix}.prune.in \
  --make-bed \
  --out ${prefix}.snp

sed -i -e 's/NW_//g' ${prefix}.snp.bim
sed -i -e 's/NC_//g' ${prefix}.snp.bim
sed -i -e 's/\.1//g' ${prefix}.snp.bim

rm ${prefix}.prune*
