#!/bin/bash
#SBATCH -o /dev/null
#SBATCH -e /dev/null


shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias vcftools="apptainer exec /usr/local/biotools/v/vcftools:0.1.16--h9a82719_5 vcftools"

proj=~/vocal-learning
vcf=${proj}/vcf/lonchura.snp.vcf.gz

workdir=${proj}/stats
[ ! -e ${workdir} ] && mkdir -p ${workdir}
cd ${workdir}

bcftools query -l ${vcf} | grep -E 'BF' > bf-ambf.txt
bcftools query -l ${vcf} | grep -E 'TB' > tb.txt
bcftools query -l ${vcf} | grep -E 'TO' > to.txt
bcftools query -l ${vcf} | grep -E 'WP' > wp.txt

## Heterozygosity
vcftools --gzvcf ${vcf} \
  --het \
  --maf 0.01 \
  --max-missing 0.1 \
  --chr NC_042565.1 \
  --chr NC_042566.1 \
  --chr NC_042567.1 \
  --chr NC_042568.1 \
  --chr NC_042569.1 \
  --chr NC_042570.1 \
  --chr NC_042571.1 \
  --chr NC_042572.1 \
  --chr NC_042573.1 \
  --chr NC_042574.1 \
  --chr NC_042575.1 \
  --chr NC_042576.1 \
  --chr NC_042577.1 \
  --chr NC_042578.1 \
  --chr NC_042579.1 \
  --chr NC_042580.1 \
  --chr NC_042581.1 \
  --chr NC_042582.1 \
  --chr NC_042583.1 \
  --chr NC_042584.1 \
  --chr NC_042585.1 \
  --chr NC_042586.1 \
  --chr NC_042587.1 \
  --chr NC_042588.1 \
  --chr NC_042589.1 \
  --chr NC_042590.1 \
  --chr NC_042591.1 \
  --chr NC_042592.1 \
  --chr NC_042593.1 \
  --chr NC_042594.1 \
  --out hetero

## Tajima's D (for subpopulations; wrm and bf are conducted in g07)
vcftools --gzvcf ${vcf} --keep bf-ambf.txt --TajimaD 10000 --out bf-ambf
vcftools --gzvcf ${vcf} --keep tb.txt --TajimaD 10000 --out tb
vcftools --gzvcf ${vcf} --keep to.txt --TajimaD 10000 --out to
vcftools --gzvcf ${vcf} --keep wp.txt --TajimaD 10000 --out wp
