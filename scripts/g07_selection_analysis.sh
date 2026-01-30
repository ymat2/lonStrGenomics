#!/bin/bash
#SBATCH -o /dev/null
#SBATCH -e /dev/null


shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias vcftools="apptainer exec /usr/local/biotools/v/vcftools:0.1.16--h9a82719_5 vcftools"

proj=~/vocal-learning
vcf=${proj}/vcf/lonchura.snp.vcf.gz

workdir=${proj}/selection
[ ! -e ${workdir} ] && mkdir -p ${workdir}
cd ${workdir}

bcftools query -l ${vcf} | grep -E 'WRM' > wrm.txt
bcftools query -l ${vcf} | grep -v -E 'SBM|WRM' > bf.txt


##### Pi #####

[ ! -e ${workdir}/pi ] && mkdir ${workdir}/pi

vcftools --gzvcf ${vcf} --keep wrm.txt \
  --window-pi 10000 --window-pi-step 5000 --out ${workdir}/pi/wrm

vcftools --gzvcf ${vcf} --keep bf.txt \
  --window-pi 10000 --window-pi-step 5000 --out ${workdir}/pi/bf


##### Tajima's D #####

[ ! -e ${workdir}/tajimasD ] && mkdir ${workdir}/tajimasD

vcftools --gzvcf ${vcf} --keep wrm.txt --TajimaD 10000 --out ${workdir}/tajimasD/wrm
vcftools --gzvcf ${vcf} --keep bf.txt --TajimaD 10000 --out ${workdir}/tajimasD/bf


##### Fst #####

[ ! -e ${workdir}/fst ] && mkdir ${workdir}/fst

vcftools --gzvcf ${vcf} --weir-fst-pop bf.txt --weir-fst-pop wrm.txt \
  --fst-window-size 10000 --fst-window-step 5000 --out ${workdir}/fst/bf_vs_wrm


##### DAF #####

[ ! -e ${workdir}/daf ] && mkdir ${workdir}/daf

bcftools view -m2 -M2 --min-af 0.05 -Oz ${proj}/vcf/lonchura.vcf.gz > ${workdir}/daf/tmp.vcf.gz
bcftools index ${workdir}/daf/tmp.vcf.gz

vcfutil daf --vcf ${workdir}/daf/tmp.vcf.gz --window_size 20 --stats median --pop1 bf.txt --pop2 wrm.txt > ${workdir}/daf/bf_vs_wrm.windowed.daf
vcfutil daf --vcf ${workdir}/daf/tmp.vcf.gz --extract_daf 0.5 --pop1 bf.txt --pop2 wrm.txt > ${workdir}/daf/bf_vs_wrm.daf50.daf

cat ${workdir}/daf/bf_vs_wrm.daf50.daf | awk 'NR>1 { if ($3 > 0.9) { print $1 "\t" $2 }}' > ${workdir}/daf/daf90.txt
bcftools view --regions-file ${workdir}/daf/daf90.txt ${workdir}/daf/tmp.vcf.gz -Ov | snpEff lonStrDom2 -stats daf90 > ${workdir}/daf/lonchura.daf90.annot.vcf
rm ${workdir}/daf/daf90 ${workdir}/daf/tmp.vcf.gz*
