#!/bin/bash
#SBATCH -o /dev/null
#SBATCH -e /dev/null


shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias vcftools="apptainer exec /usr/local/biotools/v/vcftools:0.1.16--h9a82719_5 vcftools"

proj=~/vocal-learning
vcf=${proj}/vcf/lonchura.snp.vcf.gz

workdir=${proj}/selection/fst_supple
[ ! -e ${workdir} ] && mkdir -p ${workdir}
cd ${workdir}

bcftools query -l ${vcf} | grep -E 'WRM' > wrm.txt
bcftools query -l ${vcf} | grep -v -E 'SBM|WRM|WP|TB' > bf_noartp.txt
cat ${proj}/out/sample_information.tsv | grep -E 'WRM' | awk '{if ($2 == "M") { print "bam/" $1 "/" $1 ".cfsm.bam" }}' > wrm_male.txt
cat ${proj}/out/sample_information.tsv | grep -v -E 'SBM|WRM' | awk '{if ($2 == "M") { print "bam/" $1 "/" $1 ".cfsm.bam" }}' > bf_male.txt
cat ${proj}/out/sample_information.tsv | grep -E 'WRM' | awk '{if ($2 == "F") { print "bam/" $1 "/" $1 ".cfsm.bam" }}' > wrm_female.txt
cat ${proj}/out/sample_information.tsv | grep -v -E 'SBM|WRM' | awk '{if ($2 == "F") { print "bam/" $1 "/" $1 ".cfsm.bam" }}' > bf_female.txt

vcftools --gzvcf ${vcf} --weir-fst-pop bf_male.txt --weir-fst-pop wrm_male.txt \
  --fst-window-size 10000 --fst-window-step 5000 --out ${workdir}/bf_vs_wrm_male

vcftools --gzvcf ${vcf} --weir-fst-pop bf_female.txt --weir-fst-pop wrm_female.txt \
  --fst-window-size 10000 --fst-window-step 5000 --out ${workdir}/bf_vs_wrm_female

vcftools --gzvcf ${vcf} --weir-fst-pop bf_noartp.txt --weir-fst-pop wrm.txt \
  --fst-window-size 10000 --fst-window-step 5000 --out ${workdir}/bf_vs_wrm_noartp
