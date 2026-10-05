#!/bin/bash
#SBATCH -o /dev/null
#SBATCH -e /dev/null


proj=~/lonchura
workdir=${proj}/snpeff
vcf=${proj}/vcf/lonchura.vcf.gz
snpvcf=${proj}/vcf/lonchura.snp.vcf.gz

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

bed=${proj}/snpeff/zfst001.bed
bcftools view --min-af 0.05 -m2 -M2 --regions-file ${bed} -Ov ${vcf} | snpEff lonStrDom2 - -stats zfst001 > zfst001.maf005.annot.vcf
rm tmp.vcf zfst001


# WEIGHTED_FST > 0.5

bed=${proj}/selection/fst/bf_vs_wrm.fst05.txt
sbm=/home/ymat2/lonchura/selection/sbm.txt
bcftools view --min-ac 1 --samples-file ^${sbm} -Ov ${vcf} > tmp.vcf
bedtools intersect -a tmp.vcf -b ${bed} -wa | snpEff lonStrDom2 - -stats fst05 > /dev/null
rm tmp.vcf


# MEAN_FST > 0.4

bed=${proj}/selection/fst/bf_vs_wrm.meanfst04.txt
sbm=/home/ymat2/lonchura/selection/sbm.txt
bcftools view --min-ac 1 --samples-file ^${sbm} -Ov ${vcf} > tmp.vcf
bedtools intersect -a tmp.vcf -b ${bed} -wa | snpEff lonStrDom2 - -stats meanfst04 > /dev/null
rm tmp.vcf


# FOXP2 = NC_042566.1: 38,453,818..38,853,445

gene_id=FOXP2
chr=NC_042566.1
bp_start=$((38453818-5000))
bp_end=$((38853445+5000))

vcftools --gzvcf ${vcf} \
  --chr ${chr} --from-bp ${bp_start} --to-bp ${bp_end} \
  --recode --out ${gene_id}

snpEff lonStrDom2 ${gene_id}.recode.vcf -stats ${gene_id} > ${gene_id}.annot.vcf
rm ${gene_id}.recode.vcf
