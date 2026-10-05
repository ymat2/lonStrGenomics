#!/bin/bash
#SBATCH --job-name g07_supple_dxy
#SBATCH -o log/%x.log
#SBATCH -e log/%x.log


shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"

proj=~/lonchura
workdir=${proj}/selection/dxy
reference=~/ref/lonStrDom2/GCF_005870125.1.fa
vcf1=${workdir}/NC_042574.1_30M_last.vcf.gz
vcf2=${workdir}/NC_042570.1_19M_20M.vcf.gz

[ ! -e ${workdir} ] && mkdir -p ${workdir}
cd ${workdir}


##### Joint call #####

ls ${proj}/bam/*/*.cfsm.bam > bam.list

bcftools mpileup -f ${reference} --bam-list bam.list --regions NC_042574.1:30000000- |\
  bcftools call -m -Oz > ${vcf1}
bcftools index ${vcf1}

bcftools mpileup -f ${reference} --bam-list bam.list --regions NC_042570.1:19000000-20000000 |\
  bcftools call -m -Oz > ${vcf2}
bcftools index ${vcf2}


##### Pixy #####

bcftools query -l ${vcf1} | grep -E 'WRM' | awk '{print $1 "\tWRM" }' > ${workdir}/populations_file.txt
bcftools query -l ${vcf1} | grep -v -E 'WRM|SBM' | awk '{print $1 "\tBF" }' >> ${workdir}/populations_file.txt

. ${proj}/.venv/bin/activate

pixy --stats dxy \
  --vcf ${vcf1} \
  --populations ${workdir}/populations_file.txt \
  --window_size 10000 \
  --output_folder ${workdir} \
  --output_prefix NC_042574.1_30M_last

pixy --stats dxy \
  --vcf ${vcf2} \
  --populations ${workdir}/populations_file.txt \
  --window_size 10000 \
  --output_folder ${workdir} \
  --output_prefix NC_042570.1_19M_20M

deactivate
