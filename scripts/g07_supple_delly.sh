#!/bin/bash
#SBATCH --mem 63G
#SBATCH --job-name g07_supple_delly
#SBATCH -o log/%x.log
#SBATCH -e log/%x.log


shopt -s expand_aliases
alias samtools="apptainer exec /usr/local/biotools/s/samtools:1.18--h50ea8bc_1 samtools"
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias delly="apptainer exec /usr/local/biotools/d/delly:1.7.3--hd6466ae_0 delly"
alias vcftools="apptainer exec /usr/local/biotools/v/vcftools:0.1.16--h9a82719_5 vcftools"

reference=~/ref/lonStrDom2/GCF_005870125.1.fa
proj=~/lonchura
workdir=${proj}/sv
samples=($(ls ${proj}/bam | sort -V))

cd ${proj}

[ ! -e ${workdir}/bam ] && mkdir -p ${workdir}/bam

## Region of interests
## - NC_042574.1:30000000-last
## - NC_042570.1:19000000-20000000

##### Call SV using Delly #####

for sample in ${samples[@]}; do
  samtools view -b ${proj}/bam/${sample}/${sample}.cfsm.bam NC_042574.1:30000000 NC_042570.1:19000000-20000000 -o ${workdir}/bam/${sample}.partial.bam
  samtools index ${workdir}/bam/${sample}.partial.bam
done

delly call -g ${reference} -q 30 $(ls ${workdir}/bam/*.partial.bam) -o ${workdir}/tmp.bcf
bcftools filter -i 'FILTER="PASS"' ${workdir}/tmp.bcf -Ob -o ${workdir}/lonchura.sv.bcf
bcftools index ${workdir}/lonchura.sv.bcf && rm ${workdir}/tmp.bcf*


##### SV comparison between BF and WRM #####

bcftools query -l ${workdir}/lonchura.sv.bcf | grep -E 'WRM' > ${workdir}/wrm.txt
bcftools query -l ${workdir}/lonchura.sv.bcf | grep -v -E 'SBM|WRM' > ${workdir}/bf.txt

vcftools --bcf ${workdir}/lonchura.sv.bcf \
  --maf 0.01 \
  --max-missing 0.1 \
  --weir-fst-pop ${workdir}/bf.txt \
  --weir-fst-pop ${workdir}/wrm.txt \
  --out ${workdir}/bf_vs_wrm.sv


##### DEscribe read depth of the regions #####

samtools depth -r NC_042574.1:30000000- -a --min-MQ 30 $(ls ${proj}/bam/*/*.cfsm.bam | grep -E "WRM") |\
  awk '{sum = 0; for (i = 3; i <= NF; i++) {sum += $i} print $1 "\t" $2 "\t" sum / (NF - 2) "\tWRM"}' \
  > ${workdir}/NC_042574.1_30M_last_depth.tsv
samtools depth -r NC_042574.1:30000000- -a --min-MQ 30 $(ls ${proj}/bam/*/*.cfsm.bam | grep -E -v "WRM|SBM") |\
  awk '{sum = 0; for (i = 3; i <= NF; i++) {sum += $i} print $1 "\t" $2 "\t" sum / (NF - 2) "\tBF"}' \
  >> ${workdir}/NC_042574.1_30M_last_depth.tsv

samtools depth -r NC_042570.1:19000000-20000000 -a --min-MQ 30 $(ls ${proj}/bam/*/*.cfsm.bam | grep -E "WRM") |\
  awk '{sum = 0; for (i = 3; i <= NF; i++) {sum += $i} print $1 "\t" $2 "\t" sum / (NF - 2) "\tWRM"}' \
  > ${workdir}/NC_042570.1_19M_20M_depth.tsv
samtools depth -r NC_042570.1:19000000-20000000 -a --min-MQ 30 $(ls ${proj}/bam/*/*.cfsm.bam | grep -E -v "WRM|SBM") |\
  awk '{sum = 0; for (i = 3; i <= NF; i++) {sum += $i} print $1 "\t" $2 "\t" sum / (NF - 2) "\tBF"}' \
  >> ${workdir}/NC_042570.1_19M_20M_depth.tsv
