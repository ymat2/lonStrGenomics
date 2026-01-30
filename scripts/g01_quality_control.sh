#!/bin/bash
#SBATCH -a 1-59
#SBATCH -o /dev/null
#SBATCH -e /dev/null


raw_data=~/raw_data/songbird2

samples=($(ls ${raw_data} | sort -V))
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

shopt -s expand_aliases
alias fastp="apptainer exec /usr/local/biotools/f/fastp:0.23.4--h5f740d0_0 fastp"

proj=~/vocal-learning
cd ${proj}

[ ! -e uncleaned/${sample} ] && mkdir -p uncleaned/${sample}
[ ! -e clean_data/${sample} ] && mkdir -p clean_data/${sample}
[ ! -e out/qc ] && mkdir -p out/qc

paired_1=($(ls ${raw_data}/${sample}/*_1.fq.gz))
paired_2=($(ls ${raw_data}/${sample}/*_2.fq.gz))

cat ${paired_1[@]} > uncleaned/${sample}/${sample}_1.fq.gz
cat ${paired_2[@]} > uncleaned/${sample}/${sample}_2.fq.gz

fastp \
  -i uncleaned/${sample}/${sample}_1.fq.gz \
  -I uncleaned/${sample}/${sample}_2.fq.gz \
  -o clean_data/${sample}/${sample}_clean_1.fq.gz \
  -O clean_data/${sample}/${sample}_clean_2.fq.gz \
  -h out/qc/${sample}_qc.html \
  -q 30 -u 30 -f 1 -F 1

rm -r uncleaned/${sample}
