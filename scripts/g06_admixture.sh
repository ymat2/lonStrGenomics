#!/bin/bash
#SBATCH -a 1-8
#SBATCH --mem 32G
#SBATCH -o /dev/null
#SBATCH -e /dev/null


shopt -s expand_aliases
alias admixture="apptainer exec /usr/local/biotools/a/admixture:1.3.0--0 admixture"

workdir=~/vocal-learning/structure
bed=lonchura.snp.bed

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

admixture --cv ${bed} ${SLURM_ARRAY_TASK_ID} | tee log${SLURM_ARRAY_TASK_ID}.out
grep -h CV log*.out > CV-error.txt
