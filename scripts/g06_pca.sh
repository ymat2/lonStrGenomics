#!/bin/bash
#SBATCH -o /dev/null
#SBATCH -e /dev/null

shopt -s expand_aliases
alias plink="apptainer exec /usr/local/biotools/p/plink2:2.00a5--h4ac6f70_0 plink2"

workdir=~/vocal-learning/structure
prefix=lonchura.snp

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

plink --bfile ${prefix} --pca --allow-extra-chr --double-id --out ${prefix}.pca
