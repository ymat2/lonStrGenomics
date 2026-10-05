#!/bin/bash
#SBATCH --mem 32G
#SBATCH -a 1-100
#SBATCH --job-name g06_fastsimcoal2
#SBATCH -o /dev/null
#SBATCH -e /dev/null

#shopt -s expand_aliases
#alias fsc2="apptainer exec /usr/local/biotools/f/fastsimcoal2:27093--hdfd78af_0 fsc27093"

sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

proj=~/lonchura
sfs=${proj}/fsc2/output/fastsimcoal2/lonchura_jointMAFpop1_0.obs
tpl=${proj}/fsc2/lonchura.tpl
est=${proj}/fsc2/lonchura.est
prefix=run${SLURM_ARRAY_TASK_ID}
#prefix=test

workdir=${proj}/fsc2/${prefix}
[ ! -e ${workdir} ] && mkdir -p ${workdir}
cd ${workdir}

cp ${sfs} ${workdir}/${prefix}_jointMAFpop1_0.obs
cp ${tpl} ${workdir}/${prefix}.tpl
cp ${est} ${workdir}/${prefix}.est
fsc28 -t ${prefix}.tpl -n 10000 -m -e ${prefix}.est -L 40 -q -c0 -E 1 -M
