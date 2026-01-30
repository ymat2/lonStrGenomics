#!/bin/bash

#SBATCH -o /dev/null
#SBATCH -e /dev/null

ref_dir=~/ref/lonStrDom2
accession=GCF_005870125.1

cd ${ref_dir}
apptainer exec /usr/local/biotools/h/hisat2:2.2.1--h1b792b2_3 hisat2-build ${accession}.fa lonStrDom2_hs2
