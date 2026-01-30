#!/bin/bash
#SBATCH -o /dev/null
#SBATCH -e /dev/null


shopt -s expand_aliases
alias bwa="apptainer exec /usr/local/biotools/b/bwa:0.7.17--h5bf99c6_8 bwa"

ref_dir=~/ref/lonStrDom2
accession=GCF_005870125.1

mkdir -p ${ref_dir}
cd ${ref_dir}
datasets download genome accession ${accession} --include gff3,gtf,genome
unzip ncbi_dataset.zip
genome=$(ls ncbi_dataset/data/${accession}/*.fna)
gff=$(ls ncbi_dataset/data/${accession}/*.gff)
gtf=$(ls ncbi_dataset/data/${accession}/*.gtf)

mv ${genome} ./${accession}.fa
mv ${gff} ./${accession}.gff
mv ${gtf} ./${accession}.gtf

rm ncbi_dataset.zip README.md
rm -r ncbi_dataset

apptainer exec /usr/local/biotools/b/bwa:0.7.17--h5bf99c6_8 bwa index ${accession}.fa
