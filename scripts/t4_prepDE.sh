#!/bin/bash
#SBATCH -o /dev/null
#SBATCH -e /dev/null


proj=~/vocal-learning
cd ${proj}

/usr/bin/python3 src/summary_qc_result.py -i out/qc_rna

/usr/bin/python3 src/prepDE.py \
  -i rnaseq_bam \
  -p "D$" \
  -g out/diencephalon_gene_count_matrix.csv \
  -t out/diencephalon_transcript_count_matrix.csv \
  -l 148  # average read length

/usr/bin/python3 src/extractTPM.py \
  -i rnaseq_bam \
  -p "D$" \
  -g out/diencephalon_gene_TPM_matrix.csv
