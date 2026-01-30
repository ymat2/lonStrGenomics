library(conflicted)
library(tidyverse)
source("./rstats/vcf_util.R")


## chr. 4A, HTR2C --------------------------------------------------------------

vep = readr::read_tsv("out/snpeff/HTR2C.genes.txt", skip = 1)
vcf = read_vcf("out/snpeff/HTR2C.annot.vcf")
af = calc_allele_freq(vcf) |>
  dplyr::mutate(dAF = abs(BF-WRM)) |>
  dplyr::left_join(vcf |> dplyr::distinct(`#CHROM`, POS, INFO), by = c("#CHROM", "POS"))

## chr. 8, ALC44A5, ACADM ------------------------------------------------------

vep = readr::read_tsv("out/snpeff/SLC44A5.genes.txt", skip = 1)
vcf = read_vcf("out/snpeff/SLC44A5.annot.vcf")
