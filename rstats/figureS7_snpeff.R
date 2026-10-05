library(conflicted)
library(tidyverse)
source("./rstats/vcf_util.R")


## chr. 4A, HTR2C --------------------------------------------------------------

vep_4a = readr::read_tsv("out/snpeff/HTR2C.genes.txt", skip = 1)
vcf_4a = read_vcf("out/snpeff/HTR2C.annot.vcf")
af_4a = calc_allele_freq(vcf_4a) |>
  dplyr::mutate(dAF = abs(BF-WRM)) |>
  dplyr::left_join(vcf_4a |> dplyr::distinct(`#CHROM`, POS, INFO), by = c("#CHROM", "POS")) |>
  dplyr::ungroup() |>
  dplyr::mutate(INFO = stringr::str_split_i(INFO, "ANN=", 2)) |>
  tidyr::separate_rows(INFO, sep = ",") |>
  tidyr::separate(INFO, sep = "\\|", into = ann_field)
count_4a = af_4a |>
  dplyr::filter(dAF > 0.5) |>
  dplyr::count(Annotation) |>
  dplyr::mutate(chr = "chr. 4A")

## chr. 8, ALC44A5, ACADM ------------------------------------------------------

vep_8 = readr::read_tsv("out/snpeff/SLC44A5.genes.txt", skip = 1)
vcf_8 = read_vcf("out/snpeff/SLC44A5.annot.vcf")
af_8 = calc_allele_freq(vcf_8) |>
  dplyr::mutate(dAF = abs(BF-WRM)) |>
  dplyr::left_join(vcf_8 |> dplyr::distinct(`#CHROM`, POS, INFO), by = c("#CHROM", "POS")) |>
  dplyr::ungroup() |>
  dplyr::mutate(INFO = stringr::str_split_i(INFO, "ANN=", 2)) |>
  tidyr::separate_rows(INFO, sep = ",") |>
  tidyr::separate(INFO, sep = "\\|", into = ann_field)
count_8 = af_8 |>
  dplyr::filter(dAF > 0.5) |>
  dplyr::count(Annotation) |>
  dplyr::mutate(chr = "chr. 8")


## Visualization ---------------------------------------------------------------

snpeff = dplyr::bind_rows(count_4a, count_8) |>
  dplyr::mutate(Annotation = dplyr::case_when(
    stringr::str_detect(Annotation, "UTR") ~ "UTR",
    stringr::str_detect(Annotation, "downstream|intergenic|upstream") ~ "Intergenic",
    .default = stringr::str_replace_all(Annotation, "_", " ") |>
      stringr::str_to_sentence() |>
      stringr::str_remove_all(" variant")
  )) |>
  ggplot() +
  aes(x = chr, y = n) +
  geom_col(position = "fill", aes(fill = Annotation)) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_fill_viridis_d(option = "H") +
  labs(y = "Proportion of variant annotation") +
  theme_classic() +
  theme(
    axis.title.x = element_blank()
  )
snpeff
ggsave("images/snpeff.png", snpeff, h = 3, w = 4)
