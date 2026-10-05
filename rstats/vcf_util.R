library(conflicted)
library(tidyverse)


read_vcf = function(file) {
  vcf = readr::read_tsv(file, comment = "##") |>
    tidyr::pivot_longer(dplyr::starts_with("bam/"), names_to = "sample", values_to = "GT") |>
    dplyr::mutate(
      sample = stringr::str_split(sample, "/", simplify = TRUE)[,2],
      GT = stringr::str_split(GT, ":", simplify = TRUE)[,1],
      POS = as.character(POS),
      category = dplyr::case_when(
        stringr::str_detect(sample, "WRM") ~ "WRM",
        stringr::str_detect(sample, "SBM") ~ "SBM",
        .default = "BF") |> forcats::as_factor()
    ) |>
    dplyr::select(`#CHROM`, POS, REF, ALT, QUAL, INFO, sample, GT, category)
  return(vcf)
}


calc_allele_freq = function(vcf) {
  vcf2af = vcf |>
    dplyr::select(`#CHROM`, POS, sample, GT, category) |>
    dplyr::mutate(GT = dplyr::na_if(GT, "./.")) |>
    dplyr::mutate(AF = dplyr::case_when(GT == "1/1" ~ 1, GT == "0/1" ~ 0.5, .default = 0)) |>
    dplyr::group_by(`#CHROM`, POS, category) |>
    dplyr::summarise(mean_AF = mean(AF)) |>
    tidyr::pivot_wider(names_from = category, values_from = mean_AF)
  return(vcf2af)
}


ggaf = function(af_table) {
  ggplot(af_table) +
    aes(x = POS, y = diff_AF) +
    geom_col(aes(fill = VAR)) +
    scale_fill_manual(
      values = c("O" = "#999999", "M" = "#D73027"),
      labels = c("O" = "Other", "M" = "Missense")
    ) +
    labs(
      x = "",
      y = expression(paste(Delta, "Allele Freq.", sep=" ")),
      fill = "Variant"
    ) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_blank(),
      panel.background = element_blank()
    )
}


ggtile = function(vcf) {
  ggplot(vcf) +
    aes(x = POS, y = sample) +
    geom_tile(aes(fill = GT)) +
    scale_fill_viridis_d(option = "cividis") +
    labs(x = "", y = "", fill = "GT") +
    theme_test(base_size = 16) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 6)
    )
}


ann_field = c(
  "Allele",
  "Annotation",
  "Putative_impact",
  "Gene_Name",
  "Gene_ID",
  "Feature_type",
  "Feature_ID",
  "Transcript_biotype",
  "Rank_total",
  "HGVS_c",
  "HGVS_p",
  "cDNA_position_cDNA_length",
  "CDS_position_CDS_length",
  "Protein_position_Protein_length",
  "Distance_to_feature",
  "Log"
)
