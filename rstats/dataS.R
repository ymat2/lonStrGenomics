library(conflicted)
library(tidyverse)
library(readxl)
source("./rstats/annot.R")


## Sample information of genome sequence ---------------------------------------

sample_info = readr::read_tsv("out/sample_information_excel.tsv") |>
  dplyr::rename(sample = sample_id) |>
  dplyr::mutate(sex = dplyr::if_else(is.na(sex), "Un", sex))

sex_chrom_depth = readr::read_tsv("out/sex_chrom_depth.tsv") |>
  dplyr::inner_join(sample_info, by = "sample") |>
  dplyr::mutate(ratio = meandp_sex_chrom/meandp) |>
  dplyr::mutate(sex_from_depth = dplyr::if_else(ratio > .75, "M", "F")) |>
  dplyr::rename(sex_from_PCR = sex) |>
  dplyr::select(sample, sex_from_depth, species, colony)

flag_summary = readr::read_tsv("out/flag_summary.tsv")
summary = readr::read_tsv("out/mapping_summary.tsv") |>
  dplyr::inner_join(flag_summary, by = "sample") |>
  dplyr::inner_join(sex_chrom_depth, by = "sample") |>
  dplyr::select(sample, sex_from_depth, species, colony, num_reads, prop_mapped, meandp, meancov) |>
  dplyr::rename(
    "Sample" = sample,
    "Sex" = sex_from_depth,
    "Species" = species,
    "Breeding colony" = colony,
    "Number of reads" = num_reads,
    "Properly mapped reads (%)" = prop_mapped, 
    "Mean depth" = meandp, 
    "Mean coverage (%)" = meancov
  ) |>
  dplyr::arrange(Sample)

sample_info |>
  dplyr::mutate(subpop = stringr::str_remove(sample_id, "[0-9]+$")) |>
  dplyr::group_by(subpop, sex) |>
  dplyr::summarise(N = dplyr::n()) |>
  tidyr::pivot_wider(names_from = sex, values_from = N, values_fill = 0)

summary |>
  dplyr::group_by(species) |>
  dplyr::summarise(across(where(is.numeric), mean)) |>
  print()


## FST top 0.1% genes list -----------------------------------------------------

acc2chr = readr::read_tsv("out/sequence_report.tsv") |>
  dplyr::rename(CHROM = `RefSeq seq accession`, chr = `Chromosome name`) |>
  dplyr::select(CHROM, chr)

chr_levels = c("1", "1A", "2", "3", "4", "4A", as.character(seq(5,15)), 
               as.character(seq(17,29)), "Z", "MT")

fst01 = readr::read_tsv("out/fst/bf_vs_wrm.windowed.weir.fst") |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_relevel(chr, chr_levels)) |>
  dplyr::filter(!chr %in% c("Un", "MT")) |>
  dplyr::group_by(chr) |>
  dplyr::mutate(Z_FST = (WEIGHTED_FST - mean(WEIGHTED_FST))/sd(WEIGHTED_FST)) |>
  dplyr::ungroup() |>
  dplyr::slice_max(order_by = Z_FST, prop = 0.001) |>
  annot(seqnames = "CHROM", start = "BIN_START", end = "BIN_END") |>
  dplyr::relocate(chr, .after = CHROM) |>
  dplyr::arrange(chr, BIN_START)


## PI --------------------------------------------------------------------------

bf_pi = readr::read_tsv("out/pi/bf.windowed.pi") |>
  dplyr::rename(BF = PI) |>
  dplyr::select(!N_VARIANTS)

wrm_pi = readr::read_tsv("out/pi/wrm.windowed.pi") |>
  dplyr::rename(WRM = PI) |>
  dplyr::select(!N_VARIANTS)

pi_ratio = bf_pi |>
  dplyr::inner_join(wrm_pi, by = c("CHROM", "BIN_START", "BIN_END")) |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(
    chr = forcats::fct_inorder(chr),
    log_pi_ratio = -log10(BF/WRM)
  ) |>
  dplyr::filter(!chr %in% c("Un", "MT"))

pi0025_genes_bf = pi_ratio |>
  dplyr::slice_max(order_by = log_pi_ratio, prop = .0025) |>
  annot(seqnames = "CHROM", start = "BIN_START", end = "BIN_END") |>
  tidyr::drop_na(gene) |>
  dplyr::group_by(gene) |>
  dplyr::slice_max(order_by = log_pi_ratio, n = 1)

pi0025_genes_wrm = pi_ratio |>
  dplyr::slice_min(order_by = log_pi_ratio, prop = .0025) |>
  annot(seqnames = "CHROM", start = "BIN_START", end = "BIN_END") |>
  tidyr::drop_na(gene) |>
  dplyr::group_by(gene) |>
  dplyr::slice_min(order_by = log_pi_ratio, n = 1)

pi005_genes = dplyr::bind_rows(pi0025_genes_bf, pi0025_genes_wrm) |>
  dplyr::arrange(chr, BIN_START) |>
  dplyr::relocate(chr, .after = CHROM) |>
  dplyr::rename(PI_BF = BF, PI_WRM = WRM)


pi0025_genes_bf$gene |> unique() |> length()
pi0025_genes_wrm$gene |> unique() |> length()
readr::write_csv(pi005_genes, "docs/tables/pi005_genes.csv")


## RNA-seq ---------------------------------------------------------------------

degs = readr::read_tsv("out/diencephalon_stringtie_TCC.tsv") |>
  dplyr::mutate(m.value = -m.value) |>
  dplyr::mutate(upregulated_in = dplyr::case_when(
    q.value < 0.05 & m.value < -1 ~ "WRM",
    q.value < 0.05 & m.value > 1 ~ "BF",
    .default = "N"
  )) |>
  dplyr::mutate(symbol = stringr::str_split(gene_id, "\\|", simplify = TRUE)[,2]) |>
  dplyr::filter(upregulated_in != "N") |>
  dplyr::select(symbol, a.value, m.value, p.value, q.value, upregulated_in)


## Write excel ---\\\-----------------------------------------------------------

writexl::write_xlsx(
  list(
    "DataS1" = summary, 
    "DataS2" = fst01, 
    "DataS3" = pi005_genes,
    "DataS4" = degs
    ),
  path = "docs/manuscripts/supplementary_data.xlsx",
  format_headers = FALSE
)
