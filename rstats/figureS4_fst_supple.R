library(conflicted)
library(tidyverse)
source("./rstats/common_settings.R")
source("./rstats/ggman.R")


acc2chr= readr::read_tsv("out/sequence_report.tsv") |>
  dplyr::rename(CHROM = `RefSeq seq accession`, chr = `Chromosome name`) |>
  dplyr::select(CHROM, chr)

chr_levels = c("1", "1A", "2", "3", "4", "4A", as.character(seq(5,15)),
               as.character(seq(17,29)), "Z", "MT")

fst = readr::read_tsv("out/fst/bf_vs_wrm.windowed.weir.fst") |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_relevel(chr, chr_levels)) |>
  dplyr::filter(!chr %in% c("Un", "MT")) |>
  dplyr::group_by(chr) |>
  dplyr::mutate(Z_FST = (WEIGHTED_FST - mean(WEIGHTED_FST))/sd(WEIGHTED_FST)) |>
  dplyr::ungroup()

fst_male = readr::read_tsv("out/fst/bf_vs_wrm_male.windowed.weir.fst") |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_relevel(chr, chr_levels)) |>
  dplyr::filter(!chr %in% c("Un", "MT")) |>
  dplyr::group_by(chr) |>
  dplyr::mutate(Z_FST = (WEIGHTED_FST - mean(WEIGHTED_FST))/sd(WEIGHTED_FST)) |>
  dplyr::ungroup()

fst_female = readr::read_tsv("out/fst/bf_vs_wrm_female.windowed.weir.fst") |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_relevel(chr, chr_levels)) |>
  dplyr::filter(!chr %in% c("Un", "MT")) |>
  dplyr::group_by(chr) |>
  dplyr::mutate(Z_FST = (WEIGHTED_FST - mean(WEIGHTED_FST))/sd(WEIGHTED_FST)) |>
  dplyr::ungroup()


ga = ggmanh(
  fst,
  chr = "chr", bp = "BIN_START", p = "Z_FST",
  col = c("#666666", "#BBBBBB"),
  chrlabs = c(1, "1A", 2, 3, 4, "4A", 5:10, 12, 14, 17, 20, 29, "Z"),
  xlab = "Chromosome",
  ylab = expression(paste("Z", italic(F)[ST])),
  ylim = c(NA, NA),
  title = "All (BF=52 vs WRM=31)",
  base_size = BASESIZE
)

gm = ggmanh(
  fst_male,
  chr = "chr", bp = "BIN_START", p = "Z_FST",
  col = c("#666666", "#BBBBBB"),
  chrlabs = c(1, "1A", 2, 3, 4, "4A", 5:10, 12, 14, 17, 20, 29, "Z"),
  xlab = "Chromosome",
  ylab = expression(paste("Z", italic(F)[ST])),
  ylim = c(NA, NA),
  title = "Male (BF=34 vs WRM=14)",
  base_size = BASESIZE
)

gf = ggmanh(
  fst_female,
  chr = "chr", bp = "BIN_START", p = "Z_FST",
  col = c("#666666", "#BBBBBB"),
  chrlabs = c(1, "1A", 2, 3, 4, "4A", 5:10, 12, 14, 17, 20, 29, "Z"),
  xlab = "Chromosome",
  ylab = expression(paste("Z", italic(F)[ST])),
  ylim = c(NA, NA),
  title = "Female (BF=18 vs WRM=17)",
  base_size = BASESIZE
)

g1 = cowplot::plot_grid(ga, gm, gf, nrow = 3, align = "v", axis = "lr")
ggsave("images/fst_supple_male_female.png", g1, w = 7, h = 6)
