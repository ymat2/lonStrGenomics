library(conflicted)
library(tidyverse)
source("./rstats/common_settings.R")


acc2chr= readr::read_tsv("out/sequence_report.tsv") |>
  dplyr::rename(CHROM = `RefSeq seq accession`, chr = `Chromosome name`) |>
  dplyr::select(CHROM, chr)

chr_levels = c("1", "1A", "2", "3", "4", "4A", as.character(seq(5,15)), 
               as.character(seq(17,29)), "Z", "MT")


## FST -------------------------------------------------------------------------

fst = readr::read_tsv("out/fst/bf_vs_wrm.windowed.weir.fst") |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_relevel(chr, chr_levels)) |>
  dplyr::filter(!chr %in% c("Un", "MT")) |>
  dplyr::group_by(chr) |>
  dplyr::mutate(Z_FST = (WEIGHTED_FST - mean(WEIGHTED_FST))/sd(WEIGHTED_FST)) |>
  dplyr::ungroup()


## Tajima's D ------------------------------------------------------------------

bf_tD = readr::read_tsv("out/tajimasD/bf.Tajima.D") |>
  dplyr::rename(BF = TajimaD) |>
  dplyr::select(!N_SNPS)
wrm_tD = readr::read_tsv("out/tajimasD/wrm.Tajima.D") |>
  dplyr::rename(WRM = TajimaD) |>
  dplyr::select(!N_SNPS)
tD = bf_tD |>
  dplyr::inner_join(wrm_tD, by = c("CHROM", "BIN_START")) |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_inorder(chr)) |>
  dplyr::filter(!chr %in% c("Un", "MT"))


## Comparison ------------------------------------------------------------------

fst_td_high_fst = tD |>
  dplyr::mutate(BIN_START = BIN_START + 1) |>
  dplyr::inner_join(fst, by = c("CHROM", "chr", "BIN_START")) |>
  dplyr::mutate(cat = dplyr::if_else(Z_FST >= quantile(fst$Z_FST, .999), "high", "low"))

wil_wrm = wilcox.test(
  fst_td_high_fst |> dplyr::filter(cat == "high") |> dplyr::pull(WRM),
  fst_td_high_fst |> dplyr::filter(cat == "low") |> dplyr::pull(WRM),
  paired = FALSE
)
wil_wrm$p.value

wil_bf = wilcox.test(
  fst_td_high_fst |> dplyr::filter(cat == "high") |> dplyr::pull(BF),
  fst_td_high_fst |> dplyr::filter(cat == "low") |> dplyr::pull(BF),
  paired = FALSE
)
wil_bf$p.value

fst_td_wrm = ggplot(fst_td_high_fst) +
  aes(x = cat, y = WRM, color = cat) +
  geom_boxplot(outliers = TRUE, fill = NA) +
  ggplot2::annotate("text", x = 1.5, y = 4.5, label = expression(paste(italic(P), " = 1.207 × ", , 10^{-7})), color = "#444444") +
  scale_x_discrete(labels = c(expression(paste(italic(F)[ST], " top 0.1%")), "Others")) +
  scale_y_continuous(limits = c(-3, 5)) +
  scale_color_manual(values = c("high" = "#5ab4ac", "low" = "#999999")) +
  labs(
    x = "Genomic region",
    y = expression(paste("Tajima's ", italic(D))),
    title = "WRM"
  ) +
  theme_test(base_size = BASESIZE) +
  theme(legend.position = "none", plot.title = element_text(hjust = .5))
fst_td_wrm

fst_td_bf = ggplot(fst_td_high_fst) +
  aes(x = cat, y = BF, color = cat) +
  geom_boxplot(outliers = TRUE, fill = NA) +
  ggplot2::annotate("text", x = 1.5, y = 4.5, label = expression(paste(italic(P), " = 1.338 × ", , 10^{-24})), color = "#444444") +
  scale_x_discrete(labels = c(expression(paste(italic(F)[ST], " top 0.1%")), "Others")) +
  scale_y_continuous(limits = c(-3, 5)) +
  scale_color_manual(values = c("high" = "#d8b365", "low" = "#999999")) +
  labs(
    x = "Genomic region",
    title = "BF"
  ) +
  theme_test(base_size = BASESIZE) +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = .5),
    axis.title.y = element_blank()
  )
fst_td_bf

.df = fst_td_high_fst |> dplyr::filter(cat == "high")
res = lm(data = .df, BF ~ WRM) |> summary()
round(res$adj.r.squared, digits = 3)
res$coefficients[2, 4]

td_wrm_bf = ggplot(.df) + 
  aes(BF, WRM) + 
  geom_point(size = 3, alpha = .5, shape = 16) + 
  stat_smooth(method = "lm", formula = y~x, se = FALSE, color = "#333333") +
  ggplot2::annotate(geom = "text", x = 2, y = 1, hjust = .8, vjust = -1, label = expression(paste(italic(r)^2, " = ", "0.162"))) +
  ggplot2::annotate(geom = "text", x = 2, y = 1, hjust = .8, vjust = 1, label = expression(paste(italic(P), " = ", "1.495 × ", 10^{-5}))) +
  labs(
    x = expression(paste("Tajima's ", italic(D), " in BF")),
    y = expression(paste("Tajima's ", italic(D), " in WRM")),
    title = expression(paste(italic(F)[ST], " top 0.1% regions"))
  ) +
  theme_test(base_size = BASESIZE) +
  theme(plot.title = element_text(hjust = .5))
td_wrm_bf

p = cowplot::plot_grid(
  fst_td_wrm, fst_td_bf, td_wrm_bf,
  nrow = 1, rel_widths = c(2, 2, 3),
  align = "h", axis = "tb",
  labels = c("a", "", "b"), label_size = LABELSIZE,
  scale = .95
)
ggsave("images/comparison_TajimaD_by_Fst.png", p, w = 9, h = 4, bg = "#FFFFFF")

