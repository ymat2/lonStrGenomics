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
  dplyr::rename(BF_tajD = TajimaD) |>
  dplyr::select(!N_SNPS)
wrm_tD = readr::read_tsv("out/tajimasD/wrm.Tajima.D") |>
  dplyr::rename(WRM_tajD = TajimaD) |>
  dplyr::select(!N_SNPS)
tD = bf_tD |>
  dplyr::inner_join(wrm_tD, by = c("CHROM", "BIN_START")) |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_inorder(chr)) |>
  dplyr::filter(!chr %in% c("Un", "MT"))


## Nucleotide diversity --------------------------------------------------------

bf_pi = readr::read_tsv("out/pi/bf.windowed.pi") |>
  dplyr::rename(BF_PI = PI) |>
  dplyr::select(!N_VARIANTS)
wrm_pi = readr::read_tsv("out/pi/wrm.windowed.pi") |>
  dplyr::rename(WRM_PI = PI) |>
  dplyr::select(!N_VARIANTS)
pi = bf_pi |>
  dplyr::inner_join(wrm_pi, by = c("CHROM", "BIN_START", "BIN_END")) |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_inorder(chr)) |>
  dplyr::filter(!chr %in% c("W", "Z", "Un", "MT"))

## Comparison ------------------------------------------------------------------

fst_td_high_fst_pi = tD |>
  dplyr::mutate(BIN_START = BIN_START + 1) |>
  dplyr::inner_join(fst, by = c("CHROM", "chr", "BIN_START")) |>
  dplyr::inner_join(pi, by = c("CHROM", "chr", "BIN_START", "BIN_END")) |>
  dplyr::mutate(cat = dplyr::if_else(Z_FST >= quantile(fst$Z_FST, .999), "high", "low"))

### comparison of Tajima's D

wil_wrm = wilcox.test(
  fst_td_high_fst_pi |> dplyr::filter(cat == "high") |> dplyr::pull(WRM_tajD),
  fst_td_high_fst_pi |> dplyr::filter(cat == "low") |> dplyr::pull(WRM_tajD),
  paired = FALSE
)
wil_wrm$p.value

fst_td_wrm = ggplot(fst_td_high_fst_pi) +
  aes(x = cat, y = WRM_tajD, color = cat) +
  geom_boxplot(outliers = TRUE, fill = NA) +
  scale_x_discrete(labels = c(expression(paste(italic(F)[ST], " top 0.1%")), "Others")) +
  scale_color_manual(values = c("high" = "#5ab4ac", "low" = "#999999")) +
  labs(
    x = "Genomic region",
    y = expression(paste("Tajima's ", italic(D))),
    title = "WRM",
    subtitle = expression(paste(italic(P), " = 1.207 × ", , 10^{-7}))
  ) +
  theme_test(base_size = BASESIZE) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5),
  )
fst_td_wrm


wil_bf = wilcox.test(
  fst_td_high_fst_pi |> dplyr::filter(cat == "high") |> dplyr::pull(BF_tajD),
  fst_td_high_fst_pi |> dplyr::filter(cat == "low") |> dplyr::pull(BF_tajD),
  paired = FALSE
)
wil_bf$p.value

fst_td_bf = ggplot(fst_td_high_fst_pi) +
  aes(x = cat, y = BF_tajD, color = cat) +
  geom_boxplot(outliers = TRUE, fill = NA) +
  scale_x_discrete(labels = c(expression(paste(italic(F)[ST], " top 0.1%")), "Others")) +
  scale_color_manual(values = c("high" = "#d8b365", "low" = "#999999")) +
  labs(
    x = "Genomic region",
    title = "BF",
    subtitle = expression(paste(italic(P), " = 1.319 × ", , 10^{-24}))
  ) +
  theme_test(base_size = BASESIZE) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5),
    axis.title.y = element_blank(),
  )
fst_td_bf

### comparison of Tajima's D and Fst top 0.1% regions

.df = fst_td_high_fst_pi |> dplyr::filter(cat == "high")
res = lm(data = .df, BF_tajD ~ WRM_tajD) |> summary()
round(res$adj.r.squared, digits = 3)
res$coefficients[2, 4]

.stat_result = expression(paste(italic(r)^2, " = ", "0.162, ", italic(P), " = ", "1.495 × ", 10^{-5}))
td_wrm_bf = ggplot(.df) +
  aes(BF_tajD, WRM_tajD) +
  geom_point(size = 2, alpha = .5, shape = 16) +
  stat_smooth(method = "lm", formula = y~x, se = FALSE, color = "#333333") +
  labs(
    x = expression(paste("Tajima's ", italic(D), " in BF")),
    y = expression(paste("Tajima's ", italic(D), " in WRM")),
    title = expression(paste(italic(F)[ST], " top 0.1% regions")),
    subtitle = .stat_result
  ) +
  theme_test(base_size = BASESIZE) +
  theme(
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5)
  )
td_wrm_bf

# all regions
res = lm(data = fst_td_high_fst_pi, BF_tajD ~ WRM_tajD) |> summary()
round(res$adj.r.squared, digits = 3)
res$coefficients[2, 4]

.stat_result = expression(paste(italic(r)^2, " = ", "0.041, ", italic(P), " < ", "2.2 × ", 10^{-16}))
td_wrm_bf_all = ggplot(fst_td_high_fst_pi) +
  aes(BF_tajD, WRM_tajD) +
  geom_point(size = 2, alpha = .5, shape = 16) +
  stat_smooth(method = "lm", formula = y~x, se = FALSE, color = "#eeeeee") +
  labs(
    x = expression(paste("Tajima's ", italic(D), " in BF")),
    y = expression(paste("Tajima's ", italic(D), " in WRM")),
    title = "All regions",
    subtitle = .stat_result
  ) +
  theme_test(base_size = BASESIZE) +
  theme(
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5)
  )
td_wrm_bf_all

## Tajima's D vs Pi in WRM

res = lm(data = .df, WRM_tajD ~ WRM_PI) |> summary()
round(res$adj.r.squared, digits = 3)
res$coefficients[2, 4]

.stat_result = expression(paste(italic(r)^2, " = ", "0.291, ", italic(P), " = ", "2.538 × ", 10^{-9}))
td_pi_wrm = ggplot(.df) +
  aes(WRM_tajD, WRM_PI) +
  geom_hline(yintercept = mean(pi$WRM_PI), linetype = "dashed") +
  geom_point(size = 2, color = colWRM, alpha = .5, shape = 16) +
  stat_smooth(method = "lm", formula = y~x, se = FALSE, color = "#333333") +
  labs(
    x = expression(paste("Tajima's ", italic(D), " in WRM")),
    y = expression(paste({pi}, " in WRM")),
    title = expression(paste(italic(F)[ST], " top 0.1% regions")),
    subtitle = .stat_result
  ) +
  theme_test(base_size = BASESIZE) +
  theme(
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5)
  )
td_pi_wrm


## Tajima's D vs Pi in BF

res = lm(data = .df, BF_tajD ~ BF_PI) |> summary()
round(res$adj.r.squared, digits = 3)
res$coefficients[2, 4]

.stat_result = expression(paste(italic(r)^2, " = ", "0.254, ", italic(P), " = ", "3.396 × ", 10^{-8}))
td_pi_bf = ggplot(.df) +
  aes(BF_tajD, BF_PI) +
  geom_hline(yintercept = mean(pi$BF_PI), linetype = "dashed") +
  geom_point(size = 2, color = colBF, alpha = .5, shape = 16) +
  stat_smooth(method = "lm", formula = y~x, se = FALSE, color = "#333333") +
  labs(
    x = expression(paste("Tajima's ", italic(D), " in BF")),
    y = expression(paste({pi}, " in BF")),
    title = expression(paste(italic(F)[ST], " top 0.1% regions")),
    subtitle = .stat_result
  ) +
  theme_test(base_size = BASESIZE) +
  theme(
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5)
  )
td_pi_bf


## Cowplot ---------------------------------------------------------------------

p1 = cowplot::plot_grid(
  fst_td_wrm, fst_td_bf, td_wrm_bf, td_wrm_bf_all, nrow = 1,
  rel_widths = c(2, 2, 3, 3),
  align = "vh",
  axis = "lrtb",
  labels = c("a", "", "b", "c"),
  label_size = LABELSIZE,
  scale = .95
)
p2 = cowplot::plot_grid(
  td_pi_wrm, td_pi_bf, NULL, nrow = 1,
  rel_widths = c(3, 3, 4),
  align = "vh",
  axis = "lrtb",
  labels = c("d", "e", ""),
  label_size = LABELSIZE,
  scale = .95
)
p = cowplot::plot_grid(p1, p2, nrow = 2)
ggsave("images/comparison_TajimaD_by_Fst.png", p, w = 7, h = 5, bg = "#FFFFFF")
