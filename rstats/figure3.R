library(conflicted)
library(tidyverse)
source("./rstats/common_settings.R")
source("./rstats/annot.R")


tpm = readr::read_csv("out/diencephalon_gene_TPM_matrix.csv")
degs = readr::read_tsv("out/diencephalon_stringtie_TCC.tsv") |>
  dplyr::mutate(m.value = -m.value) |>
  dplyr::mutate(estimatedDEG = dplyr::case_when(
    q.value < 0.05 & m.value < -1 ~ "WRM",
    q.value < 0.05 & m.value > 1 ~ "BF",
    .default = "N"
  )) |>
  dplyr::mutate(logP = -log10(p.value)) |>
  dplyr::mutate(symbol = stringr::str_split(gene_id, "\\|", simplify = TRUE)[,2]) |>
  dplyr::left_join(tpm, by = "gene_id")

degs |> dplyr::filter(estimatedDEG == "BF") |> nrow()
degs |> dplyr::filter(estimatedDEG == "WRM") |> nrow()


## Volcano plot ----------------------------------------------------------------

degs_goi = degs |> dplyr::filter(symbol %in% c(
  "CFAP54",
  "LHX8",
  "SLC44A5",
  "ACADM",
  "RABGGTB",
  "NAA20",
  "HTR2C",
  "IL13RA2"
))
  
vol = degs |>
  dplyr::arrange(desc(rank)) |>
  ggplot() +
  aes(x = m.value, y = logP) +
  geom_point(aes(color = estimatedDEG)) +
  ggrepel::geom_text_repel(
    data = degs_goi |> dplyr::filter(m.value < 0),
    aes(label = symbol),
    fontface = "italic",
    size = BASESIZE/2.83465,
    segment.size = .25,
    nudge_x = -8, nudge_y = .5, direction = "y", max.overlaps = Inf
    ) +
  ggrepel::geom_text_repel(
    data = degs_goi |> dplyr::filter(m.value > 0),
    aes(label = symbol),
    fontface = "italic",
    size = BASESIZE/2.83465,
    segment.size = .25,
    nudge_x = 9, nudge_y = .5, direction = "y", max.overlaps = Inf
  ) +
  geom_point(data = degs_goi, color = "#333333", shape = 1) +
  scale_color_manual(
    values = c("BF" = colBF, "N" = "#CCCCCC", "WRM" = colWRM),
    labels = c("BF" = "Higher expression in BF,", "N" = "Not significant,", "WRM" = "Higher expression in WRM")
    ) +
  labs(
    x = expression(paste(log[2], " Fold Change")), 
    y = expression(paste(-log[10], " ", italic(P), "-value")),
    color = "Expression difference in hypothalamus (FDR < 0.05, |logFC| > 1)"
    ) +
  theme_test(base_size = BASESIZE) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(size = BASESIZE),
    legend.key.size = unit(BASESIZE * 1.33, "pt"),
    legend.text = element_text(size = BASESIZE, margin = margin(l = 0)),
    axis.text = element_text(size = BASESIZE)
  ) +
  guides(color = guide_legend(nrow = 1, title.position = "top"))

pvol = cowplot::ggdraw() +
  cowplot::draw_plot(vol) +
  cowplot::draw_image("photos/wrm.nobg.png", scale = .2, x = -.35, y = .3) +
  cowplot::draw_image("photos/bf.nobg.png", scale = .2, x = .4, y = .3)
pvol

## Overlap with ZFST top 0.1% genes --------------------------------------------

acc2chr = readr::read_tsv("out/sequence_report.tsv") |>
  dplyr::rename(CHROM = `RefSeq seq accession`, chr = `Chromosome name`) |>
  dplyr::select(CHROM, chr)

chr_levels = c("1", "1A", "2", "3", "4", "4A", as.character(seq(5,15)), 
               as.character(seq(17,29)), "Z", "MT")

zfst_top01 = readr::read_tsv("out/fst/bf_vs_wrm.windowed.weir.fst") |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_relevel(chr, chr_levels)) |>
  dplyr::filter(!chr %in% c("Un", "MT")) |>
  dplyr::group_by(chr) |>
  dplyr::mutate(Z_FST = (WEIGHTED_FST - mean(WEIGHTED_FST))/sd(WEIGHTED_FST)) |>
  dplyr::ungroup() |>
  dplyr::slice_max(order_by = Z_FST, prop = .001) |>
  annot(seqnames = "CHROM", start = "BIN_START", end = "BIN_END")

genes = list(
  wrm = degs |> dplyr::filter(estimatedDEG == "WRM") |> dplyr::pull(symbol) |> unique(),
  bf = degs |> dplyr::filter(estimatedDEG == "BF") |> dplyr::pull(symbol) |> unique(),
  fst = zfst_top01 |> dplyr::pull(gene) |> unique()
)

pven = ggplot(data = dplyr::tibble(x = c(1,10), y = c(1,10))) +
  aes(x, y) +
  ggforce::geom_circle(aes(x0 = 9, y0 = 5, r = 2.7), fill = lighter(colBF), color = NA) +
  ggforce::geom_circle(aes(x0 = 2.5, y0 = 5, r = 1.1), fill = lighter(colWRM), color = NA) +
  ggforce::geom_circle(aes(x0 = 5, y0 = 5, r = 1.8), fill = "#CCCCCC66", color = NA) +
  ggplot2::annotate("text", x = 9, y = 5, label = length(genes$bf), size = BASESIZE, size.unit = "pt") +
  ggplot2::annotate("text", x = 2.5, y = 5, label = length(genes$wrm), size = BASESIZE, size.unit = "pt") +
  ggplot2::annotate("text", x = 5, y = 5, label = length(genes$fst), size = BASESIZE, size.unit = "pt") +
  ggplot2::annotate("text", x = 6.6, y = 5, label = "1", size = BASESIZE, size.unit = "pt") +
  ggplot2::annotate("text", x = 3.4, y = 5, label = "1", size = BASESIZE, size.unit = "pt") +
  ggplot2::annotate(
    "text", x = 9.5, y = 2, size = BASESIZE, size.unit = "pt",
    label = "Upregulated in BF",
    ) +
  ggplot2::annotate(
    "text", x = 2, y = 3.5, size = BASESIZE, size.unit = "pt",
    label = "Upregulated in WRM",
  ) +
  ggplot2::annotate(
    "text", x = 5, y = 7, size = BASESIZE, size.unit = "pt", parse = TRUE,
    label = as.character(expression(paste("Z", italic(F)[ST], " top 0.1%"))),
  ) +
  ggplot2::annotate("segment", x = 8, y = 8, xend = 6.6, yend = 5.2, linewidth = .25) +
  ggplot2::annotate("segment", x = 3, y = 8, xend = 3.4, yend = 5.2, linewidth = .25) +
  ggplot2::annotate(
    "label", x = 8, y = 8, size = BASESIZE, size.unit = "pt", parse = TRUE,
    fill = "#FFFFFF", border.color = NA,
    label = as.character(expression(italic(CFAP54)))
  ) +
  ggplot2::annotate(
    "label", x = 3, y = 8, size = BASESIZE, size.unit = "pt", parse = TRUE,
    fill = "#FFFFFF", border.color = NA,
    label = as.character(expression(italic(LHX8)))
  ) +
  coord_fixed(clip = "off") +
  theme_void()
pven

exp_cfap54 = degs |>
  dplyr::filter(symbol == "CFAP54") |>
  tidyr::pivot_longer(dplyr::matches("^B[0-9]+D$|^W[0-9]+D$"), names_to = "sample", values_to = "tpm") |>
  dplyr::mutate(group = dplyr::if_else(stringr::str_detect(sample, "^W"), "WRM", "BF")) |>
  dplyr::mutate(group = forcats::fct_relevel(group,  c("WRM", "BF"))) |>
  ggplot() +
  aes(x = group, y = tpm) +
  geom_bar(aes(fill = group), width = .5, alpha = .8, position = "dodge", stat = "summary", fun = "mean") +
  geom_point(shape = 21, fill = NA) +
  ggplot2::annotate("segment", x = 1, xend = 2, y = 3.9, yend = 3.9) +
  ggplot2::annotate("text", x = 1.5, y = 4, label = "*", size = BASESIZE, size.unit = "pt") +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_fill_manual(values = colors) +
  labs(y = "", title = "CFAP54") +
  theme_classic(base_size = BASESIZE) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = BASESIZE, face = "italic", hjust = .5),
    axis.title.x = element_blank(),
    axis.text = element_text(size = BASESIZE)
  )

exp_lhx8 = degs |>
  dplyr::filter(symbol == "LHX8") |>
  tidyr::pivot_longer(dplyr::matches("^B[0-9]+D$|^W[0-9]+D$"), names_to = "sample", values_to = "tpm") |>
  dplyr::mutate(group = dplyr::if_else(stringr::str_detect(sample, "^W"), "WRM", "BF")) |>
  dplyr::mutate(group = forcats::fct_relevel(group,  c("WRM", "BF"))) |>
  ggplot() +
  aes(x = group, y = tpm) +
  geom_bar(aes(fill = group), width = .5, alpha = .8, position = "dodge", stat = "summary", fun = "mean") +
  geom_point(shape = 21, fill = NA) +
  ggplot2::annotate("segment", x = 1, xend = 2, y = 49, yend = 49) +
  ggplot2::annotate("text", x = 1.5, y = 50, label = "*", size = BASESIZE, size.unit = "pt") +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_fill_manual(values = colors) +
  labs(y = "TPM", title = "LHX8") +
  theme_classic(base_size = BASESIZE) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = BASESIZE, face = "italic", hjust = .5),
    axis.title.x = element_blank(),
    axis.text = element_text(size = BASESIZE)
  )


## Serotonin-related genes -----------------------------------------------------

serotonin_related_genes = c(
  "TPH1", 
  "TPH2",
  "DDC",
  # "SLC18A1",  # lost
  "SLC18A2", 
  "LOC110474872",  # MAOA 
  "MAOB", 
  "SLC6A4",  # SERT
  "HTR6",
  "HTR1A",
  "HTR1F",
  "HTR5A",
  "HTR2C",
  "HTR1E",
  "HTR7",
  "HTR4",
  "HTR1B",
  "HTR2A",
  "HTR2B",
  "LOC110470479",  # HTR3A
  "HTR1D"
)

serotonin_genes_expression = degs |> 
  dplyr::filter(symbol %in% serotonin_related_genes) |>
  dplyr::mutate(
    symbol = symbol |>
      stringr::str_replace("LOC110474872", "LOC110474872 (MAOA)") |>
      stringr::str_replace("SLC18A2", "SLC18A2 (VMAT2)") |>
      stringr::str_replace("SLC6A4", "SLC6A4 (SERT)") |>
      stringr::str_replace("LOC110470479", "LOC110470479 (HTR3A)"),
    color = dplyr::if_else(m.value > 0, "BF", "WRM")
  )

pser = ggplot(serotonin_genes_expression) +
  aes(x = m.value, y = symbol) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(color = estimatedDEG, size = -log10(q.value))) +
  geom_point(
    data = serotonin_genes_expression |> dplyr::filter(estimatedDEG != "N"), 
    shape = "*",
    size = 5
    ) +
  scale_color_manual(values = c("BF" = colBF, "N" = "#CCCCCC", "WRM" = colWRM)) +
  labs(
    x = expression(paste({log[2]}, " Fold Change", sep = "")),
    size = expression(paste({-log[10]}, " FDR", sep = "")),
  ) +
  theme_bw(base_size = BASESIZE) +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = BASESIZE),
    axis.text.y = element_text(size = BASESIZE, face = "italic"),
    legend.position = c(-.7, 0),
    legend.justification = c(0, 0),
    legend.text = element_text(size = BASESIZE),
    panel.grid.minor = element_blank()
  ) +
  guides(color = "none")
pser


## cowplot ---------------------------------------------------------------------

d = cowplot::plot_grid(exp_lhx8, exp_cfap54, scale = .95)
ab = cowplot::plot_grid(pvol, pser, ncol = 2, scale = .99, rel_widths = c(3, 2), labels = c("a", "b"), label_size = LABELSIZE)
cd = cowplot::plot_grid(pven, d, ncol = 2, scale = c(.95, 1), rel_widths = c(3, 2), labels = c("c", "d"), label_size = LABELSIZE)
abcd = cowplot::plot_grid(ab, cd, nrow = 2, rel_heights = c(2, 1))
ggsave("images/figure3.png", abcd, w = 183, h = 160, units = "mm", bg = "#FFFFFF")

