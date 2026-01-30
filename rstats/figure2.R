library(conflicted)
library(tidyverse)
source("./rstats/common_settings.R")
source("./rstats/annot.R")


acc2chr= readr::read_tsv("out/sequence_report.tsv") |>
  dplyr::rename(CHROM = `RefSeq seq accession`, chr = `Chromosome name`) |>
  dplyr::select(CHROM, chr)

chr_levels = c("1", "1A", "2", "3", "4", "4A", as.character(seq(5,15)), 
               as.character(seq(17,29)), "Z", "MT")
chr_labels = c(1, "1A", 2, 3, 4, "4A", 5:10, 12, 14, 17, 20, 29, "Z")
chr_colors = c("#666666", "#BBBBBB")


## FST -------------------------------------------------------------------------

fst = readr::read_tsv("out/fst/bf_vs_wrm.windowed.weir.fst") |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_relevel(chr, chr_levels)) |>
  dplyr::filter(!chr %in% c("Un", "MT")) |>
  dplyr::group_by(chr) |>
  dplyr::mutate(Z_FST = (WEIGHTED_FST - mean(WEIGHTED_FST))/sd(WEIGHTED_FST)) |>
  dplyr::ungroup()

.fst_bp_cum = fst |>
  dplyr::group_by(chr) |>
  dplyr::summarise(max_bp = max(BIN_START)) |>
  dplyr::mutate(bp_add = dplyr::lag(cumsum(max_bp), default = 0)) |>
  dplyr::select(chr, bp_add)

fst = fst |>
  dplyr::inner_join(.fst_bp_cum, by = "chr") |>
  dplyr::mutate(bp4axis = BIN_START + bp_add)

df4axis = fst |>
  dplyr::group_by(chr) |>
  dplyr::summarize(center = mean(bp4axis)) |>
  dplyr::mutate(chr = dplyr::if_else(chr %in% chr_labels, as.character(chr), ""))

top10peaks = fst |>
  dplyr::select(CHROM, chr, BIN_START, BIN_END, Z_FST, bp4axis) |>
  dplyr::slice_max(order_by = Z_FST, prop = .001) |>
  annot(seqnames = "CHROM", start = "BIN_START", end = "BIN_END") |>
  tidyr::drop_na(gene) |>
  dplyr::group_by(gene) |>
  dplyr::slice_max(order_by = Z_FST, n = 1, with_ties = FALSE) |>
  dplyr::filter(stringr::str_detect(gene, "^LOC", negate = TRUE)) |>
  dplyr::group_by(chr) |>
  dplyr::slice_max(order_by = Z_FST, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::slice_max(order_by = Z_FST, n = 10) |>
  dplyr::filter(gene != "ANGPTL2")  # remove from visualization to avoid label ovorlap

goi = c("CFAP54", "OCA2", "FOXP2", "SOD2", "CSMD1")
gene4fst = fst |>
  dplyr::select(CHROM, chr, BIN_START, BIN_END, Z_FST, bp4axis) |>
  dplyr::slice_max(order_by = Z_FST, prop = .001) |>
  annot(seqnames = "CHROM", start = "BIN_START", end = "BIN_END") |>
  tidyr::drop_na(gene) |>
  dplyr::group_by(gene) |>
  dplyr::slice_max(order_by = Z_FST, n = 1, with_ties = FALSE) |>
  dplyr::filter(gene %in% goi)


pfst = ggplot(fst) +
  aes(x = bp4axis, y = Z_FST, color = forcats::as_factor(chr)) +
  geom_point() +  # size  =.5
  geom_text(
    data = top10peaks,
    aes(y = Z_FST + .5, label = gene), 
    fontface = "italic", 
    color = "#333333", 
    size = BASESIZE - 1, 
    size.unit = "pt"
  ) +
  ggrepel::geom_text_repel(
    data = gene4fst, 
    aes(y = Z_FST, label = gene), 
    fontface = "italic", 
    color = "#333333", 
    size = (BASESIZE-1)/2.83465,
    segment.size = .25,
    min.segment.length = 0,
    #nudge_y = 10 - gene4fst$Z_FST, 
    #direction = "x",
    #angle = 90,
    #hjust = 0,
    ylim = c(5.5, NA),
    max.overlaps = Inf
  ) +  
  geom_hline(yintercept = quantile(fst$Z_FST, .999), linetype = "longdash") +
  scale_x_continuous(expand = c(.01, .00), label = df4axis$chr, breaks = df4axis$center) +
  scale_y_continuous(
    expand = c(.02, .02), 
    limit = c(0, 9),
    sec.axis = dup_axis(breaks = quantile(fst$Z_FST, .999), labels = c("99.9%"))
  ) +
  scale_color_manual(values = rep(chr_colors, unique(length(df4axis$chr)))) +
  labs(y = expression(paste("Z", italic(F)[ST]))) +
  theme_classic(base_size = BASESIZE) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y.right = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = BASESIZE),
    axis.ticks.x = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.title.y.right = element_blank()
  )


## PI-ratio --------------------------------------------------------------------

bf_pi = readr::read_tsv("out/pi/bf.windowed.pi") |>
  dplyr::rename(BF = PI) |>
  dplyr::select(!N_VARIANTS)
wrm_pi = readr::read_tsv("out/pi/wrm.windowed.pi") |>
  dplyr::rename(WRM = PI) |>
  dplyr::select(!N_VARIANTS)

pi_ratio = bf_pi |>
  dplyr::inner_join(wrm_pi, by = c("CHROM", "BIN_START", "BIN_END")) |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_inorder(chr)) |>
  dplyr::mutate(log_pi_ratio = -log10(BF/WRM)) |>
  dplyr::filter(!chr %in% c("Un", "MT"))

.pi_ratio_bp_cum = pi_ratio |>
  dplyr::group_by(chr) |>
  dplyr::summarise(max_bp = max(BIN_START)) |>
  dplyr::mutate(bp_add = dplyr::lag(cumsum(max_bp), default = 0)) |>
  dplyr::select(chr, bp_add)

pi_ratio = pi_ratio |>
  dplyr::inner_join(.pi_ratio_bp_cum, by = "chr") |>
  dplyr::mutate(bp4axis = BIN_START + bp_add)

df4axis = pi_ratio |>
  dplyr::group_by(chr) |>
  dplyr::summarize(center = mean(bp4axis)) |>
  dplyr::mutate(chr = dplyr::if_else(chr %in% chr_labels, as.character(chr), ""))

gene4pi_bf = pi_ratio |>
  dplyr::select(CHROM, chr, BIN_START, BIN_END, log_pi_ratio, bp4axis) |>
  dplyr::filter(log_pi_ratio > quantile(pi_ratio$log_pi_ratio, .9995)) |>
  #dplyr::filter(log_pi_ratio > 1) |>
  annot(seqnames = "CHROM", start = "BIN_START", end = "BIN_END") |>
  tidyr::drop_na(gene) |>
  dplyr::group_by(gene) |>
  dplyr::slice_max(order_by = log_pi_ratio, n = 1, with_ties = FALSE) |>
  dplyr::group_by(chr) |>  #
  dplyr::slice_max(order_by = log_pi_ratio, n = 1, with_ties = FALSE) |>
  dplyr::filter(!chr %in% c(3, 15, 18, 22, 27))

gene4pi_wrm = pi_ratio |>
  dplyr::select(CHROM, chr, BIN_START, BIN_END, log_pi_ratio, bp4axis) |>
  dplyr::filter(log_pi_ratio < quantile(pi_ratio$log_pi_ratio, .0005)) |>
  #dplyr::filter(log_pi_ratio < -1.2) |>
  annot(seqnames = "CHROM", start = "BIN_START", end = "BIN_END") |>
  tidyr::drop_na(gene) |>
  dplyr::group_by(gene) |>
  dplyr::slice_min(order_by = log_pi_ratio, n = 1, with_ties = FALSE) |>
  dplyr::group_by(chr) |>  #
  dplyr::slice_min(order_by = log_pi_ratio, n = 1, with_ties = FALSE) |>
  dplyr::filter(!chr %in% c(3, 18, 24, 25))

ppi = ggplot(pi_ratio) +
  aes(x = bp4axis, y = log_pi_ratio, color = forcats::as_factor(chr)) +
  geom_point(size = .5) +
  geom_text(data = gene4pi_bf, aes(y = log_pi_ratio+.2, label = gene), fontface="italic", color = "#333333", size = BASESIZE-1, size.unit = "pt") +
  geom_text(data = gene4pi_wrm, aes(y = log_pi_ratio-.2, label = gene), fontface="italic", color = "#333333", size = BASESIZE-1, size.unit = "pt") +
  geom_hline(yintercept = quantile(pi_ratio$log_pi_ratio, .9995), linetype = "longdash") +
  geom_hline(yintercept = quantile(pi_ratio$log_pi_ratio, .0005), linetype = "longdash") +
  scale_x_continuous(expand = c(.02, .02), label = df4axis$chr, breaks = df4axis$center) +
  scale_y_continuous(
    expand = c(.05, .05),
    sec.axis = dup_axis(breaks = quantile(pi_ratio$log_pi_ratio, c(.0005, .9995)), labels = c("0.05%", "99.95%"))
  ) +
  scale_color_manual(values = rep(chr_colors, unique(length(df4axis$chr)))) +
  labs(
    x = "Chromosome",
    y = expression(paste(-log[10], " (", pi[BF], "/", pi[WRM], ")"))
  ) +
  theme_classic(base_size = BASESIZE) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y.right = element_blank(),
    axis.text = element_text(size = BASESIZE),
    axis.ticks.y.right = element_blank(),
    axis.title.y.right = element_blank()
  )

ggsave("images/manhattan_plot_pi_ratio.png", ppi, w = 183, h = 50, units = "mm", bg = "#FFFFFF")


## PI distribution of FST top 0.1% regions -------------------------------------

.fst = fst |> dplyr::select(CHROM, chr, BIN_START, BIN_END, Z_FST)
overlap_top01 = pi_ratio |> 
  dplyr::select(CHROM, chr, BIN_START, BIN_END, log_pi_ratio) |>
  dplyr::inner_join(.fst, by = c("CHROM", "chr", "BIN_START", "BIN_END")) |>
  dplyr::filter(Z_FST > quantile(fst$Z_FST, .999))

.fst_s = fst |> 
  dplyr::select(CHROM, chr, BIN_START, BIN_END, Z_FST) |>
  dplyr::filter(Z_FST > quantile(fst$Z_FST, .999))
.pi_bf_s = pi_ratio |> 
  dplyr::select(CHROM, chr, BIN_START, BIN_END, log_pi_ratio) |>
  dplyr::filter(log_pi_ratio > quantile(pi_ratio$log_pi_ratio, .9995))
.pi_wrm_s = pi_ratio |> 
  dplyr::select(CHROM, chr, BIN_START, BIN_END, log_pi_ratio) |>
  dplyr::filter(log_pi_ratio < quantile(pi_ratio$log_pi_ratio, .0005))

overlap_select = dplyr::bind_rows(.pi_bf_s, .pi_wrm_s) |>
  dplyr::inner_join(.fst_s, by = c("CHROM", "chr", "BIN_START", "BIN_END")) |>
  dplyr::filter(Z_FST > 6) |>
  annot(seqnames = "CHROM", start = "BIN_START", end = "BIN_END")

overlap_select_label = overlap_select |>
  dplyr::group_by(chr, gene) |>
  dplyr::summarise(log_pi_ratio = mean(log_pi_ratio)) |>
  dplyr::mutate(label = stringr::str_c(gene, " (chr. ", chr, ")"))

p0 = ggplot(overlap_top01) +
  aes(x = log_pi_ratio) +
  geom_point(aes(y = 1), size = 1.5, color = "#888888", fill = "#888888", shape = 21, alpha = .2) +
  geom_point(data = overlap_select, size = 1.5, shape = 21, aes(y = 1, fill = log_pi_ratio)) +
  geom_text(data = overlap_select_label, aes(y = 1, label = label), vjust = -1.5, size = BASESIZE - 1, size.unit = "pt") +
  geom_vline(xintercept = quantile(pi_ratio$log_pi_ratio, c(.0005, .9995)), linetype = "dashed") +
  scale_fill_gradient(high = colBF, low = colWRM) +
  labs(
    x = expression(paste(-log[10], " (", pi[BF], "/", pi[WRM], ")")),
    title = expression(paste(pi, "-ratio distribution of Z", italic(F)[ST], " top 0.1% regions"))
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = BASESIZE) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = BASESIZE),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(t = 15, r = 30, b = 5, l = 30)
  )

pi = cowplot::ggdraw() +
  cowplot::draw_plot(p0) +
  cowplot::draw_image("photos/wrm.nobg.png", scale = .25, x = -.45, y = -.2) +
  cowplot::draw_image("photos/bf.nobg.png", scale = .25, x = .45, y = -.2)
pi


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
  dplyr::mutate(chr = forcats::fct_inorder(chr), BIN_START = BIN_START + 1, BIN_END = BIN_START + 10000) |>
  dplyr::filter(!chr %in% c("Un", "MT"))


## Chr. 4A ---------------------------------------------------------------------

.chr4a_start = 19500000
.chr4a_end = 19860000

fst_chr4a = fst |>
  dplyr::filter(chr == "4A" & BIN_START >= .chr4a_start & BIN_START <= .chr4a_end) |>
  ggplot() +
  aes(x = (BIN_START + BIN_END)/2, y = Z_FST) +
  #geom_rect(aes(xmin = 19710001, xmax = 19725000, ymin = -Inf, ymax = Inf), fill = "#eeeeee") +
  geom_line() +
  scale_x_continuous(limits = c(.chr4a_start, .chr4a_end), expand = expansion(mult = c(0, 0))) +
  labs(title = "Chromosome 4A") +
  theme_classic(base_size = BASESIZE) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = BASESIZE),
    axis.ticks.x = element_blank()
  )
fst_chr4a

pi_chr4a = pi_ratio |>
  dplyr::filter(chr == "4A" & BIN_START >= .chr4a_start & BIN_START <= .chr4a_end) |>
  tidyr::pivot_longer(c("BF", "WRM"), names_to = "species", values_to = "PI") |>
  ggplot() +
  aes(x = (BIN_START + BIN_END)/2, y = log10(PI)) +
  #geom_rect(aes(xmin = 19710001, xmax = 19725000, ymin = -Inf, ymax = Inf), fill = "#eeeeee") +
  geom_line(aes(color = species), linewidth = .5) +
  scale_x_continuous(limits = c(.chr4a_start, .chr4a_end), expand = expansion(mult = c(0, 0))) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = BASESIZE) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = BASESIZE),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.justification = c(1, .5),
    legend.background = element_rect(fill = "NA"),
    legend.key = element_rect(fill = "NA")
  ) +
  guides(color = guide_legend(nrow = 2))
pi_chr4a

tD_chr4a = tD |>
  dplyr::filter(chr == "4A" & BIN_START >= .chr4a_start & BIN_START <= .chr4a_end) |>
  tidyr::pivot_longer(c("BF", "WRM"), names_to = "species", values_to = "TajimasD") |>
  ggplot() +
  aes(x = (BIN_START + BIN_END)/2, y = TajimasD) +
  geom_line(aes(color = species)) +
  labs(y = expression(paste("Tajima's ", italic(D)))) +
  scale_x_continuous(limits = c(.chr4a_start, .chr4a_end), expand = expansion(mult = c(0, 0))) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = BASESIZE) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = BASESIZE),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
tD_chr4a

genes4a = lonStrDom2_gene_list |>
  dplyr::left_join(acc2chr, by = dplyr::join_by(Accession == CHROM)) |>
  dplyr::filter(chr == "4A" & Begin >= .chr4a_start & Begin <= .chr4a_end) |>
  dplyr::filter(stringr::str_detect(Symbol, "^LOC", negate = TRUE)) |>
  dplyr::mutate(vj = c(2.5, 2.5))

gene_chr4a = fst |>
  dplyr::filter(chr == "4A" & BIN_START >= .chr4a_start & BIN_START <= .chr4a_end) |>
  ggplot() +
  aes(x = (BIN_START + BIN_END)/2) +
  geom_segment(data = genes4a, aes(x = Begin, xend = End), y = .2, linewidth = 2, color = "#08519c") +
  geom_point(data = genes4a, aes(x = (Begin+End)/2, shape = Orientation), y = .2, size = 3, color = "#FFFFFF") +
  geom_text(data = genes4a, aes(x = (Begin+End)/2, label = Symbol, vjust = vj), y = .2, 
            size = BASESIZE-1, size.unit = "pt", fontface = "italic") +
  scale_x_continuous(
    labels = scales::label_number(scale = 1/1000000, suffix = "M"),
    limits = c(.chr4a_start, .chr4a_end),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(limits = c(.1, .3)) +
  scale_shape_manual(values = c("plus" = 62, "minus" = 60)) +
  theme_classic(base_size = BASESIZE) +
  theme(
    plot.title = element_text(size = BASESIZE),
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = BASESIZE),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
gene_chr4a

p4a = cowplot::plot_grid(fst_chr4a, pi_chr4a, gene_chr4a, ncol = 1, align = "v", axis = "lr")


## Chr. 8 ----------------------------------------------------------------------

.chr8_start = 30840000

fst_chr8 = fst |>
  dplyr::filter(chr == 8 & BIN_START >= .chr8_start) |>
  ggplot() +
  aes(x = (BIN_START + BIN_END)/2, y = Z_FST) +
  #geom_rect(aes(xmin = 30940001, xmax = 30955000, ymin = -Inf, ymax = Inf), fill = "#eeeeee") +
  geom_line() +
  scale_x_continuous(limits = c(.chr8_start, NA), expand = expansion(mult = c(0, .05))) +
  labs(
    y = expression(paste("Z", italic(F)[ST])),
    title = "Chromosome 8"
  ) +
  theme_classic(base_size = BASESIZE) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = BASESIZE),
    axis.ticks.x = element_blank()
  )
fst_chr8

pi_chr8 = pi_ratio |>
  dplyr::filter(chr == 8 & BIN_START >= .chr8_start) |>
  tidyr::pivot_longer(c("BF", "WRM"), names_to = "species", values_to = "PI") |>
  ggplot() +
  aes(x = (BIN_START + BIN_END)/2, y = log10(PI)) +
  #geom_rect(aes(xmin = 30940001, xmax = 30955000, ymin = -Inf, ymax = Inf), fill = "#eeeeee") +
  geom_line(aes(color = species), linewidth = .5) +
  labs(y = expression(paste(log[10], " ", pi))) +
  scale_x_continuous(limits = c(.chr8_start, NA), expand = expansion(mult = c(0, .05))) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = BASESIZE) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = BASESIZE),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
pi_chr8

tD_chr8 = tD |>
  dplyr::filter(chr == 8 & BIN_START >= .chr8_start) |>
  tidyr::pivot_longer(c("BF", "WRM"), names_to = "species", values_to = "TajimasD") |>
  ggplot() +
  aes(x = (BIN_START + BIN_END)/2, y = TajimasD) +
  geom_line(aes(color = species)) +
  #labs(y = expression(paste("Tajima's ", italic(D)))) +
  scale_x_continuous(limits = c(.chr8_start, NA), expand = expansion(mult = c(0, .05))) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = BASESIZE) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = BASESIZE),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
tD_chr8

genes8 = lonStrDom2_gene_list |>
  dplyr::left_join(acc2chr, by = dplyr::join_by(Accession == CHROM)) |>
  dplyr::filter(chr ==8 & Begin >= .chr8_start) |>
  dplyr::filter(stringr::str_detect(Symbol, "^LOC", negate = TRUE)) |>
  dplyr::mutate(vj = c(2.5, -1.5, 2.5, -1.5, 2.5))

gene_chr8 = fst |>
  dplyr::filter(chr == 8 & BIN_START >= .chr8_start) |>
  ggplot() +
  aes(x = (BIN_START + BIN_END)/2) +
  geom_segment(data = genes8, aes(x = Begin, xend = End), y = .2, linewidth = 2, color = "#08519c") +
  geom_point(data = genes8, aes(x = (Begin+End)/2, shape = Orientation), y = .2, size = 3, color = "#FFFFFF") +
  geom_text(data = genes8, aes(x = (Begin+End)/2, label = Symbol, vjust = vj), y = .2, 
            size = BASESIZE-1, size.unit = "pt", fontface = "italic") +
  scale_x_continuous(
    labels = scales::label_number(scale = 1/1000000, suffix = "M"),
    limits = c(.chr8_start, NA),
    expand = expansion(mult = c(0, .05))
  ) +
  scale_y_continuous(limits = c(.1, .3)) +
  scale_shape_manual(values = c("plus" = 62, "minus" = 60)) +
  labs(y = "Genes") +
  theme_classic(base_size = BASESIZE) +
  theme(
    plot.title = element_text(size = BASESIZE),
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = BASESIZE),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
gene_chr8

p8 = cowplot::plot_grid(fst_chr8, pi_chr8, gene_chr8, ncol = 1, align = "v", axis = "lr")


## Align plots -----------------------------------------------------------------

cd = cowplot::plot_grid(p8, p4a, ncol = 2, labels = c("c", "d"), label_size = LABELSIZE, scale = .95, rel_widths = c(4, 5))
abcd = cowplot::plot_grid(pfst, pi, cd, labels = c("a", "b", ""), label_size = LABELSIZE, nrow = 3, rel_heights = c(5, 4, 6))
p = cowplot::ggdraw() +
  cowplot::draw_plot(abcd) +
  cowplot::draw_line(x = c(.09, .31, .31), y = c(.52, .38, .15), linetype = "dashed", linewidth = .2) +
  cowplot::draw_line(x = c(.15, .35, .35), y = c(.52, .38, .15), linetype = "dashed", linewidth = .2) +
  cowplot::draw_line(x = c(.86, .72, .72), y = c(.52, .38, .15), linetype = "dashed", linewidth = .2) +
  cowplot::draw_line(x = c(.91, .74, .74), y = c(.52, .38, .15), linetype = "dashed", linewidth = .2)
ggsave("images/figure2.png", p, w = 183, h = 150, units = "mm", bg = "#FFFFFF")

