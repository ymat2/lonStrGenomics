library(conflicted)
library(tidyverse)
library(cowplot)

source("./rstats/common_settings.R")


## Photos ----------------------------------------------------------------------

wrm_photo = cowplot::ggdraw() +
  cowplot::draw_image("photos/wrm.nobg.png", scale = .95, valign = .25) +
  #ggplot2::annotate(geom = "text", x = .2, y = .9, label = "●", size = .size, hjust = 1.5, color = "#5ab4ac") +
  ggplot2::annotate(geom = "text", x = .5, y = .9, label = "White-rumped munia\n(WRM; Wild)", size = BASESIZE, size.unit = "pt", hjust = .5) +
  theme_void()

bf_photo = cowplot::ggdraw() +
  cowplot::draw_image("photos/bf.nobg.png", scale = .95, valign = .25) +
  #ggplot2::annotate(geom = "text", x = .2, y = .9, label = "●", size = .size, hjust = 1.5, color = "#d8b365") +
  ggplot2::annotate(geom = "text", x = .5, y = .9, label = "Bengalese finch\n(BF; Domesticated)", size = BASESIZE, size.unit = "pt", hjust = .5) +
  theme_void()

photos = cowplot::plot_grid(wrm_photo, bf_photo, ncol = 2, labels = c("a", ""), label_size = LABELSIZE, scale = .95)


## PCA -------------------------------------------------------------------------

sample_info = readr::read_tsv("out/sample_information.tsv") |> dplyr::rename(sex = sex_from_depth)

eigenvec = readr::read_delim("out/pca/lonchura.snp.pca.eigenvec") |>
  dplyr::mutate(sample = stringr::str_split(IID, "/", simplify = TRUE)[,2]) |>
  dplyr::left_join(sample_info, by = "sample") |>
  dplyr::mutate(group = forcats::fct_relevel(group, names(colors)))
eigenval = readr::read_csv("out/pca/lonchura.snp.pca.eigenval", col_names = "V1")
df = data.frame(pc = 1:nrow(eigenval), eigenval/sum(eigenval)*100)

pca12 = ggplot(eigenvec) +
  aes(PC1, PC2, color = group) +
  geom_point(shape = 16, size = 2) +
  scale_color_manual(values = colors) +
  labs(
    x = paste0("PC1 (", round(df[1,2], digits = 2), "%)"),
    y = paste0("PC2 (", round(df[2,2], digits = 2), "%)")
  ) +
  theme_test(base_size = BASESIZE) +
  theme(
    legend.position = "inside",
    legend.title = element_blank(),
    legend.justification = c(.01, .99),
    axis.text = element_text(size = BASESIZE)
  )
pca12


## ADMIXTURE -------------------------------------------------------------------

fam = readr::read_tsv("out/admixture/lonchura.snp.fam", col_names = c("FID", "IID", "X1", "X2", "X3", "X4")) |> dplyr::select(IID)

### K=2 -----

q2 = readr::read_table("out/admixture/lonchura.snp.2.Q", col_names = c("Pop0", "Pop1")) |>
  dplyr::bind_cols(fam) |>
  dplyr::left_join(eigenvec, by = "IID") |>
  dplyr::select(sample, colony, Pop0, Pop1, group) |>
  dplyr::group_by(group, colony) |>
  dplyr::mutate(sample = forcats::fct_inorder(sample)) |>
  tidyr::pivot_longer(dplyr::starts_with("Pop"), names_to = "anc", values_to = "prop") |>
  dplyr::mutate(k = "K = 2")

bar = ggplot(q2) +
  aes(x = sample, y = prop, fill = anc) +
  geom_bar(stat = "identity", position = "fill", fill = "transparent") +
  geom_segment(x = 1, xend = 31, y = 1, yend = 1, linewidth = .2) +
  ggplot2::annotate(geom = "text", x = (1+31)/2, y = 0, label = "WRM (n = 31)", vjust = -1, size = BASESIZE, size.unit = "pt") +
  geom_segment(x = 32, xend = 83, y = 1, yend = 1, linewidth = .2) +
  ggplot2::annotate(geom = "text", x = (32+83)/2, y = 0, label = "BF (n = 52)", vjust = -1, size = BASESIZE, size.unit = "pt") +
  #scale_x_continuous(limits = c(1, 83)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  labs(y = "") +
  theme_minimal(base_size = BASESIZE) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_blank()
  )
bar

padm2 = ggplot(q2) +
  aes(x = sample, y = prop, fill = anc) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Pop0" = colBF, "Pop1" = colWRM)) +
  labs(y = "K = 2") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_minimal(base_size = BASESIZE) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text = element_blank()
  )
padm2

### K=3 -----

q3 = readr::read_table("out/admixture/lonchura.snp.3.Q", col_names = c("Pop0", "Pop1", "Pop2")) |>
  dplyr::bind_cols(fam) |>
  dplyr::left_join(eigenvec, by = "IID") |>
  dplyr::select(sample, colony, Pop0, Pop1, Pop2, group) |>
  dplyr::group_by(group, colony) |>
  dplyr::mutate(sample = forcats::fct_inorder(sample)) |>
  tidyr::pivot_longer(dplyr::starts_with("Pop"), names_to = "anc", values_to = "prop") |>
  dplyr::mutate(k = "K = 3")

padm3 = ggplot(q3) +
  aes(x = sample, y = prop, fill = anc) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Pop0" = colWRM, "Pop1" = colBF, "Pop2" = lighter(colBF))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(y = "K = 3") +
  theme_minimal(base_size = BASESIZE) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text = element_blank()
  )
padm3


## Nucleotide diversity --------------------------------------------------------

acc2chr = readr::read_tsv("out/sequence_report.tsv") |>
  dplyr::rename(CHROM = `RefSeq seq accession`, chr = `Chromosome name`) |>
  dplyr::select(CHROM, chr)

bf_pi = readr::read_tsv("out/pi/bf.windowed.pi") |>
  dplyr::rename(BF = PI) |>
  dplyr::select(!N_VARIANTS)
wrm_pi = readr::read_tsv("out/pi/wrm.windowed.pi") |>
  dplyr::rename(WRM = PI) |>
  dplyr::select(!N_VARIANTS)

wind_pi = bf_pi |>
  dplyr::inner_join(wrm_pi, by = c("CHROM", "BIN_START", "BIN_END")) |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_inorder(chr)) |>
  dplyr::filter(!chr %in% c("W", "Z", "Un", "MT"))

wilcox_bf_vs_wrm = wilcox.test(wind_pi$BF, wind_pi$WRM, paired = TRUE)
wilcox_bf_vs_wrm

violin_pi = wind_pi |>
  tidyr::pivot_longer(c("BF", "WRM"), names_to = "pop", values_to = "pi") |>
  dplyr::mutate(pop = forcats::fct_relevel(pop, names(colors))) |>
  ggplot() +
  aes(x = pop, y = pi) +
  geom_violin(aes(fill = pop), linewidth = 0, alpha = .75) +
  geom_boxplot(
    width = .2,
    outlier.alpha = .5,
    outlier.shape = 16,
    outlier.size = 1
    ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = .03,
    label = as.character(expression(paste(italic(P), "< 2.2 × ", 10^{-16}))),
    parse = TRUE,
    size = BASESIZE,
    size.unit = "pt"
    ) +
  ggplot2::annotate("text", x = 1.5, y = .003, label = ">", size = 5) +
  scale_fill_manual(values = colors) +
  ylab(expression(paste("Nucleotide diversity (", {pi}, ")" , sep = ""))) +
  theme_test(base_size = BASESIZE) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text = element_text(size = BASESIZE)
  )
violin_pi


## Tajima's D ------------------------------------------------------------------

tajimasD_bf = readr::read_tsv("out/tajimasD/bf.Tajima.D") |>
  dplyr::rename(BF = TajimaD) |>
  dplyr::select(!N_SNPS)

tajimasD_wrm = readr::read_tsv("out/tajimasD/wrm.Tajima.D") |>
  dplyr::rename(WRM = TajimaD) |>
  dplyr::select(!N_SNPS)

wind_D = tajimasD_bf |>
  dplyr::inner_join(tajimasD_wrm, by = c("CHROM", "BIN_START")) |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_inorder(chr)) |>
  dplyr::filter(!chr %in% c("W", "Z", "Un", "MT")) |>
  dplyr::filter(!is.na(BF) & !is.na(WRM))

wind_D |> dplyr::summarise(BF = mean(BF), WRM = mean(WRM))
wilcox_bf_vs_wrm_D = wilcox.test(wind_D$BF, wind_D$WRM, paired = TRUE)
wilcox_bf_vs_wrm_D

violin_D = wind_D |>
  tidyr::pivot_longer(c("BF", "WRM"), names_to = "pop", values_to = "tajimaD") |>
  dplyr::mutate(pop = forcats::fct_relevel(pop, names(colors))) |>
  ggplot() +
  aes(x = pop, y = tajimaD) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin(aes(fill = pop), linewidth = 0, alpha = .75) +
  geom_boxplot(
    width = .2,
    outlier.alpha = .5,
    outlier.shape = 16,
    outlier.size = 1
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = 4,
    label = as.character(expression(paste(italic(P), "< 2.2 × ", 10^{-16}))),
    parse = TRUE,
    size = BASESIZE,
    size.unit = "pt"
  ) +
  scale_fill_manual(values = colors) +
  ylab(expression(paste("Tajima's ", italic(D) , sep = ""))) +
  theme_test(base_size = BASESIZE) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text = element_text(size = BASESIZE)
  )
violin_D


## Heterozygosity --------------------------------------------------------------

het = readr::read_tsv("out/hetero.het") |>
  dplyr::mutate(sample = stringr::str_split(INDV, "/", simplify = TRUE)[,2]) |>
  dplyr::inner_join(sample_info, by = "sample") |>
  dplyr::mutate(H_o = (N_SITES - `O(HOM)`)/N_SITES) |>
  dplyr::filter(sex != "Un")

glm_group_sex = glm(H_o ~ group * sex, data = het)
summary(glm_group_sex)

hetF = het |>
  dplyr::filter(sex == "F") |>
  dplyr::mutate(group = forcats::fct_relevel(group, names(colors))) |>
  ggplot() +
  aes(x = group, y = H_o) +
  geom_boxplot(
    aes(fill = group),
    outliers = FALSE
  ) +
  geom_jitter(
    aes(fill = group),
    color = "#FFFFFF",
    shape = 21,
    alpha = .9,
    size = 2,
    width = .2,
    height = 0
  ) +
  scale_fill_manual(values = colors) +
  labs(title = "Female", y = "Heterozygosity") +
  theme_classic(base_size = BASESIZE) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = .5),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title = element_text(size = BASESIZE),
    axis.text = element_text(size = BASESIZE)
  )

hetM = het |>
  dplyr::filter(sex == "M") |>
  dplyr::mutate(group = forcats::fct_relevel(group, names(colors))) |>
  ggplot() +
  aes(x = group, y = H_o) +
  geom_boxplot(
    aes(fill = group),
    outliers = FALSE
  ) +
  geom_jitter(
    aes(fill = group),
    color = "#FFFFFF",
    shape = 21,
    alpha = .9,
    size = 2,
    width = .2,
    height = 0
  ) +
  scale_fill_manual(values = colors) +
  labs(title = "Male") +
  theme_classic(base_size = BASESIZE) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = .5),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_line(linetype = "dashed"),
    axis.text.x = element_text(size = BASESIZE)
  )

phet = cowplot::plot_grid(hetF, hetM)
phet


## Align plots -----------------------------------------------------------------

# First align fig.c and fig.d by explicit call to align_plots().
cd = cowplot::align_plots(padm2, padm3, bar, violin_pi, align = "v", axis = "l")

padm = cowplot::plot_grid(cd[[1]], cd[[2]], cd[[3]], nrow = 3, rel_heights = c(2, 2, 1))
de = cowplot::plot_grid(cd[[4]], violin_D, ncol = 2, labels = c("d", "e"), label_size = LABELSIZE)
acde = cowplot::plot_grid(photos, padm, de, nrow = 3, rel_heights = c(3, 3, 4), labels = c("", "c", ""), label_size = LABELSIZE, align = "v", axis = "lr", scale = .97)
bf = cowplot::plot_grid(pca12, phet, nrow = 2, rel_heights = c(3, 2), labels = c("b", "f"), label_size = LABELSIZE, scale = .97)
p = cowplot::plot_grid(acde, bf, ncol = 2, rel_widths = c(3, 2))

ggsave("images/figure1.png", p, w = 183, h = 120, units = "mm", bg = "#FFFFFF")
