library(conflicted)
library(tidyverse)
source("./rstats/common_settings.R")


acc2chr= readr::read_tsv("out/sequence_report.tsv") |>
  dplyr::rename(CHROM = `RefSeq seq accession`, chr = `Chromosome name`) |>
  dplyr::select(CHROM, chr)

chr_levels = c("1", "1A", "2", "3", "4", "4A", as.character(seq(5,15)), 
               as.character(seq(17,29)), "Z", "MT")
chr_labels = c(1, "1A", 2, 3, 4, "4A", 5:10, 12, 14, 17, 20, 29, "Z")
chr_colors = c("#666666", "#BBBBBB")

fst = readr::read_tsv("out/fst/bf_vs_wrm.windowed.weir.fst") |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::mutate(chr = forcats::fct_relevel(chr, chr_levels)) |>
  dplyr::filter(!chr %in% c("Un", "MT")) |>
  dplyr::group_by(chr) |>
  dplyr::mutate(Z_FST = (WEIGHTED_FST - mean(WEIGHTED_FST))/sd(WEIGHTED_FST)) |>
  dplyr::ungroup()

zfst001 = readr::read_tsv("out/zfst001.genes.txt", skip = 1) |>
  tidyr::pivot_longer(dplyr::starts_with("variants_"), names_to = "variants_type", values_to = "count") |>
  #dplyr::filter(stringr::str_detect(variants_type, "^variants_impact")) |>
  dplyr::filter(stringr::str_detect(variants_type, "^variants_effect")) |>
  dplyr::group_by(variants_type) |>
  dplyr::summarise(sum(count))

zfst001vcf = readr::read_tsv("out/zfst001.maf005.annot.vcf", comment = "##") |>
  dplyr::rename(CHROM = `#CHROM`)
daf = zfst001vcf |>
  dplyr::select(CHROM, POS, dplyr::starts_with("bam")) |>
  tidyr::pivot_longer(dplyr::starts_with("bam"), names_to = "sample", values_to = "GT") |>
  dplyr::mutate(sample = stringr::str_split(sample, "/", simplify = TRUE)[,2]) |>
  dplyr::mutate(group = dplyr::if_else(stringr::str_detect(sample, "WRM"), "WRM", "BF")) |>
  dplyr::mutate(GT = stringr::str_split(GT, ":", simplify = TRUE)[,1]) |>
  dplyr::mutate(AF = dplyr::case_when(GT == "1/1" ~ 1, GT == "0/1" ~ .5, GT == "0/0" ~ 0, .default = NA)) |>
  dplyr::group_by(CHROM, POS, group) |>
  dplyr::summarise(meanAF = mean(AF, na.rm = TRUE)) |>
  tidyr::pivot_wider(names_from = group, values_from = meanAF) |>
  dplyr::mutate(dAF = abs(BF - WRM))
annot = zfst001vcf |>
  dplyr::select(CHROM, POS, INFO) |>
  dplyr::mutate(INFO = stringr::str_split(INFO, "ANN=", simplify = TRUE)[,2]) |>
  tidyr::separate_rows(INFO, sep = ",") |>
  tidyr::separate(INFO, sep = "\\|", into = ann_field, convert = TRUE) |>
  dplyr::right_join(daf, by = c("CHROM", "POS"))

## chr. 8 ----------------------------------------------------------------------

.chr8_start = 30840000

annot_chr8 = annot |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::filter(chr == 8 & POS >= .chr8_start)

fst_chr8 = fst |>
  dplyr::filter(chr == 8 & BIN_START >= .chr8_start) |>
  ggplot() +
  aes(x = (BIN_START + BIN_END)/2, y = Z_FST) +
  geom_hline(yintercept = quantile(fst$Z_FST, .999), linetype = "dashed") +
  geom_line() +
  geom_point(data = annot_chr8, aes(x = POS, color = dAF, y = 1), shape = "|", size = 2) +
  scale_x_continuous(
    labels = scales::label_number(scale = 1/1000000, suffix = "M"),
    limits = c(.chr8_start, NA),
    expand = expansion(mult = c(0, .05))
  ) +
  scale_color_viridis_c(option = "cividis") +
  labs(
    y = expression(paste("Z", italic(F)[ST])),
    title = "Chromosome 8"
  ) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
fst_chr8

cols = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attribute")
gff_chr8 = readr::read_tsv("out/lonStrDom2_chr8_genomic.gff", comment = "##", col_names = cols) |>
  dplyr::mutate(gene = stringr::str_extract(attribute, "(?<=gene=)[^;]+")) |>
  dplyr::mutate(gene = dplyr::if_else(stringr::str_detect(gene, "^LOC"), "", gene)) |>
  dplyr::filter(start >= .chr8_start)

genes_chr8 = gff_chr8 |> dplyr::filter(type == "gene")
exons_chr8 = gff_chr8 |> dplyr::filter(type == "exon")

gene_chr8 = ggplot(genes_chr8) +
  aes(y = 1) +
  geom_segment(aes(x = start, xend = end, y = 1, yend = 1), linewidth = 1, color = "#9B9477") +
  geom_segment(data = exons_chr8, aes(x = start, xend = end, y = 1, yend = 1), linewidth = 5, color = "#9B9477") +
  geom_text(aes(x = (start + end)/2, label = gene, vjust = c(2, -1, 2, -1, 2, 2)), fontface = "italic", size = 2.5) +
  scale_x_continuous(
    labels = scales::label_number(scale = 1/1000000, suffix = "M"),
    limits = c(.chr8_start, NA),
    expand = expansion(mult = c(0, .05))
  ) +
  scale_color_viridis_c(option = "cividis") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_chr8 = cowplot::plot_grid(fst_chr8, gene_chr8, nrow = 2, align = "v", axis = "lr", rel_heights = c(2, 1))
p_chr8

## chr. 4A ---------------------------------------------------------------------

.chr4a_start = 19500000
.chr4a_end = 19860000

annot_chr4a = annot |>
  dplyr::left_join(acc2chr, by = "CHROM") |>
  dplyr::filter(chr == "4A" & POS >= .chr4a_start & POS <= .chr4a_end)

fst_chr4a = fst |>
  dplyr::filter(chr == "4A" & BIN_START >= .chr4a_start & BIN_START <= .chr4a_end) |>
  ggplot() +
  aes(x = (BIN_START + BIN_END)/2, y = Z_FST) +
  geom_hline(yintercept = quantile(fst$Z_FST, .999), linetype = "dashed") +
  geom_line() +
  geom_point(data = annot_chr4a, aes(x = POS, color = dAF, y = -2.5), shape = "|", size = 2) +
  scale_x_continuous(
    labels = scales::label_number(scale = 1/1000000, suffix = "M"),
    limits = c(.chr4a_start, .chr4a_end),
    expand = expansion(mult = c(0, .05))
  ) + 
  scale_color_viridis_c(option = "cividis") +
  labs(
    y = expression(paste("Z", italic(F)[ST])),
    title = "Chromosome 4A"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "none"
  )
fst_chr4a

gff_chr4a = readr::read_tsv("out/lonStrDom2_chr4a_genomic.gff", comment = "##", col_names = cols) |>
  dplyr::mutate(gene = stringr::str_extract(attribute, "(?<=gene=)[^;]+")) |>
  dplyr::mutate(gene = dplyr::if_else(stringr::str_detect(gene, "^LOC"), "", gene)) |>
  dplyr::filter(start >= .chr4a_start & end <= .chr4a_end)

genes_chr4a = gff_chr4a |> dplyr::filter(type == "gene")
exons_chr4a = gff_chr4a |> dplyr::filter(type == "exon")

gene_chr4a = ggplot(genes_chr4a) +
  aes(y = 1) +
  geom_segment(data = genes_chr4a, aes(x = start, xend = end, y = 1, yend = 1), linewidth = 1, color = "#9B9477") +
  geom_segment(data = exons_chr4a, aes(x = start, xend = end, y = 1, yend = 1), linewidth = 5, color = "#9B9477") +
  geom_text(aes(x = (start + end)/2, label = gene), fontface = "italic", vjust = 2, size = 2.5) +
  scale_x_continuous(
    labels = scales::label_number(scale = 1/1000000, suffix = "M"),
    limits = c(.chr4a_start, .chr4a_end),
    expand = expansion(mult = c(0, .05))
  ) +
  scale_color_viridis_c(option = "cividis") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_chr4a = cowplot::plot_grid(fst_chr4a, gene_chr4a, nrow = 2, align = "v", axis = "lr", rel_heights = c(2, 1))
p_chr4a

## Annotation ------------------------------------------------------------------

snpeff = dplyr::bind_rows(annot_chr4a, annot_chr8) |>
  dplyr::mutate(
    Annotation = stringr::str_replace_all(Annotation, "_", " ") |> 
      stringr::str_to_sentence() |>
      stringr::str_replace_all(" prime utr", "' UTR") |>
      stringr::str_remove(" variant$"),
    chr = stringr::str_c("chr. ", chr)
  ) |>
  ggplot() +
  aes(x = chr) +
  geom_bar(stat = "count", position = "fill", aes(fill = Annotation)) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_fill_viridis_d(option = "A") +
  theme_classic() +
  theme(
    axis.title = element_blank()
  )
snpeff

## Cowplot ---------------------------------------------------------------------

p = cowplot::plot_grid(p_chr4a, p_chr8, snpeff, nrow = 1, labels = c("a", "", "b"), scale = .95)
ggsave("images/snpeff.png", p, w = 12, h = 3, bg = "#FFFFFF")
