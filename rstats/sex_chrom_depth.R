library(conflicted)
library(tidyverse)


sample_info = readr::read_tsv("out/sample_information_excel.tsv") |>
  dplyr::rename(sample = sample_id) |>
  dplyr::mutate(sex = dplyr::if_else(is.na(sex), "Un", sex)) |>
  dplyr::mutate(group = dplyr::case_when(
    species == "Bengalese finch" ~ "BF",
    species == "White-rumped munia" ~ "WRM"
  ))

## Sex identification from read depth ------------------------------------------

sex_chrom_depth = readr::read_tsv("out/sex_chrom_depth.tsv") |>
  dplyr::inner_join(sample_info, by = "sample") |>
  dplyr::mutate(ratio = meandp_sex_chrom/meandp) |>
  dplyr::mutate(sex_from_depth = dplyr::if_else(ratio > .75, "M", "F")) |>
  dplyr::rename(sex_from_PCR = sex) |>
  dplyr::arrange(sample)

sex_chrom_depth |> dplyr::group_by(group, sex_from_depth) |> dplyr::summarise(sex = n())

count = ggplot(sex_chrom_depth) +
  aes(x = ratio, color = sex_from_depth) +
  #geom_density() +
  geom_point(aes(y = 0), shape = 21, size = 3) +
  scale_x_continuous(limit = c(0, NA), breaks = c(0, .5, 1), label = c(0, .5, 1)) +
  scale_color_manual(
    values = c("F" = "#d6604d", "M" = "#4393c3"),
    labels = c("F" = "Female: WRM, n = 17; BF, n = 18", "M" = "Male: WRM, n = 14; BF, n = 34")
    ) +
  labs(x = "Ratio of read depth (Z chromosome / Autosome)") +
  theme_test(base_size = 12) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

indv_bf = sex_chrom_depth |>
  dplyr::filter(group == "BF") |>
  tidyr::pivot_longer(dplyr::starts_with("mean"), names_to = "chrom", values_to = "depth") |>
  ggplot() +
  aes(x = chrom, y = depth, fill = chrom) +
  geom_col() +
  scale_fill_manual(
    values = c("#666666", "#7b329466"),
    labels = c("Autosomes", "Z")
  ) +
  labs(title = "Bengalese finch", y = "Average read depth", fill = "Chromosome") +
  facet_wrap(vars(sample), ncol = 9, scales = "free_y") +
  theme_test(base_size = 12) +
  theme(
    #plot.title = element_text(margin = margin(t = 40,b = -30)),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )

indv_wrm = sex_chrom_depth |>
  dplyr::filter(group == "WRM") |>
  tidyr::pivot_longer(dplyr::starts_with("mean"), names_to = "chrom", values_to = "depth") |>
  ggplot() +
  aes(x = chrom, y = depth, fill = chrom) +
  geom_col() +
  scale_fill_manual(
    values = c("#666666", "#7b329466"),
    labels = c("Autosomes", "Z")
  ) +
  labs(title = "White-rumped munia", y = "Average read depth", fill = "Chromosome") +
  facet_wrap(vars(sample), ncol = 9, scales = "free_y") +
  theme_test(base_size = 12) +
  theme(
    #plot.title = element_text(margin = margin(t = 40,b = -30)),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = "inside",
    legend.justification = c(1, 0)
  ) +
  guides(fill = guide_legend(nrow = 1))

p = cowplot::plot_grid(
  count, NULL, indv_bf, indv_wrm, 
  nrow = 4, rel_heights = c(1, .1, 3.6, 2.4), 
  align = "v", axis = "l", 
  labels = c("a", "", "b", ""), label_size = 16
)
ggsave("images/sex_chrom_depth.png", p, w = 11, h = 12, bg = "#FFFFFF")

sex_chrom_depth |>
  dplyr::select(sample, species, group, colony, sex_from_depth) |>
  readr::write_tsv("out/sample_information.tsv")
