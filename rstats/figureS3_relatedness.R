library(conflicted)
library(tidyverse)


indv2sex = readr::read_tsv("out/sample_information.tsv") |>
  dplyr::select(sample, sex_from_depth)

### WRM

relatedness_wrm = readr::read_tsv("out/relatedness/wrm.relatedness") |>
  dplyr::mutate(INDV1 = stringr::str_split_i(INDV1, "/", 2)) |>
  dplyr::mutate(INDV2 = stringr::str_split_i(INDV2, "/", 2)) |>
  dplyr::left_join(indv2sex, by = dplyr::join_by(INDV1 == sample)) |>
  dplyr::rename(INDV1_sex = sex_from_depth) |>
  dplyr::left_join(indv2sex, by = dplyr::join_by(INDV2 == sample)) |>
  dplyr::rename(INDV2_sex = sex_from_depth) |>
  dplyr::filter(INDV1 != INDV2) |>
  dplyr::mutate(sex_comb = dplyr::case_when(
    INDV1_sex == "M" & INDV2_sex == "M" ~ "Male-Male",
    INDV1_sex == "F" & INDV2_sex == "F" ~ "Female-Female",
    .default = "Female-Male"
  )) |>
  dplyr::mutate(group = "WRM")

kwt = kruskal.test(
  RELATEDNESS_AJK ~ sex_comb,
  data = relatedness_wrm
) |> print()

pwt = pairwise.wilcox.test(
  x = relatedness_wrm$RELATEDNESS_AJK,
  g = relatedness_wrm$sex_comb,
  p.adjust.method = "BH",
  paired = FALSE
) |> print()

p = ggplot(relatedness_wrm) +
  aes(x = sex_comb, y = 4) +
  geom_point(color = NA) +
  annotate("segment", x = 1, xend = 2, y = 1, yend = 1) +
  annotate("segment", x = 1, xend = 3, y = 2, yend = 2) +
  annotate("segment", x = 2, xend = 3, y = 3, yend = 3) +
  annotate("text", size = 3, x = 1.5, y = 1.5, label = "italic(P) == 4.85 %*% 10^-11", parse = TRUE) +
  annotate("text", size = 3, x = 2, y = 2.5, label = "italic(P) == 1.01 %*% 10^-7", parse = TRUE) +
  annotate("text", size = 3, x = 2.5, y = 3.5, label = "italic(P) < 2.0 %*% 10^-16", parse = TRUE) +
  theme_void()

r = ggplot(relatedness_wrm) +
  aes(x = sex_comb, y = RELATEDNESS_AJK) +
  geom_violin(color = NA, fill = "#CCCCCC") +
  geom_boxplot(width = .2, outliers = FALSE) +
  labs(
    x = "WRM",
    y = expression(paste("Relatedness (", A[italic(jk)], ")"))
  ) +
  theme_classic()

gwrm = cowplot::plot_grid(p, r, nrow = 2, rel_heights = c(1, 4), align = "v", axis = "lr")

### BF

relatedness_bf = readr::read_tsv("out/relatedness/bf.relatedness") |>
  dplyr::mutate(INDV1 = stringr::str_split_i(INDV1, "/", 2)) |>
  dplyr::mutate(INDV2 = stringr::str_split_i(INDV2, "/", 2)) |>
  dplyr::left_join(indv2sex, by = dplyr::join_by(INDV1 == sample)) |>
  dplyr::rename(INDV1_sex = sex_from_depth) |>
  dplyr::left_join(indv2sex, by = dplyr::join_by(INDV2 == sample)) |>
  dplyr::rename(INDV2_sex = sex_from_depth) |>
  dplyr::filter(INDV1 != INDV2) |>
  dplyr::mutate(sex_comb = dplyr::case_when(
    INDV1_sex == "M" & INDV2_sex == "M" ~ "Male-Male",
    INDV1_sex == "F" & INDV2_sex == "F" ~ "Female-Female",
    .default = "Female-Male"
  )) |>
  dplyr::mutate(group = "BF")

kwt = kruskal.test(
  RELATEDNESS_AJK ~ sex_comb,
  data = relatedness_bf
) |> print()

pwt = pairwise.wilcox.test(
  x = relatedness_bf$RELATEDNESS_AJK,
  g = relatedness_bf$sex_comb,
  p.adjust.method = "BH",
  paired = FALSE
) |> print()

p = ggplot(relatedness_bf) +
  aes(x = sex_comb, y = 4) +
  geom_point(color = NA) +
  annotate("segment", x = 1, xend = 2, y = 1, yend = 1) +
  annotate("segment", x = 1, xend = 3, y = 2, yend = 2) +
  annotate("segment", x = 2, xend = 3, y = 3, yend = 3) +
  annotate("text", size = 3, x = 1.5, y = 1.5, label = "italic(P) == 1.06 %*% 10^-6", parse = TRUE) +
  annotate("text", size = 3, x = 2, y = 2.5, label = "italic(P) == 0.83", parse = TRUE) +
  annotate("text", size = 3, x = 2.5, y = 3.5, label = "italic(P) < 2.0 %*% 10^-16", parse = TRUE) +
  theme_void()

r = ggplot(relatedness_bf) +
  aes(x = sex_comb, y = RELATEDNESS_AJK) +
  geom_violin(color = NA, fill = "#CCCCCC") +
  geom_boxplot(width = .2, outliers = FALSE) +
  labs(
    x = "BF",
    y = expression(paste("Relatedness (", A[italic(jk)], ")"))
  ) +
  theme_classic()

gbf = cowplot::plot_grid(p, r, nrow = 2, rel_heights = c(1, 4), align = "v", axis = "lr")

cowplot::plot_grid(gwrm, gbf, scale = .95)
ggsave("images/relatedness.png", w = 9, h = 4, bg = "#FFFFFF")
