library(conflicted)
library(tidyverse)
source("./rstats/common_settings.R")


## PCA -----

sample_info = readr::read_tsv("out/sample_information.tsv") |> dplyr::rename(sex = sex_from_depth)

eigenvec_bf = readr::read_delim("out/pca/lonchura.snp.pca.eigenvec") |>
  dplyr::mutate(sample = stringr::str_split(IID, "/", simplify = TRUE)[,2]) |>
  dplyr::left_join(sample_info, by = "sample") |>
  dplyr::mutate(group = forcats::fct_relevel(group, names(colors))) |>
  dplyr::filter(group == "BF")
eigenval = readr::read_csv("out/pca/lonchura.snp.pca.eigenval", col_names = "V1")
df = data.frame(pc = 1:nrow(eigenval), eigenval/sum(eigenval)*100)

pca12_bf = eigenvec_bf |>
  dplyr::mutate(colony = forcats::fct_relevel(colony, c("Azabu University", "Okanoya Lab", "Tomakomai", "Others"))) |>
  ggplot() +
  aes(PC1, PC2, shape = colony) +
  geom_point(size = 2, alpha = 0.9, color = "#d8b365") +
  #scale_color_manual(values = colors) +
  labs(
    x = paste0("PC1 (", round(df[1,2], digits = 2), "%)"),
    y = paste0("PC2 (", round(df[2,2], digits = 2), "%)"),
    shape = "Colony"
  ) +
  theme_test(base_size = BASESIZE+1) +
  theme(
    legend.position = "inside",
    legend.justification = c(.01, .99),
    legend.background = element_blank()
  )
pca12_bf


## ADMIXTURE -----

### CV-eror -----

cve = readr::read_table("out/admixture/CV-error.txt", col_names = c("CV", "error", "K", "cv_error")) |>
  dplyr::select(K, cv_error) |>
  dplyr::mutate(K = stringr::str_extract(K, "[0-9]+") |> as.numeric()) |>
  ggplot() +
  aes(K, cv_error) +
  geom_point(color = "#333333") +
  geom_line(color = "#333333") +
  scale_x_continuous(labels = seq(1,8), breaks = seq(1,8)) +
  labs(y = "CV-error") +
  theme_test(base_size = BASESIZE+1)
cve


### K=3 -----

fam = readr::read_tsv("out/admixture/lonchura.snp.fam", col_names = c("FID", "IID", "X1", "X2", "X3", "X4")) |> dplyr::select(IID)

q3_bf = readr::read_table("out/admixture/lonchura.snp.3.Q", col_names = c("Pop0", "Pop1", "Pop2")) |>
  dplyr::bind_cols(fam) |>
  dplyr::inner_join(eigenvec_bf, by = "IID") |>
  dplyr::select(sample, colony, Pop0, Pop1, Pop2, group) |>
  dplyr::group_by(group, colony) |>
  dplyr::mutate(sample = forcats::fct_inorder(sample)) |>
  tidyr::pivot_longer(dplyr::starts_with("Pop"), names_to = "anc", values_to = "prop")

padm3_bf = q3_bf |>
  dplyr::ungroup() |>
  dplyr::mutate(colony = forcats::fct_relevel(colony, c("Tomakomai", "Azabu University", "Okanoya Lab", "Others"))) |>
  ggplot() +
  aes(y = prop, x = sample, fill = anc) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Pop0" = colWRM, "Pop1" = colBF, "Pop2" = lighter(colBF))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  facet_grid(cols = vars(colony), space = "free_x", scales = "free_x") +
  labs(x = "", title = "K = 3 (BF only)") +
  theme_minimal(base_size = BASESIZE+1) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = .5),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )
padm3_bf

p = cowplot::plot_grid(
  cowplot::plot_grid(
    pca12_bf,
    cve,
    labels = c("a", "b"),
    align = "h",
    axis = "tb",
    scale = .95
  ),
  padm3_bf,
  nrow = 2,
  rel_heights = c(3, 2),
  labels = c("", "c"),
  align = "v",
  axis = "lr",
  scale = c(1, .95)
)
ggsave("images/structure_result_bf.png", p, w = 7, h = 5, bg = "#FFFFFF")
